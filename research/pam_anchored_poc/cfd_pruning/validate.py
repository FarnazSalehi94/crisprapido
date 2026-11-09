#!/usr/bin/env python3
"""A: integrated downstream filtering; B: isolated native pruning; C: forward oracle.

Only bounded references. Native performance uses real FM-index traversal, not
model states. Existing production binaries/indexes are opened read-only.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import resource
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'research/pam_anchored_poc/cfd_feasibility'))
from study import fixtures, keep, load_tables, oracle, rc, representatives


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def command(cmd, timeout=60):
    r = subprocess.run(list(map(str,cmd)), cwd=ROOT, capture_output=True, text=True, timeout=timeout)
    assert r.returncode == 0, (cmd, r.stderr[-3000:])
    return r


def native_rows(text, expected):
    lines = text.splitlines()
    start = next(i for i,line in enumerate(lines) if line.startswith('guide\tcontig\t'))
    rows = list(csv.DictReader(lines[start:],delimiter='\t'))
    found = set()
    for row in rows:
        key = row['contig'],row['strand'],int(row['anchor'])
        assert key not in found, ('duplicate',key)
        found.add(key)
        hit = expected[key]
        assert (int(row['start']),int(row['end']),row['cigar'],int(row['mismatches'])) == (hit.start,hit.end,hit.cigar,hit.mm), (row,hit)
        assert int(row['inserted_guide_bases'])==hit.ins and int(row['deleted_guide_bases'])==hit.dele
        actual = None if row['cfd']=='NA' else float(row['cfd'])
        assert (actual is None)==(hit.score is None), (row,hit)
        if actual is not None: assert abs(actual-hit.score)<1e-12,(row,hit)
    assert found == set(expected), ('missing anchors',set(expected)-found)
    return len(rows)


def paf_rows(text, expected):
    found = set()
    for line in text.splitlines():
        f = line.split('\t'); tags={x.split(':',2)[0]:x.split(':',2)[2] for x in f[12:]}
        key = f[5],f[4],int(f[8]) if f[4]=='+' else int(f[7])
        assert key not in found;found.add(key)
        hit = expected[key]
        assert (int(f[7]),int(f[8]),tags['cg'],int(tags['nm'])) == (hit.start,hit.end,hit.cigar,hit.mm)
        assert ('cf' in tags)==(hit.score is not None)
        if hit.score is not None: assert abs(float(tags['cf'])-hit.score)<=.000051
    assert found == set(expected), ('missing anchors',set(expected)-found)
    return len(found)


def metrics(stderr):
    values=[json.loads(x[len('CFD_METRICS '):]) for x in stderr.splitlines() if x.startswith('CFD_METRICS ')]
    assert len(values)==1
    return values[0]


def read_fasta(path):
    refs={};name=None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):name=line[1:].split()[0];assert name not in refs;refs[name]=''
        else:refs[name]+=''.join(line.split()).upper()
    return refs


def main():
    p=argparse.ArgumentParser()
    for key in ('cli','old-cli','native','old-native','builder','out'):
        p.add_argument('--'+key,required=True)
    p.add_argument('--real-panel-root',type=Path,help='optional retained /tmp/crispr-v2-3mb; only 50 kb prefixes used')
    a=p.parse_args();out=Path(a.out);out.mkdir(exist_ok=False,parents=True)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    os.environ['PYTHONDONTWRITEBYTECODE']='1'
    tables=load_tables();g,resolved,nrefs=fixtures()
    datasets=[('resolved',g,resolved),('N',g,nrefs)]
    source_info=[]
    if a.real_panel_root:
        regions=json.loads((ROOT/'research/pam_anchored_poc/scalability_v2/results/3mb/regions.json').read_text())
        for region in regions:
            path=a.real_panel_root/region['name']/'ref.fa'
            assert digest(path)==region['ref_sha256']
            seq=next(iter(read_fasta(path).values()))[:50000]
            # Exact intended target from this bounded source, not an invented sequence.
            end=next(i for i in range(40,len(seq)-3) if seq[i] in 'ACGT' and seq[i+1:i+3]=='GG')
            guide=seq[end-20:end]
            datasets.append((region['name']+'_50kb',guide,{region['name']+'_prefix_50kb':seq}))
            source_info.append(dict(source=str(path),source_sha256=digest(path),region=region['region'],
                                    extracted_bases=len(seq),guide=guide,oriented_exact_pam_anchor=end))
    configs=[(5,0,.75,'GG'),(5,1,.75,'GG'),(5,2,.75,'GG'),(8,2,.5,'GG'),(20,0,0,'GG'),(0,2,.75,'GG')]
    summary=dict(baseline_commit='252a6bb',host=platform.node(),python=platform.python_version(),
                 bounded_real_sources=source_info,comparisons=[],performance=[],special_checks=[],
                 resource_scope='native ru_maxrss plus 20 ms sampled summed process tree (includes supervisor); not unique physical memory')
    paths=[Path(a.cli),Path(a.old_cli),Path(a.native),Path(a.old_native),Path(a.builder),Path(__file__),
           Path(__file__).with_name('native.cpp'),ROOT/'mismatch_scores.txt',ROOT/'pam_scores.txt']
    summary['sha256']={str(p):digest(p) for p in paths}
    nativeargs=lambda idx,ref,guides,conf,mode,tau,mm=ROOT/'mismatch_scores.txt',pm=ROOT/'pam_scores.txt':[
        a.native,idx,ref,guides,str(conf[0]),str(conf[1]),str(conf[2]),conf[3],str(int(conf[1]>0)),mode,str(tau),mm,pm]
    saved={};began=time.perf_counter()
    with tempfile.TemporaryDirectory(prefix='crispr-cfd-native-') as temporary:
        tmp=Path(temporary)
        for label,guide,refs in datasets:
            panel=tmp/label;panel.mkdir();ref=panel/'ref.fa';guides=panel/'guides.fa';idx=panel/'index'
            ref.write_text(''.join(f'>{c}\n{s}\n' for c,s in refs.items()));guides.write_text('>g\n'+guide+'\n')
            with (out/(label+'-index.log')).open('x') as log:
                subprocess.run([sys.executable,ROOT/'research/pam_anchored_poc/compatibility/index_manifest.py',
                    'build','--builder',a.builder,'--reference',ref,'--index',idx],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=30)
            summary.setdefault('panels',[]).append(dict(label=label,guide=guide,contigs=len(refs),bases=sum(map(len,refs.values())),reference_sha256=digest(ref)))
            # Native/CLI inputs remain in temporary storage only; no genomic files committed.
            saved[label]=(panel,ref,guides,idx,refs,guide)
            for conf in configs if label in ('resolved','N') else configs[:3]:
                m,z,frac,pam=conf
                full=representatives(oracle(refs,guide,pam,m,z,frac,tables))
                basecmd=['--search-backend','pam-anchored','--pam-search-bin',a.old_native,'--columba-index',idx,
                         '-r',ref,'-G',guides,'-t','1','-m',str(m),'-b',str(int(z>0)),'-z',str(z),'-f',str(frac),'-p',pam]
                before=command([a.old_cli]+basecmd).stdout
                current=command([a.cli]+basecmd).stdout
                assert before==current, ('default output changed',label,conf)
                paf_rows(current,full)
                if m==0 and z==2 and label=='resolved':
                    assert any(x.score is None for x in full.values()), 'mixed-only unavailable representatives must be covered'
                row=dict(panel=label,config=conf,unfiltered_anchors=len(full),default_byte_identical=True,
                         unfiltered_paf_sha256=hashlib.sha256(current.encode()).hexdigest(),thresholds=[])
                for tau in (0,.1,.2,.5):
                    expected={k:v for k,v in full.items() if keep(v,tau)}
                    filtered=command([a.cli]+basecmd+['--min-cfd',str(tau)]).stdout
                    paf_rows(filtered,expected)
                    control=command(nativeargs(idx,ref,guides,conf,'unpruned',tau))
                    pruned=command(nativeargs(idx,ref,guides,conf,'pruned',tau))
                    native_rows(control.stdout,expected);native_rows(pruned.stdout,expected)
                    am,bm=metrics(control.stderr),metrics(pruned.stderr)
                    if label=='N':assert bm['pruned']==0 and bm['states']==am['states']
                    if tau==0:assert filtered==current and bm['states']==am['states']
                    row['thresholds'].append(dict(threshold=tau,anchors=len(expected),
                        unavailable=sum(x.score is None for x in expected.values()),
                        final_anchor_or_representative_differences=0,control=am,pruned=bm))
                summary['comparisons'].append(row)
                print(label,conf,'anchors',len(full),'A/B/C PASS',flush=True)
            # Performance: same native binary/matrices/index, modes alternate over 3 repetitions.
            for conf in configs[:3]:
                full=representatives(oracle(refs,guide,'GG',conf[0],conf[1],conf[2],tables))
                expected={k:v for k,v in full.items() if keep(v,.1)}
                record=dict(panel=label,config=conf,threshold=.1,runs=[])
                for repetition in range(3):
                    order=['unpruned','overhead','pruned']
                    if repetition%2:order.reverse()
                    for mode in order:
                        prefix=out/f'perf-{label}-z{conf[1]}-{repetition}-{mode}'
                        cmd=[sys.executable,ROOT/'research/pam_anchored_poc/monitored_pilot/supervise.py',
                             '--limit-kib','1572864','--interval','.02','--log',str(prefix)+'.rss',
                             '--summary',str(prefix)+'.monitor.json','--']+nativeargs(idx,ref,guides,conf,mode,.1)
                        r=command(cmd)
                        native_rows(r.stdout[:r.stdout.rfind('\n{')+1] if '\n{' in r.stdout else r.stdout,expected)
                        measurement=metrics(r.stderr)
                        monitor=json.loads(Path(str(prefix)+'.monitor.json').read_text())
                        record['runs'].append(dict(repetition=repetition,mode=mode,native=measurement,monitor=monitor))
                summary['performance'].append(record)
        # Additional all-PAM, non-20, mixed-only, and matrix-fallback checks on tiny panels.
        panel,ref,guides,idx,refs,guide=saved['resolved']
        for pam in [x+y for x in 'ACGT' for y in 'ACGT' if x+y!='GG']:
            conf=(5,2,.75,pam);full=representatives(oracle(refs,guide,pam,5,2,.75,tables))
            for tau in (0,.1,.2,.5):
                expected={k:v for k,v in full.items() if keep(v,tau)}
                r=command(nativeargs(idx,ref,guides,conf,'pruned',tau));native_rows(r.stdout,expected)
            summary['special_checks'].append(dict(kind='PAM',pam=pam,thresholds=4,status='pass'))
        # A guide length unsupported by the simple bound must use a complete fallback.
        guides.write_text('>g\n'+guide[:19]+'\n')
        full=representatives(oracle(refs,guide[:19],'GG',5,2,.75,tables))
        r=command(nativeargs(idx,ref,guides,(5,2,.75,'GG'),'pruned',.2))
        native_rows(r.stdout,{k:v for k,v in full.items() if keep(v,.2)})
        assert metrics(r.stderr)['pruned']==0
        summary['special_checks'].append(dict(kind='non20_fallback',status='pass'))
        guides.write_text('>g\n'+guide+'\n')
        # Contract-violating but valid custom matrices must disable bounds, not cap mismatch settings.
        for name,change in [('above_one',lambda s:s.replace('rA:dG,20 0.22727272699999998','rA:dG,20 1.2')),
                            ('missing_key',lambda s:'\n'.join(x for x in s.splitlines() if not x.startswith('rA:dG,20 '))+'\n')]:
            original=(ROOT/'mismatch_scores.txt').read_text();modified=change(original);assert modified!=original
            custom=panel/(name+'.txt');custom.write_text(modified)
            conf=(8,2,.5,'GG')
            aa=command(nativeargs(idx,ref,guides,conf,'unpruned',.2,custom))
            bb=command(nativeargs(idx,ref,guides,conf,'pruned',.2,custom))
            # Identical TSV payloads (constructor diagnostic logs are ignored).
            header='guide\tcontig\t';assert aa.stdout[aa.stdout.index(header):]==bb.stdout[bb.stdout.index(header):]
            assert not metrics(bb.stderr)['matrix_contract'] and metrics(bb.stderr)['pruned']==0
            summary['special_checks'].append(dict(kind='matrix_fallback',case=name,status='pass'))
        # Fail closed for malformed invocation or broken binding, never publish a TSV header.
        bad=subprocess.run([a.native],capture_output=True,text=True);assert bad.returncode!=0
        with ref.open('a') as f:f.write('A\n')
        bad=subprocess.run(list(map(str,nativeargs(idx,ref,guides,(5,2,.75,'GG'),'pruned',.1))),capture_output=True,text=True)
        assert bad.returncode!=0 and 'guide\tcontig\t' not in bad.stdout
        summary['special_checks'].append(dict(kind='invalid_protocol_and_changed_reference_rejected',status='pass'))
    summary['elapsed_seconds']=time.perf_counter()-began
    for path,expected in summary['sha256'].items():assert digest(path)==expected
    (out/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print('ALL PASS',len(summary['comparisons']),'configurations',len(summary['performance']),'performance cases',summary['elapsed_seconds'])


if __name__=='__main__':
    main()
