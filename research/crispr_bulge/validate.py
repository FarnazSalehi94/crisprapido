#!/usr/bin/env python3
"""Bounded scoring-only validation and CPU overhead benchmark. Never searches."""
from collections import Counter
from itertools import groupby
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time

from adapter import IndexedFasta, convert, digest, encode_tokens, guides_from_fasta, parse_paf

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent


def run(cmd,cwd=ROOT,stdin=None,timeout=300):
    p=subprocess.run(list(map(str,cmd)),cwd=cwd,input=stdin,capture_output=True,text=True,timeout=timeout)
    if p.returncode:raise AssertionError((cmd,p.stderr[-4000:]))
    return p


def compressed(ops):
    return ''.join(str(len(list(v)))+k for k,v in groupby(ops))


def synthetic_panel(out,probe):
    g='GAAACAGTCGATTTTATCAC';cases=[]
    def add(name,ops,pam='AGG',dna='A',guide=g,expected='scored'):
        cases.append((name,ops,pam,dna,guide,expected))
    for pam in ('AGG','CGG','GGG','TGG'):add('exact_'+pam,'='*20,pam)
    for pos in range(20):
        for nt in 'ACGT':
            if nt!=g[pos]:add(f'mismatch_{pos}_{nt}','='*pos+'X'+'='*(19-pos),dna=nt)
    for pos in range(21):
        for nt in 'ACGT':add(f'D1_{pos}_{nt}','='*pos+'D'+'='*(20-pos),dna=nt)
    for pos in range(20):add(f'I1_{pos}','='*pos+'I'+'='*(19-pos))
    add('mm4_D1','XXXXD'+'='*16)
    add('mm4_I1','XXXXI'+'='*15)
    add('mm5_D1','XXXXX'+'D'+'='*15,expected='outside_documented_mismatch_domain')
    add('mm6','XXXXXX'+'='*14)
    add('mm7','XXXXXXX'+'='*13,expected='outside_documented_mismatch_domain')
    add('D2','='*5+'DD'+'='*15,expected='alignment_exceeds_24_columns')
    add('I2','='*5+'II'+'='*13,expected='unsupported_multibase_bulge')
    add('ID','='*5+'ID'+'='*14,expected='unsupported_mixed_bulge')
    add('DI','='*5+'DI'+'='*14,expected='unsupported_mixed_bulge')
    add('events2','==='+'I'+'===='+'I'+'='*11,expected='unsupported_multiple_bulge_events')
    add('neutral_gap_N','='+'D'+'='*19,dna='N',expected='ambiguous_reference')
    add('distal_gap_N','D'+'='*20,dna='N',expected='ambiguous_reference')
    add('guide19','='*19,guide=g[:-1],expected='unsupported_guide_length')
    add('guideN','='*20,guide='N'+g[1:],expected='ambiguous_guide')
    add('nonNRR','='*20,pam='AAC',expected='unsupported_pam')
    for pam in ('AAG', 'AGA', 'AAA'):
        add('nrr_'+pam,'='*20,pam=pam,expected='scored')
    add('pamN','='*20,pam='NGG',expected='ambiguous_pam')
    add('incomplete','='*20,pam='GG',expected='incomplete_pam')
    # Different valid alignments of the SAME target string on separate fixture
    # contigs. They are supplied representatives, not a new search experiment.
    records=[];guides={};refs={};expect={};probe_inputs=[];probe_indices=[]
    for name,ops,pam,dna,guide,status in cases:
        qi=0;target=[]
        for op in ops:
            if op=='I':qi+=1
            elif op=='D':target.append(dna)
            else:
                base=guide[qi];qi+=1
                target.append(base if op=='=' else ('C' if base!='C' else 'A'))
        target=''.join(target)
        if name.startswith('mismatch_'):
            i=int(name.split('_')[1]);target=target[:i]+dna+target[i+1:]
        assert qi==len(guide)
        guides[name]=guide
        for strand in '+-':
            contig=name+strand
            seq=target+pam
            if strand=='-':seq=seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
            refs[contig]=seq
            lo,hi=(0,len(target)) if strand=='+' else (len(pam),len(seq))
            cg=compressed(ops);events=len(re.findall(r'[ID]+',ops))
            f=list(map(str,[name,len(guide),0,len(guide),strand,contig,len(seq),lo,hi,
                           ops.count('='),len(ops),255]))
            f += [f'nm:i:{ops.count("X")}',f'ng:i:{events}',f'bs:i:{max([len(x) for x in re.findall("[ID]+",ops)] or [0])}','cg:Z:'+cg]
            index=len(records);records.append(f);expect[index]=status
            if len(pam)==3 and set(pam)<=set('ACGT'):
                probe_indices.append(index);probe_inputs.append('\t'.join((guide,target,cg,pam[1:])))
    # Real original CFD function, not invented fixture CFD values.
    result=run([probe],stdin='\n'.join(probe_inputs)+'\n').stdout.splitlines()
    assert len(result)==len(probe_indices)
    for index,value in zip(probe_indices,result):
        if value!='NA':records[index].append('cf:f:'+format(float(value),'.4f'))
    out.mkdir()
    ref=out/'ref.fa';guidefile=out/'guides.fa';paf=out/'input.paf'
    ref.write_text(''.join(f'>{name}\n{seq}\n' for name,seq in refs.items()))
    guidefile.write_text(''.join(f'>{name}\n{seq}\n' for name,seq in guides.items()))
    paf.write_text(''.join('\t'.join(f)+'\n' for f in records))
    return ref,guidefile,paf,expect


def independent_features(row,guides,refs):
    # Independent bounded in-memory reference reader, coordinate extraction and
    # expanded-CIGAR reconstruction, separate from the indexed adapter.
    f=row.rstrip().split('\t');g=guides[f[0]];seq=refs[f[5]];lo,hi=int(f[7]),int(f[8])
    comp=dict(zip('ACGTN','TGCAN'))
    rev=lambda s:''.join(comp[c] for c in reversed(s))
    target=seq[lo:hi] if f[4]=='+' else rev(seq[lo:hi])
    pam=seq[hi:hi+3] if f[4]=='+' else rev(seq[lo-3:lo])
    cg=next(x[5:] for x in f[12:] if x.startswith('cg:Z:'))
    ops=''.join(op*int(n) for n,op in re.findall(r'(\d+)([=XMID])',cg))
    qi=ti=0;ag=[];at=[]
    for op in ops:
        if op=='I':ag.append(g[qi]);at.append('-');qi+=1
        elif op=='D':ag.append('-');at.append(target[ti]);ti+=1
        else:ag.append(g[qi]);at.append(target[ti]);qi+=1;ti+=1
    assert qi==len(g) and ti==len(target)
    return dict(guide=''.join(ag)+'NGG',target=''.join(at)+pam), '\t'.join((g,target,cg,pam[1:]))


def fasta_dict(path):
    out={};name=None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):name=line[1:].split()[0];out[name]=''
        else:out[name]+=line.strip().upper()
    return out


def main():
    p=argparse.ArgumentParser()
    for key in ('out','model-root','model-python','cfd-probe','panel-reference','panel-guides','panel-paf'):
        p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    source_files=[HERE/'adapter.py',HERE/'worker.py',HERE/'score_paf.py',HERE/'official_reference.py',HERE/'model_manifest.json',Path(__file__),
                  a.cfd_probe,a.panel_reference,a.panel_guides,a.panel_paf]
    summary=dict(baseline_commit='178de00',panels=[],parity={},benchmarks=[],
                 sha256={str(x):digest(x) for x in source_files},model_root=str(a.model_root),model_python=str(a.model_python),
                 scope='bounded post-search scoring only; not a search experiment or biological validation')
    began=time.perf_counter()
    ref,guides,paf,expected=synthetic_panel(a.out/'synthetic-inputs',a.cfd_probe)
    panels=[('synthetic',ref,guides,paf,expected),('retained172',a.panel_reference,a.panel_guides,a.panel_paf,None)]
    requests=[];cfd_requests=[]
    for name,ref,guidefile,paf,expected in panels:
        indexed=IndexedFasta(ref);guide_dict=guides_from_fasta(guidefile);refs=fasta_dict(ref)
        counts=Counter();raws=paf.read_text().splitlines(keepends=True)
        for i,raw in enumerate(raws):
            _,state,features=convert(raw,guide_dict,indexed);counts[state]+=1
            if expected is not None:assert expected[i]==state,(i,expected[i],state)
            if features is not None:
                independent,cfd=independent_features(raw,guide_dict,refs)
                assert independent==features
                requests.append(dict(id=len(requests),**features));cfd_requests.append(cfd)
        indexed.close()
        for mode in ('cfd','both','crispr-bulge'):
            dest=a.out/(name+'-'+mode)
            cmd=[sys.executable,HERE/'score_paf.py','--input',paf,'--out-dir',dest,'--scoring',mode]
            if mode!='cfd':cmd+=['--reference',ref,'--guides',guidefile,'--expected-reference-sha256',digest(ref),
                                 '--model-root',a.model_root,'--model-python',a.model_python]
            run(cmd)
            actual=(dest/'output.paf').read_text().splitlines()
            assert len(actual)==len(raws)
            for original,line in zip(raws,actual):
                f,tags=parse_paf(original);af=line.split('\t')
                # Existing fields/tags and original canonical CIGAR stay unchanged.
                if mode=='cfd':assert line==original.rstrip('\r\n');continue
                expected_fields=f if mode=='both' else f[:12]+[x for x in f[12:] if not x.startswith('cf:')]
                assert af[:len(expected_fields)]==expected_fields
                status=next(x[5:] for x in af[12:] if x.startswith('bu:Z:'))
                assert any(x.startswith('cb:f:') for x in af[12:])==(status=='scored')
            if mode=='cfd':assert (dest/'output.paf').read_bytes()==paf.read_bytes()
        summary['panels'].append(dict(name=name,rows=len(raws),statuses=dict(counts),
            cfd_default_byte_identical=True,rows_cigar_order_preserved_all_modes=True,
            reference_sha256=digest(ref),guides_sha256=digest(guidefile),paf_sha256=digest(paf)))
    requestfile=a.out/'parity-input.jsonl';requestfile.write_text(''.join(json.dumps(x)+'\n' for x in requests))
    # Two adapter batch sizes, one original ensemble call, identical input order.
    for batch in (512,32):
        run([a.model_python,'-I',HERE/'worker.py','--model-root',a.model_root,'--input',requestfile,
             '--output',a.out/f'adapter-{batch}.jsonl','--metrics',a.out/f'adapter-{batch}-metrics.json','--batch-size',batch])
    run([a.model_python,'-I',HERE/'official_reference.py','--upstream',a.model_root,'--input',requestfile,
         '--output',a.out/'official.jsonl','--metrics',a.out/'official-metrics.json'])
    original=[json.loads(x) for x in (a.out/'official.jsonl').read_text().splitlines()]
    assert len(original)==len(requests)
    for batch in (512,32):
        actual=[json.loads(x) for x in (a.out/f'adapter-{batch}.jsonl').read_text().splitlines()]
        assert [x['id'] for x in actual]==[x['id'] for x in original]
        maximum=max(abs(x['score']-y['score']) for x,y in zip(actual,original))
        assert maximum<=1e-6,maximum
        summary['parity'][str(batch)]=dict(rows=len(actual),maximum_absolute_difference=maximum,tolerance=1e-6)
    import numpy as np
    categories=np.asarray(encode_tokens(requests),dtype='<i8')
    official_metrics=json.loads((a.out/'official-metrics.json').read_text())
    assert hashlib.sha256(categories.tobytes()).hexdigest()==official_metrics['categorical_sha256']
    assert official_metrics['dna_two_base_encoder_error']=='actual_seq_size > seq_len'
    summary['feature_parity']='exact categorical byte equality to official one-hot + argmax'
    single=next(r for r in requests if '-' in r['guide'] and any(x!=y for x,y in zip(r['guide'][:-3],r['target'][:-3]) if x!='-' and y!='-'))
    singlefile=a.out/'single-input.jsonl';singlefile.write_text(json.dumps(single)+'\n')
    run([a.model_python,'-I',HERE/'worker.py','--model-root',a.model_root,'--input',singlefile,
         '--output',a.out/'single-adapter.jsonl','--metrics',a.out/'single-adapter-metrics.json'])
    run([a.model_python,'-I',HERE/'official_reference.py','--upstream',a.model_root,'--input',singlefile,
         '--output',a.out/'single-official.jsonl','--metrics',a.out/'single-official-metrics.json'])
    single_a=json.loads((a.out/'single-adapter.jsonl').read_text());single_b=json.loads((a.out/'single-official.jsonl').read_text())
    assert single_a['id']==single_b['id'] and abs(single_a['score']-single_b['score'])<=1e-6
    summary['single_record_parity']=dict(maximum_absolute_difference=abs(single_a['score']-single_b['score']))
    summary['upstream_limits']=dict(dna_two_base_error=official_metrics['dna_two_base_encoder_error'],
                                    rna_two_base_encoder_shape=official_metrics['rna_two_base_encoder_shape'])
    # Identical supported alignments repeated for inference throughput, not new sites.
    bench=[dict(requests[i%len(requests)],id=i) for i in range(1024)]
    benchmark_input=a.out/'benchmark-input.jsonl';benchmark_input.write_text(''.join(json.dumps(x)+'\n' for x in bench))
    cfd_input=a.out/'benchmark-cfd.tsv';cfd_input.write_text(''.join(cfd_requests[i%len(cfd_requests)]+'\n' for i in range(1024)))
    for repeat in range(3):
        order=('cfd','adapter','official') if repeat%2==0 else ('official','adapter','cfd')
        for mode in order:
            prefix=a.out/f'bench-{repeat}-{mode}'
            if mode=='adapter':
                cmd=[a.model_python,'-I',HERE/'worker.py','--model-root',a.model_root,'--input',benchmark_input,
                     '--output',str(prefix)+'.pred.jsonl','--metrics',str(prefix)+'.worker.json']
            elif mode=='official':
                cmd=[a.model_python,'-I',HERE/'official_reference.py','--upstream',a.model_root,'--input',benchmark_input,
                     '--output',str(prefix)+'.pred.jsonl','--metrics',str(prefix)+'.worker.json']
            else:cmd=[a.cfd_probe]
            measured=[sys.executable,HERE/'measure.py','--metrics',str(prefix)+'.time.json',
                      '--stdout',str(prefix)+'.stdout','--stderr',str(prefix)+'.stderr']
            if mode=='cfd':measured+=['--stdin',cfd_input]
            cmd=[sys.executable,ROOT/'research/pam_anchored_poc/monitored_pilot/supervise.py','--limit-kib','4194304',
                 '--interval','.02','--log',str(prefix)+'.rss','--summary',str(prefix)+'.tree.json','--']+measured+['--']+cmd
            run(cmd,timeout=360)
            data=dict(repetition=repeat,mode=mode,rows=1024,
                measurement=json.loads(Path(str(prefix)+'.time.json').read_text()),
                process_tree=json.loads(Path(str(prefix)+'.tree.json').read_text()))
            if mode!='cfd':data['worker']=json.loads(Path(str(prefix)+'.worker.json').read_text())
            else:assert len(Path(str(prefix)+'.stdout').read_text().splitlines())==1024
            summary['benchmarks'].append(data)
    summary['benchmark_prediction_parity']=[]
    for repeat in range(3):
        left=[json.loads(x) for x in (a.out/f'bench-{repeat}-adapter.pred.jsonl').read_text().splitlines()]
        right=[json.loads(x) for x in (a.out/f'bench-{repeat}-official.pred.jsonl').read_text().splitlines()]
        assert len(left)==len(right)==1024
        assert [x['id'] for x in left]==[x['id'] for x in right]==list(range(1024))
        maximum=max(abs(x['score']-y['score']) for x,y in zip(left,right))
        assert maximum<=1e-6,maximum
        summary['benchmark_prediction_parity'].append(dict(repetition=repeat,rows=1024,maximum_absolute_difference=maximum))
    for path,value in summary['sha256'].items():assert digest(path)==value
    summary['elapsed_seconds']=time.perf_counter()-began
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(status='pass',panels=summary['panels'],parity=summary['parity'],elapsed_seconds=summary['elapsed_seconds']),indent=2))


if __name__=='__main__':main()
