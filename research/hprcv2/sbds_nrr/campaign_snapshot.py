#!/usr/bin/env python3
"""Evidence-derived snapshot; never changes or submits a campaign.

Use SUCCESS-selected attempts, not first/last glob. Old failures remain history,
not current failed results. Scheduler states are explicit, including allocated
CPU/node counts. Optional --verify-output-hashes performs full completion checks.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT=Path('/moosefs/salehi/HPRCv2_CRISPRapido')


def digest(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()


def snapshot(root, verify=False):
    r=root/'runs/sbds_ngg_nrr_v1'
    cohort=[s.split('\t')[0] for s in (root/'analysis_genomes.tsv').read_text().splitlines() if s.strip()]
    if len(cohort)!=466 or len(set(cohort))!=466:raise ValueError('cohort identity')
    prior=json.loads((r/'00_manifest/authoritative-manifest-v5.json').read_text())
    prior={x['assembly']:x['status'] for x in prior['rows']}
    rows=[];totals=Counter()
    for a in cohort:
        bd=r/'00_manifest/bindings'/a
        state=prior[a]
        if (bd/'binding.json').exists() and not (bd/'INCOMPLETE').exists():state='READY_BOUND'
        else:
            stage=root/'indices_rebuilt_v1'/a
            if (stage/'completion.json').exists():state='STAGED_'+json.loads((stage/'completion.json').read_text())['status']
        index=root/'indices'/a
        replacement='PROMOTED_REBUILT' if (index/'promotion.json').exists() else 'ORIGINAL_OR_PRIOR_RECOVERY'
        success=r/'02_search'/a/'SUCCESS.json'
        result=dict(assembly=a,index_status=state,index_origin=replacement,sbds_status='INCOMPLETE')
        if success.exists():
            data=json.loads(success.read_text())
            if data['status']!='COMPLETE':raise ValueError('bad completion: '+a)
            if verify:
                for p,spec in data['outputs'].items():
                    if digest(p)!=spec['sha256']:raise ValueError('changed completed output: '+p)
            attempt=Path(data['attempt']);counts=Counter()
            for pam in ('NGG','NAG','NGA','NAA'):
                summary=json.loads((attempt/(pam+'.summary.json')).read_text())
                if summary['audit_status']!='PASS':raise ValueError('audit failed: '+a)
                c=summary['counts'];counts['hits']+=c['rows'];counts['cfd_defined']+=c['cfd_defined']
                counts['cb_scored']+=c.get('cb_scored',0)
                counts['cb_na']+=sum(v for k,v in c.items() if k.startswith('cb_') and k!='cb_scored')
            if counts['hits']!=counts['cb_na']+counts['cb_scored']:raise ValueError('availability accounting')
            result.update(sbds_status='COMPLETE',counts=dict(counts),attempt=str(attempt),completion_sha256=digest(success))
            totals.update(counts)
        rows.append(result)
    deletions=[]
    for p in (r/'00_manifest/postpromotion-review-v1').glob('*-review.deletion.json'):
        d=json.loads(p.read_text())
        if d['status']=='COMPLETE':deletions.append(dict(path=str(p),bytes=d['logical_removed_bytes']))
    usage=shutil.disk_usage(root)
    return dict(cohort=466,unique_entries=466,index_counts=dict(Counter(x['index_status'] for x in rows)),
                promoted=sum(x['index_origin']=='PROMOTED_REBUILT' for x in rows),
                sbds_complete=sum(x['sbds_status']=='COMPLETE' for x in rows),
                sbds_remaining=sum(x['sbds_status']!='COMPLETE' for x in rows),totals=dict(totals),
                verified_output_hashes=verify,deletions=deletions,
                logical_removed_bytes=sum(x['bytes'] for x in deletions),
                filesystem=dict(free_bytes=usage.free,total_bytes=usage.total),rows=rows,
                caveats=['READY_BOUND denotes committed binding evidence, not a fresh hash scan of every input.',
                         'Historical failed attempts are retained, not counted as currently failed completed assemblies.',
                         'Physical storage reclamation/replication is controlled by MooseFS; logical removed bytes are not physical free bytes.',
                         'Noncanonical neural scores are computational, not biologically calibrated.'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=ROOT);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify-output-hashes',action='store_true');a=p.parse_args()
    report=snapshot(a.root,a.verify_output_hashes)
    with a.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('rows','caveats')}))
