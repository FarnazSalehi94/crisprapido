#!/usr/bin/env python3
"""Review/apply a narrow corrupted-index deletion, never a recursive removal.

Requires independent post-promotion validation + a complete supervised SBDS
search and intact checksums. Review writes a pinned, exclusive manifest; apply
requires its SHA256, repeats all checks, and only unlinks the reviewed regular
files inside indices_quarantine_v1/<assembly>-original. Healthy files outside
that exact quarantined index directory are never eligible. A failure preserves
remaining files and logs PARTIAL; freed bytes are logical removed bytes, not
an assertion of physical MooseFS reclamation.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import time

from validate_replacement import SUFFIXES, digest, require


def checked_json(path):
    return json.loads(path.read_text())


def review(root, assembly, validation_path, monitor_path):
    require('/' not in assembly and assembly not in ('', '.', '..'), 'invalid assembly')
    index = root/'indices'/assembly
    promotion = checked_json(index/'promotion.json')
    require(promotion['status']=='PROMOTED' and promotion['assembly']==assembly, 'not promoted')
    quarantine = root/'indices_quarantine_v1'/(assembly+'-original')
    require(not quarantine.is_symlink() and quarantine.is_dir(), 'quarantine directory absent/symlink')
    require(Path(promotion['original_quarantined'])==quarantine, 'quarantine path mismatch')
    completion_path=index/'completion.json'
    completion=checked_json(completion_path)
    validation=checked_json(validation_path)
    require(validation['status']=='PASS' and validation['assembly']==assembly, 'replacement validation absent')
    require(validation['completion_sha256']==digest(completion_path), 'replacement completion changed')
    require(validation['components_sha256']==completion['components_sha256'], 'component identities differ')
    cross=completion['cross_node_recheck']
    require(cross['result']=='PASS' and cross['host']!=cross['build_host'] and bool(cross['build_host']), 'cross-host gate')
    for s in SUFFIXES:
        require(digest(index/(assembly+'.'+s))==completion['components_sha256'][s], 'replacement drift: '+s)
    runs=root/'runs/sbds_ngg_nrr_v1'
    binding=runs/'00_manifest/bindings'/assembly/'index.crispr-binding.json'
    require(digest(binding)==validation['binding_sha256'], 'binding changed after validation')
    binding_json=checked_json(binding)
    require(binding_json['components_sha256']==completion['components_sha256'], 'binding index mismatch')
    source_reference=binding.parent/'reference.fa'
    require(digest(source_reference)==binding_json['reference_sha256'], 'bound reference drift')
    success_path=runs/'02_search'/assembly/'SUCCESS.json'
    success=checked_json(success_path)
    require(success['status']=='COMPLETE' and success['assembly']==assembly, 'SBDS incomplete')
    require(set(success['pams'])=={'NGG','NAG','NGA','NAA'}, 'missing PAM stages')
    for path,spec in success['outputs'].items():
        require(digest(path)==spec['sha256'], 'search output drift: '+path)
    monitor=checked_json(monitor_path)
    require(monitor['reason']=='completed' and monitor['child_returncode']==0, 'supervisor not complete')
    # Monitor binding is explicit: filename includes scheduler array/task, supplied
    # after sacct lookup of this SUCCESS attempt's child job ID.
    job=Path(success['attempt']).name.removeprefix('attempt-')
    require(job in monitor.get('validated_child_job_ids', []), 'monitor not tied to completed attempt')
    attempt=Path(success['attempt'])
    for pam in success['pams']:
        c=checked_json(attempt/(pam+'.summary.json'))
        require(c['audit_status']=='PASS', 'audit failure')
        counts=c['counts']
        require(counts['rows']==counts['cfd_defined'], 'missing CFD')
        require(counts['rows']==sum(v for k,v in counts.items() if k.startswith('cb_')), 'neural availability accounting')
    error_path=runs/'00_manifest/verification-batch-1'/(assembly+'.json')
    error=checked_json(error_path)
    require(error['assembly']==assembly and error['status']=='FAILED' and error['failure_kind']=='EIO',
            'no confirmed corrupted-original evidence; healthy originals ineligible')
    allowed={assembly+'.'+s for s in SUFFIXES}|{'index_build.time'}
    entries=[]
    for p in sorted(quarantine.iterdir()):
        st=p.lstat()
        require(p.name in allowed and stat.S_ISREG(st.st_mode) and st.st_nlink==1, 'unexpected quarantine entry: '+str(p))
        entries.append(dict(name=p.name,bytes=st.st_size,inode=st.st_ino,mtime_ns=st.st_mtime_ns))
    require({assembly+'.'+s for s in SUFFIXES} <= {e['name'] for e in entries}, 'quarantine components missing')
    return dict(assembly=assembly,root=str(root),quarantine=str(quarantine),entries=entries,
        validation=str(validation_path),validation_sha256=digest(validation_path),
        monitor=str(monitor_path),monitor_sha256=digest(monitor_path),
        success_sha256=digest(success_path),completion_sha256=digest(completion_path),
        error_report=error,logical_removed_bytes=sum(e['bytes'] for e in entries),
        note='Only confirmed corrupted-original regular files reviewed; no FASTA or historical result eligible.')


def apply(manifest_path, expected_sha256):
    require(digest(manifest_path)==expected_sha256, 'review manifest changed')
    m=checked_json(manifest_path)
    current=review(Path(m['root']),m['assembly'],Path(m['validation']),Path(m['monitor']))
    require(current==m, 'live evidence differs from review; nothing deleted')
    log=manifest_path.with_suffix('.deletion.json')
    require(not log.exists(), 'deletion log exists')
    events=dict(status='STARTED',review_sha256=expected_sha256,deleted=[],logical_removed_bytes=0)
    def save():
        tmp=log.with_suffix('.tmp')
        with tmp.open('w') as f:
            json.dump(events,f,indent=2); f.write('\n'); f.flush(); os.fsync(f.fileno())
        os.replace(tmp,log)
    save()
    try:
        for e in m['entries']:
            p=Path(m['quarantine'])/e['name']
            st=p.lstat()
            require(stat.S_ISREG(st.st_mode) and st.st_nlink==1 and st.st_ino==e['inode']
                    and st.st_size==e['bytes'] and st.st_mtime_ns==e['mtime_ns'], 'quarantine changed during deletion')
            p.unlink()
            events['deleted'].append(e['name']); events['logical_removed_bytes']+=e['bytes'];save()
        Path(m['quarantine']).rmdir()
        events.update(status='COMPLETE',finished_unix=time.time())
        save()
    except BaseException as exc:
        events.update(status='PARTIAL',error=str(exc));save();raise
    return events


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path)
    p.add_argument('--assembly');p.add_argument('--validation',type=Path);p.add_argument('--monitor',type=Path)
    p.add_argument('--manifest',type=Path,required=True);p.add_argument('--apply-sha256')
    a=p.parse_args()
    if a.apply_sha256: print(json.dumps(apply(a.manifest,a.apply_sha256)))
    else:
        result=review(a.root,a.assembly,a.validation,a.monitor)
        with a.manifest.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
        print(json.dumps(dict(status='REVIEWED',manifest=str(a.manifest),sha256=digest(a.manifest),bytes=result['logical_removed_bytes'])))
