#!/usr/bin/env python3
"""Fully rehash all 12 resolved index sets and sources before a pilot."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

SUFFIXES=('bwt','brt','rev.brt','sa.4','sa.bv.4','txt.bin','cct','pos','meta','fsid','headerSN.bin','sna')

def digest(path, limiter):
    h=hashlib.sha256();total=0
    with open(path,'rb') as stream:
        while True:
            chunk=stream.read(4*1024*1024)
            if not chunk:break
            h.update(chunk);total+=len(chunk);limiter[0]+=len(chunk)
            if limiter[1]:
                delay=limiter[0]/limiter[1]-(time.monotonic()-limiter[2])
                if delay>0:time.sleep(delay)
    return total,h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('assemblies','checksums','binary','freeze','output'):p.add_argument('--'+n,required=True)
    p.add_argument('--rate-mib',type=float,default=128);p.add_argument('--expected-assemblies',type=int,default=12)
    a=p.parse_args()
    if a.rate_mib<0 or a.expected_assemblies<1:p.error('invalid limits')
    out=Path(a.output);out.mkdir(exist_ok=False)
    sources=list(csv.DictReader(Path(a.assemblies).read_text().splitlines(),delimiter='\t'))
    components=list(csv.DictReader(Path(a.checksums).read_text().splitlines(),delimiter='\t'))
    ids={x['assembly_id'] for x in sources}
    if len(sources)!=a.expected_assemblies or len(ids)!=len(sources):raise ValueError('wrong assembly count or duplicate ID')
    bykey={}
    for row in components:
        key=(row['assembly'],row['path'])
        if key in bykey:raise ValueError('duplicate component manifest row')
        bykey[key]=row
    freeze=json.loads(Path(a.freeze).read_text()); limiter=[0,a.rate_mib*1024*1024,time.monotonic()]
    results=[]
    def verify(kind,assembly,path,expected,expected_size=None):
        if not re_hex(expected):raise ValueError('missing or malformed SHA256 for '+path)
        size,actual=digest(path,limiter)
        record={'kind':kind,'assembly':assembly,'path':path,'bytes':size,'sha256':actual,'expected':expected,
                'status':'PASS' if actual==expected and (expected_size is None or size==int(expected_size)) else 'FAIL'}
        results.append(record)
        with (out/'checks.jsonl').open('a') as f:f.write(json.dumps(record,sort_keys=True)+'\n')
        print(kind,assembly,record['status'],size,flush=True)
        if record['status']!='PASS':raise ValueError('hash/size mismatch: '+path)
    try:
        verify('binary','-',a.binary,freeze['binary_sha256'])
        for row in sources:
            asm=row['assembly_id'];prefix=row['index_prefix']
            if row['verification_status'] not in ('REBUILT_VERIFIED','ORIGINAL_VERIFIED'):raise ValueError('unverified assembly '+asm)
            for suffix in SUFFIXES:
                path=prefix+'.'+suffix;entry=bykey.get((asm,path))
                if not entry or entry['status']!='READABLE':raise ValueError('missing verified component '+path)
                verify('index',asm,path,entry['sha256'],entry['size'])
            verify('fasta',asm,row['fasta'],row['fasta_sha256'],row['fasta_bytes'])
            verify('fai',asm,row['fasta']+'.fai',row['fai_sha256'])
            verify('gzi',asm,row['fasta']+'.gzi',row['gzi_sha256'])
        if len(bykey)!=len(sources)*len(SUFFIXES):raise ValueError('unexpected extra component rows')
        summary={'status':'PASS','assemblies':len(sources),'components':len(sources)*len(SUFFIXES),
                 'files':len(results),'bytes_read':limiter[0],'elapsed_seconds':time.monotonic()-limiter[2],
                 'rate_mib_per_second':a.rate_mib}
        (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(out/'PASS').write_text('Full current hash verification passed.\n')
        return 0
    except (OSError,ValueError,KeyError) as exc:
        (out/'summary.json').write_text(json.dumps({'status':'FAIL','error':str(exc),'files_completed':len(results),'bytes_read':limiter[0]},indent=2)+'\n')
        print(str(exc),file=sys.stderr);return 1

def re_hex(value):
    return len(value)==64 and all(c in '0123456789abcdef' for c in value)
if __name__=='__main__':sys.exit(main())
