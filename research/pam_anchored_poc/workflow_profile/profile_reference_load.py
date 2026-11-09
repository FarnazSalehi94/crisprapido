#!/usr/bin/env python3
"""Compare equivalent reference reconstruction strategies on a verified FASTA."""
import argparse,hashlib,json,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--mode',choices=['read_text','stream'],required=True);p.add_argument('--out',required=True);a=p.parse_args();seqs={};name=None;start=time.perf_counter()
if a.mode=='read_text':
 lines=iter(Path(a.reference).read_text().splitlines())
else:lines=Path(a.reference).open()
for line in lines:
 if line.startswith('>'):
  name=line[1:].split()[0]
  if name in seqs:raise ValueError('duplicate reference record')
  seqs[name]=[]
 elif name is not None:seqs[name].append(''.join(line.split()).upper())
seqs={name:''.join(chunks) for name,chunks in seqs.items()}
seconds=time.perf_counter()-start;start=time.perf_counter();h=hashlib.sha256()
for name,seq in seqs.items():h.update(name.encode('ascii')+b'\0');h.update(seq.encode('ascii')+b'\0')
result={'mode':a.mode,'load_seconds':seconds,'digest_seconds':time.perf_counter()-start,'contigs':len(seqs),'bases':sum(map(len,seqs.values())),'ordered_record_sha256':h.hexdigest()};Path(a.out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
