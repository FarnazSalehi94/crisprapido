#!/usr/bin/env python3
"""Inspect a bounded sample of serialized SA offsets without any index search."""
import argparse,bisect,json,struct
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--stage',required=True);p.add_argument('--out',required=True);a=p.parse_args();d=Path(a.stage)
if not (d/'index.crispr-binding.json').is_file() or (d/'.incomplete').exists():raise ValueError('stage is not a completed import')
raw=(d/'index.sa.4').open('rb').read(4096);v=struct.unpack('<'+'I'*(len(raw)//4),raw);pr=(d/'index.pos').read_bytes();off=struct.unpack('<'+'I'*(len(pr)//4),pr);sna=(d/'index.sna').read_bytes();names=[];i=0
while i<len(sna):
 n=struct.unpack_from('<Q',sna,i)[0];i+=8;names.append(sna[i:i+n].decode());i+=n
assert len(names)+1==len(off) and all(x<=off[-1] for x in v)
large=[x for x in v if 2**31<x<off[-1]]
result={'sampled_sa_entries':len(v),'above_INT32_MAX_inside_contigs':len(large),'examples':[{'global':x,'contig':names[bisect.bisect_right(off,x)-1],'offset':x-off[bisect.bisect_right(off,x)-1]} for x in large[:3]],'scope':'raw serialized SA values only; does not call findSA or validate search'}
Path(a.out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
