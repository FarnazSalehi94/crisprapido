#!/usr/bin/env python3
import json,sys
x=json.load(open(sys.argv[1]));keys=('label','kind','strand','contig','global_start','start','anchor','indexed_genomic','reference_genomic','guide','true_mismatches')
with open(sys.argv[2],'w') as f:
 for case in x['cases']:
  assert all('\t' not in str(case[k]) and '\n' not in str(case[k]) for k in keys)
  f.write('\t'.join(str(case[k]) for k in keys)+'\n')
