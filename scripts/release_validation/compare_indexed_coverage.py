import sys,json,re,collections
from pathlib import Path
from coverage_oracle import enumerate_sites
p=Path(sys.argv[1])
def fasta(path):
 d={};name=None
 for line in path.read_text().splitlines():
  if line.startswith('>'):name=line[1:].split()[0];d[name]=''
  else:d[name]+=line.strip()
 return d
refs=fasta(p/'ref.fa'); guides=fasta(p/'guides.fa')
truth=collections.defaultdict(set)
for name,g in guides.items():
 for chrom,seq in refs.items():
  for x in enumerate_sites(seq,g):
   category=('exact' if x['mismatches']==0 else 'mismatch') if x['bulge']=='none' else ('I' if x['bulge']=='RNA' else 'D') if x['mismatches']==0 else 'mismatch_bulge'
   truth[(name,chrom,x['start'],x['end'],x['strand'])].add(category)
def parse_sam(path):
 rows=[]
 for line in path.read_text().splitlines():
  if line.startswith('@'):continue
  f=line.split('\t');flag=int(f[1]);
  if flag&4:continue
  span=sum(int(n) for n,o in re.findall(r'(\d+)([MIDNSHP=X])',f[5]) if o in 'MDN=X')
  rows.append((f[0],f[2],int(f[3])-1,int(f[3])-1+span,'-' if flag&16 else '+'))
 return set(rows)
def paf(path):
 result=collections.defaultdict(set)
 for line in path.read_text().splitlines():
  f=line.split('\t'); ops=dict((x[:2],x[5:]) for x in f[12:]); cigar=ops['cg']; mm='X' in cigar; gap='I' in cigar or 'D' in cigar
  cls='mismatch_bulge' if mm and gap else 'I' if 'I' in cigar else 'D' if 'D' in cigar else 'mismatch' if mm else 'exact'
  result[(f[0],f[5],int(f[7]),int(f[8]),f[4])].add(cls)
 return result
pred=paf(p/'out.paf'); exact=parse_sam(p/'exact.sam'); edit=parse_sam(p/'edit.sam');
for cls in ('exact','mismatch','I','D','mismatch_bulge'):
 expected={k for k,v in truth.items() if cls in v}; actual={k for k,v in pred.items() if cls in v}; tp=expected&actual
 print(json.dumps({'category':cls,'expected':len(expected),'predicted':len(actual),'tp':len(tp),'recall':len(tp)/len(expected) if expected else None,'precision':len(tp)/len(actual) if actual else None,'missed':sorted(expected-actual),'unexpected':sorted(actual-expected)},sort_keys=True))
expected={k for k,v in truth.items() if 'exact' in v}
strict_exact = set()
for key in exact:
 name, chrom, start, end, strand = key
 if chrom not in refs or name not in guides or end-start != len(guides[name]):
  continue
 seq = refs[chrom]
 if strand == '+':
  target, pam = seq[start:end], seq[end:end+3]
 else:
  from coverage_oracle import rc
  target, pam = rc(seq[start:end]), rc(seq[max(0,start-3):start])
 if target == guides[name] and len(pam) == 3 and set(pam) <= set('ACGT') and pam[1:] == 'GG':
  strict_exact.add(key)
print(json.dumps({'stage':'raw_exact','truth':len(expected),'found':len(expected&strict_exact),'strict_predictions':len(strict_exact),'missed':sorted(expected-strict_exact),'unexpected_strict':sorted(strict_exact-expected)}))
print(json.dumps({'stage':'raw_edit_exact_interval','truth_loci':len(truth),'found':len(set(truth)&edit),'missed':sorted(set(truth)-edit)}))
(p/'truth.json').write_text(json.dumps([{'key':k,'categories':sorted(v)} for k,v in sorted(truth.items())],indent=2)+'\n')

# Record every category-level miss and each overlapping raw/final interval.
# An overlap is descriptive evidence, not automatic equivalence or proof of cause.
discrepancies = []
for key, categories in sorted(truth.items()):
 for category in sorted(categories - pred.get(key, set())):
  def overlaps(other):
   return (other[0] == key[0] and other[1] == key[1] and other[4] == key[4]
           and other[2] < key[3] and key[2] < other[3])
  discrepancies.append({'key':key, 'category':category,
    'raw_edit_exact_interval':key in edit,
    'overlap_raw_edit':[x for x in sorted(edit) if overlaps(x)],
    'overlap_final':[{'key':x,'categories':sorted(pred[x])} for x in sorted(pred) if overlaps(x)]})
(p/'discrepancies.json').write_text(json.dumps(discrepancies, indent=2)+'\n')
