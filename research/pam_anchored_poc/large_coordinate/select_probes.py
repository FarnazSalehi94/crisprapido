#!/usr/bin/env python3
"""Select bounded, exact index-text probes from known FASTA coordinates.

Uses the trusted FAI for random access to the imported byte-identical plain FASTA;
reads only small windows and does not perform a guide search.
"""
import argparse,bisect,hashlib,json,mmap,re,struct
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--stage',required=True);p.add_argument('--fai',required=True);p.add_argument('--out',required=True);p.add_argument('--window',type=int,default=12000);a=p.parse_args();stage=Path(a.stage);rows=[]
expected=json.loads((stage/'import.json').read_text())['components']['fai']['sha256']
if hashlib.sha256(Path(a.fai).read_bytes()).hexdigest()!=expected:raise ValueError('FAI SHA256 does not match imported evidence')
for line in Path(a.fai).read_text().splitlines():
 x=line.split('\t');rows.append({'name':x[0],'length':int(x[1]),'byte_offset':int(x[2]),'line_bases':int(x[3]),'line_width':int(x[4])})
raw=(stage/'index.pos').read_bytes();positions=struct.unpack('<'+'I'*(len(raw)//4),raw);assert len(positions)==len(rows)+1 and [b-a for a,b in zip(positions,positions[1:])]==[x['length'] for x in rows]
seqf=(stage/'reference.fa').open('rb');idxf=(stage/'index.txt.bin').open('rb');rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def fetch(i,start,end):
 r=rows[i];assert 0<=start<end<=r['length'];begin=r['byte_offset']+(start//r['line_bases'])*r['line_width']+start%r['line_bases'];endoff=r['byte_offset']+((end-1)//r['line_bases'])*r['line_width']+(end-1)%r['line_bases']+1;seqf.seek(begin);b=seqf.read(endoff-begin);seq=b.replace(b'\n',b'').replace(b'\r',b'').decode('ascii').upper();assert len(seq)==end-start;idxf.seek(4+positions[i]+start);indexed=idxf.read(end-start).decode('ascii');assert len(indexed)==len(seq) and all(x==y for x,y in zip(seq,indexed) if x in 'ACGT');return seq,indexed
def candidates(i,lo,hi,strand,kind):
 ref,idx=fetch(i,lo,hi);out=[]
 if kind!='resolved' and 'N' not in ref:return out
 # Only inspect candidate indexed NGG positions, using C regex scanning.
 pattern=r'(?=[ACGT]GG)' if strand=='+' else r'(?=CC[ACGT])'
 for m in re.finditer(pattern,idx):
  j=m.start()-20 if strand=='+' else m.start()
  if j<0 or j+23>len(ref):continue
  g=ref[j:j+23];t=idx[j:j+23];r=g if strand=='+' else rc(g);z=t if strand=='+' else rc(t);pam=r[20:]
  strict=pam[0] in 'ACGT' and pam[1:]=='GG'
  if kind=='resolved':ok=strict and 'N' not in r and len(set(z[k:k+2] for k in range(19)))>=8
  elif kind=='target_n':ok=strict and 1<=r[:20].count('N')<=5 and set(r[:20])<=set('ACGTN')
  elif kind=='pam_required_n':ok=not strict and 'N' in pam[1:] and pam[0] in 'ACGT'
  else:ok=not strict and pam[0]=='N' and pam[1:]=='GG'
  if ok:out.append((lo+j,g,t,r,z))
 return out

cases=[]
def add_one(label,i,lo,hi,strand,kind,goal=None):
 possible=candidates(i,lo,hi,strand,kind)
 if not possible:return False
 if goal is None:goal=(lo+hi)//2
 chosen=[]
 for item in sorted(possible,key=lambda q: (len({q[4][k:k+2] for k in range(19)})<10 if kind=='resolved' else False, abs(q[0]-goal))):
  if all(abs(item[0]-other[0])>=12 for other in chosen):chosen.append(item)
  if len(chosen)==3:break
 for start,g,t,r,z in chosen:
  anchor=start+20 if strand=='+' else start+3;cases.append({'label':label,'kind':kind,'contig':rows[i]['name'],'contig_index':i,'contig_length':rows[i]['length'],'global_start':positions[i]+start,'start':start,'end':start+23,'strand':strand,'anchor':anchor,'reference_genomic':g,'indexed_genomic':t,'oriented_reference':r,'oriented_index':z,'guide':z[:20],'pam':r[20:],'true_mismatches':sum(x!=y for x,y in zip(z[:20],r[:20])),'strict_pam':r[20] in 'ACGT' and r[21:23]=='GG'})
 return True

N=len(rows);pivot=2**31;mid=bisect.bisect_right(positions,pivot)-1
if not 0<=mid<N:raise ValueError('reference does not cross 2^31; use fixture mode separately')
for label,i,center in [('beginning',0,0),('above_2_31',mid,pivot-positions[mid]),('end',N-1,rows[-1]['length'])]:
 lo=max(0,center-a.window//2);hi=min(rows[i]['length'],center+a.window//2)
 for strand in '+-':
  if not add_one(label,i,lo,hi,strand,'resolved',center):raise ValueError('no resolved strict PAM candidate: '+label+' '+strand)
# Additional near-end probes in the preceding contig provide alternatives
# if the final contig's repetitive exact patterns exceed the SA cap.
last_alt=N-2;alt_len=rows[last_alt]['length']
for strand in '+-':
 if not add_one('end_fallback',last_alt,max(0,alt_len-20000),alt_len,strand,'resolved',alt_len):raise ValueError('no near-end fallback '+strand)
# Find a real contig junction with resolved, nontrivial PAM patterns on both
# sides in both orientations; keep the probes within 1 kb of the junction.
boundary=None
for i in range(N-1):
 if rows[i]['length']<1000 or rows[i+1]['length']<1000:continue
 L=rows[i]['length'];tests=[(i,max(0,L-1000),L,'+'),(i,max(0,L-1000),L,'-'),(i+1,0,1000,'+'),(i+1,0,1000,'-')]
 if all(candidates(j,lo,hi,strand,'resolved') for j,lo,hi,strand in tests):boundary=i;break
if boundary is None:raise ValueError('no suitable bounded contig junction')
for label,j,lo,hi,goal in [('boundary_left',boundary,rows[boundary]['length']-1000,rows[boundary]['length'],rows[boundary]['length']),('boundary_right',boundary+1,0,1000,0)]:
 for strand in '+-':assert add_one(label,j,lo,hi,strand,'resolved',goal)
# Control: an indexed 23-mer assembled across the chosen contig junction.
# It may exist in the FM index, but must never be emitted as an alignment.
for strand in '+-':
 found_cross=False
 for junction in [boundary]+[x for x in range(N-1) if x!=boundary]:
  L=rows[junction]['length']
  if L<23 or rows[junction+1]['length']<23:continue
  for left in range(1,23):
   a1,i1=fetch(junction,L-left,L);a2,i2=fetch(junction+1,0,23-left);g=a1+a2;t=i1+i2;r=g if strand=='+' else rc(g);z=t if strand=='+' else rc(t)
   if r[20] in 'ACGT' and r[21:23]=='GG' and 'N' not in r and len(set(z[k:k+2] for k in range(19)))>=8:
    start=L-left;cases.append({'label':'cross_boundary','kind':'cross_boundary','contig':rows[junction]['name'],'contig_index':junction,'contig_length':L,'global_start':positions[junction]+start,'start':start,'end':start+23,'strand':strand,'anchor':start+20 if strand=='+' else start+3,'reference_genomic':g,'indexed_genomic':t,'oriented_reference':r,'oriented_index':z,'guide':z[:20],'pam':r[20:],'true_mismatches':0,'strict_pam':True});found_cross=True;break
  if found_cross:break
 if not found_cross:print('no strict crossing-boundary test for '+strand)
for label,lo,hi in [('before_2_31',max(0,pivot-positions[mid]-4000),pivot-positions[mid]-23),('after_2_31',pivot-positions[mid],min(rows[mid]['length'],pivot-positions[mid]+4000))]:
 for strand in '+-':
  if not add_one(label,mid,lo,hi,strand,'resolved',pivot-positions[mid]):raise ValueError('no 2^31 side candidate: '+label+' '+strand)
# Locate N blocks only in the first 1 GB of the reference, then inspect bounded
# windows around each block. This is a byte scan for N, not a guide search.
wanted={(kind,strand) for kind in ('target_n','pam_required_n','pam_wildcard_n') for strand in '+-'}
byte_offsets=[r['byte_offset'] for r in rows]
with (stage/'reference.fa').open('rb') as nf:
 with mmap.mmap(nf.fileno(),0,access=mmap.ACCESS_READ) as mapped:
  cursor=byte_offsets[0];limit=min(len(mapped),1000000000);clusters=0
  while wanted and cursor<limit and clusters<100:
   found=mapped.find(b'N',cursor,limit)
   if found<0:break
   i=bisect.bisect_right(byte_offsets,found)-1
   r=rows[i];delta=found-r['byte_offset'];line=delta//r['line_width'];col=delta%r['line_width'];center=line*r['line_bases']+col
   if col>=r['line_bases'] or center>=r['length']:
    cursor=found+1;continue
   lo=max(0,center-10000);hi=min(r['length'],center+10000)
   for kind,strand in sorted(wanted):
    if add_one(kind,i,lo,hi,strand,kind):wanted.remove((kind,strand))
   clusters+=1
   next_resolved=re.compile(b'[ACGT]').search(mapped,found,limit)
   cursor=(next_resolved.start()+1) if next_resolved else limit
if wanted:print('N examples not found in bounded discovery scan: '+repr(sorted(wanted)))

Path(a.out).write_text(json.dumps({'fai':a.fai,'reference':'HG00735_pat imported plain FASTA','text_length':positions[-1]+1,'n_discovery_max_bytes':1000000000,'n_discovery_clusters':clusters,'boundary_index':boundary,'cases':cases},indent=2)+'\n');print(json.dumps({'cases':len(cases),'kinds':sorted(set(x['kind'] for x in cases)),'labels':sorted(set(x['label'] for x in cases))}))
