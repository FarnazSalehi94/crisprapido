#!/usr/bin/env python3
"""Exercise bounded probe and in-memory N controls on a small bound index."""
import argparse,csv,struct,subprocess,sys,tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[3];sys.path.insert(0,str(root/'research/pam_anchored_poc/validation'));from index_check import fasta
p=argparse.ArgumentParser();p.add_argument('--probe',required=True);p.add_argument('--index',required=True);p.add_argument('--reference',required=True);p.add_argument('--guides',required=True);p.add_argument('--native-tsv',required=True);a=p.parse_args();refs=dict(fasta(a.reference));guides=dict(fasta(a.guides));names=list(refs);raw=Path(a.index+'.pos').read_bytes();offsets=struct.unpack('<'+'I'*(len(raw)//4),raw);seq=Path(a.index+'.txt.bin').read_bytes()[4:];lines=Path(a.native_tsv).read_text().splitlines();i=next(i for i,s in enumerate(lines) if s.startswith('guide\tcontig\t'));rows=list(csv.DictReader(lines[i:],delimiter='\t'))
with tempfile.TemporaryDirectory(prefix='exact-probe-fixture-') as td:
 cases=Path(td)/'cases.tsv';selected=[]
 for strand in '+-':
  r=next(x for x in rows if x['strand']==strand and x['cigar']=='20=' and x['mismatches']=='0');ci=names.index(r['contig']);start=int(r['start']) if strand=='+' else int(r['start'])-3;genomic=refs[r['contig']][start:start+23];indexed=seq[offsets[ci]+start:offsets[ci]+start+23].decode('ascii');assert genomic==indexed
  selected.append(['before_2_31','resolved',strand,r['contig'],offsets[ci]+start,start,r['anchor'],indexed,genomic,guides[r['guide']],0])
 cases.write_text(''.join('\t'.join(map(str,x))+'\n' for x in selected));run=subprocess.run([a.probe,a.index,a.reference,str(cases)],capture_output=True,text=True,timeout=120,check=True,cwd=root);case=[s for s in run.stdout.splitlines() if s.startswith('CASE\t')];controls=[s for s in run.stdout.splitlines() if s.startswith('CONTROL\t')];assert len(case)==2 and len(controls)==6;print('two-strand exact fixture and six N replay/PAM controls passed')
