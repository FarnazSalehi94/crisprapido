#!/usr/bin/env python3
"""Generate a Rust test-only module for directly probing missed oracle spans."""
import argparse
import json
from pathlib import Path
from coverage_oracle import rc
from resolve_bulge_discrepancies import fasta
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--panels',required=True);p.add_argument('--output',required=True);a=p.parse_args()
entries=[]
for panel in ('synthetic','real'):
 d=Path(a.panels)/panel; refs=fasta(d/'ref.fa');guides=fasta(d/'guides.fa')
 for row in json.loads((d/'missing_probes.json').read_text()):
  name,contig,start,end,strand=row['key'];seq=refs[contig][start:end]
  if strand=='-':seq=rc(seq)
  ident='|'.join(map(str,(panel,name,contig,start,end,strand)))
  entries.append((ident,guides[name],seq))
lines=['#[cfg(test)]','mod broad_wfa_diagnostics {','use super::*;','#[test]','fn enumerate_frozen_wfa_results() {','for (name, guide, target) in [']
for entry in entries:lines.append('({}, {}, {}),'.format(*(json.dumps(x) for x in entry)))
lines+= ['] {','let mut aligner = AffineWavefronts::with_penalties(0,3,5,1);',
          'aligner.set_alignment_span(AlignmentSpan::End2End); aligner.align(target.as_bytes(), guide.as_bytes());',
          'let cigar = String::from_utf8_lossy(aligner.cigar()).to_string(); let v = validate_cigar_accounting(&cigar, guide.as_bytes(), false);',
          'println!("WFA_RAW\\t{}\\t{}\\t{}\\t{}\\t{}\\t{}", name, cigar, v.matches, v.mismatches, v.gaps, v.max_gap_size);',
          'let result = anchored_candidate_alignment(&mut aligner, guide.as_bytes(), target.as_bytes(), 4,1,1,0.75,false);',
          'println!("WFA_PROBE\\t{}\\t{:?}", name, result);','}','}','}']
Path(a.output).write_text('\n'.join(lines)+'\n')
