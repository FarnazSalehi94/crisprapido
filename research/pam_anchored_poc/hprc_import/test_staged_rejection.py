#!/usr/bin/env python3
"""Test an imported binding with separate fake paths; never mutate the import."""
import argparse,os,sys,tempfile
from pathlib import Path
root=Path(__file__).resolve().parents[3];sys.path.insert(0,str(root/'research/pam_anchored_poc/compatibility'));from index_manifest import SUFFIXES,verify
p=argparse.ArgumentParser();p.add_argument('--stage',required=True);a=p.parse_args();stage=Path(a.stage).resolve()
with tempfile.TemporaryDirectory(prefix='staged-binding-negative-') as td:
 pfx=Path(td)/'index'
 for s in SUFFIXES:os.symlink(stage/('index.'+s),str(pfx)+'.'+s)
 os.symlink(stage/'index.crispr-binding.json',str(pfx)+'.crispr-binding.json')
 ref=Path(td)/'mismatched.fa';ref.write_text('>bad\nACGT\n')
 try:verify(ref,pfx)
 except ValueError as e:assert 'reference FASTA does not match' in str(e)
 else:raise AssertionError('mismatched reference accepted')
 # Restore the exact staged reference and replace the first component by a
 # truncated stand-in. It fails before the other large components are read.
 bwt=Path(str(pfx)+'.bwt');bwt.unlink();bwt.write_bytes(b'corrupt')
 try:verify(stage/'reference.fa',pfx)
 except ValueError as e:assert 'index component does not match binding: bwt' in str(e)
 else:raise AssertionError('corrupt index accepted')
 print('published reference and component mismatches rejected; stage unchanged')
