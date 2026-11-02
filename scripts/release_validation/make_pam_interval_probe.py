#!/usr/bin/env python3
"""Create disposable profiled source copy with a strict-PAM interval precheck.

Research prototype only. Never modifies the checkout or production binary.
Uses already instrumented source from make_profile_build.py --no-cache.
"""
import argparse,shutil
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True);p.add_argument('--output',required=True);a=p.parse_args()
src=Path(a.source);out=Path(a.output);shutil.copytree(src,out)
f=out/'src/verification.rs';s=f.read_text();needle='    let candidate_seq = seq.get(reference_start..reference_end)?;'
if s.count(needle)!=1:raise ValueError('source layout differs')
insert='''    // Experimental: only the strict default mode; never discard a candidate,
    // and preserve the existing overlapping-interval expansion policy.
    if args.max_bulges == 1 && args.max_bulge_size == 2
        && !args.no_filter && !args.allow_noncanonical_pam
    {
        let strand = if reverse { '-' } else { '+' };
        if !extract_adjacent_pam(seq, reference_start, reference_end, strand)
            .as_deref().map_or(false, |pam| pam.eq_ignore_ascii_case(&args.pam))
        {
            return None;
        }
    }
'''
f.write_text(s.replace(needle,insert+needle,1))
