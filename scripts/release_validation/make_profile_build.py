#!/usr/bin/env python3
"""Create a disposable instrumented source copy. Does not modify the worktree."""
import argparse
from pathlib import Path
import shutil
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True);p.add_argument('--output',required=True);p.add_argument('--no-cache',action='store_true');a=p.parse_args()
src=Path(a.source);out=Path(a.output);out.mkdir(exist_ok=False);shutil.copytree(src/'src',out/'src')
for x in ('Cargo.toml','Cargo.lock'):shutil.copyfile(src/x,out/x)
def replace(path,old,new):
 f=out/'src'/path;s=f.read_text();n=s.count(old)
 if n!=1:raise ValueError(f'{path}: expected one match, got {n}: {old!r}')
 f.write_text(s.replace(old,new,1))
replace('main.rs','mod cfd_score;','mod profile;\nmod cfd_score;')
replace('main.rs','fn main() {','fn main() {\n    let _profile = profile::Session::new();')
replace('columba.rs','let output = command.output().map_err(|e| {','let output = { let _p = crate::profile::guard(0); command.output() }.map_err(|e| {')
replace('verification.rs',') -> Option<String> {\n    // Cas9 has a three-base PAM;',') -> Option<String> {\n    let _p = crate::profile::guard(4);\n    // Cas9 has a three-base PAM;')
replace('verification.rs','fn hit_matches_requested_pam(hit: &Hit, requested_pam: &str) -> bool {','fn hit_matches_requested_pam(hit: &Hit, requested_pam: &str) -> bool {\n    let _p = crate::profile::guard(4);')
replace('reporting.rs','fn pam_matches_requested(hit: &Hit, requested_pam: &str) -> bool {','fn pam_matches_requested(hit: &Hit, requested_pam: &str) -> bool {\n    let _p = crate::profile::guard(4);')
if a.no_cache:
 replace('verification.rs', "    let candidate_seq = seq.get(reference_start..reference_end)?;", "    crate::profile::event(8); crate::profile::event(10);\n    let candidate_seq = seq.get(reference_start..reference_end)?;")
else:
 replace('verification.rs','''    let key = (record_id, reference_start, reference_end, reverse);
     if let Some(result) = cache.get(&key) {
         return result.clone();
     }
     let result = verify_columba_interval_uncached(''','''    crate::profile::event(8);
     let key = (record_id, reference_start, reference_end, reverse);
     if let Some(result) = cache.get(&key) {
         crate::profile::event(9);
         return result.clone();
     }
     crate::profile::event(10);
     let result = verify_columba_interval_uncached(''')
replace('verification.rs',') -> Vec<(usize, usize)> {\n    // Any target overlapping',') -> Vec<(usize, usize)> {\n    let _p = crate::profile::guard(1);\n    // Any target overlapping')
replace('verification.rs','''    aligner.align(candidate_span, guide);
    let score = aligner.score();''','''    { let _p = crate::profile::guard(2); aligner.align(candidate_span, guide); }
    let score = aligner.score();''')
replace('verification.rs',') -> Option<(i32, String, u32, u32, u32)> {\n    let (op, positions)',') -> Option<(i32, String, u32, u32, u32)> {\n    let _p = crate::profile::guard(3);\n    let (op, positions)')
replace('reporting.rs','''    let mut by_anchor: BTreeMap<(String, char, usize), Hit> = BTreeMap::new();''','''    let selection_guard = crate::profile::guard(6);
    let mut by_anchor: BTreeMap<(String, char, usize), Hit> = BTreeMap::new();''')
replace('reporting.rs','''    let mut output = String::new();
    for best_hit in by_anchor.values()''','''    drop(selection_guard);
    let mut output = String::new();
    for best_hit in by_anchor.values()''')
replace('reporting.rs','''        cfd_score::get_cfd_score(guide, target_seq, cigar, pam)''','''        let _p = crate::profile::guard(5);
        cfd_score::get_cfd_score(guide, target_seq, cigar, pam)''')
(out/'src/profile.rs').write_text('''use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;
static TIMES: [AtomicU64; 11] = [const { AtomicU64::new(0) }; 11];
static COUNTS: [AtomicU64; 11] = [const { AtomicU64::new(0) }; 11];
const NAMES: [&str; 11] = ["columba_wait", "expansion_enumeration", "wfa_align", "eligible_enumeration", "pam_extract", "cfd_score", "representative_selection", "unused", "cache_lookup", "cache_hit", "cache_miss"];
pub struct Guard(usize, Instant);
pub fn guard(i: usize) -> Guard { COUNTS[i].fetch_add(1, Ordering::Relaxed); Guard(i, Instant::now()) }
pub fn event(i: usize) { COUNTS[i].fetch_add(1, Ordering::Relaxed); }
impl Drop for Guard { fn drop(&mut self) { TIMES[self.0].fetch_add(self.1.elapsed().as_nanos().min(u64::MAX as u128) as u64, Ordering::Relaxed); } }
pub struct Session(Instant);
impl Session { pub fn new() -> Self { Session(Instant::now()) } }
impl Drop for Session { fn drop(&mut self) { for i in 0..NAMES.len() { eprintln!("PROFILE\\t{}\\t{}\\t{}",NAMES[i],COUNTS[i].load(Ordering::Relaxed),TIMES[i].load(Ordering::Relaxed)); } eprintln!("PROFILE\\ttotal_wall\\t1\\t{}",self.0.elapsed().as_nanos()); } }
''')
