#!/usr/bin/env python3
"""Produce disposable timing sources; no production source/binary is modified."""
import argparse,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--out',required=True);a=p.parse_args();repo=Path(a.repo).resolve();out=Path(a.out).resolve();out.mkdir(exist_ok=False)
subprocess.run(['python3',str(repo/'scripts/release_validation/make_profile_build.py'),'--no-cache','--source',str(repo),'--output',str(out/'baseline')],check=True)
f=out/'baseline/src/columba.rs';s=f.read_text();old='pub(crate) fn parse_columba_sam_file(path: &PathBuf) -> Result<Vec<ColumbaCandidate>, String> {';assert s.count(old)==1;s=s.replace(old,old+'\n    let _p = crate::profile::guard(11);');f.write_text(s)
f=out/'baseline/src/profile.rs';s=f.read_text().replace('AtomicU64; 11','AtomicU64; 12').replace('}; 11]', '}; 12]').replace('&str; 11','&str; 12').replace('"cache_miss"]','"cache_miss", "sam_read_parse"]');f.write_text(s)
f=repo/'research/pam_anchored_poc/poc.cpp';s=f.read_text();s=s.replace('struct Search {','''struct EmitTimer {
 long long& ns; chrono::steady_clock::time_point t;
 explicit EmitTimer(long long& value):ns(value),t(chrono::steady_clock::now()) {}
 ~EmitTimer(){ns+=chrono::duration_cast<chrono::nanoseconds>(chrono::steady_clock::now()-t).count();}
};
struct Search {''');s=s.replace('map<tuple<string,char,size_t>,Hit> hits;','map<tuple<string,char,size_t>,Hit> hits; long long emit_ns=0;');sig='void emit(const SARangePair& range, int span, int /* indexed_mm */, int ins, int del, string ops) {';assert s.count(sig)==1;s=s.replace(sig,sig+'\n   EmitTimer timer(emit_ns);');s=s.replace('long long states=0; double search_seconds=0;', 'long long states=0, emit_ns=0; double search_seconds=0;').replace('states+=s.states; for(auto& kv:s.hits)', 'states+=s.states; emit_ns+=s.emit_ns; for(auto& kv:s.hits)').replace('<<search_seconds<<\'\\n\';','<<search_seconds<<" locate_replay_ns="<<emit_ns<<\'\\n\';');(out/'poc.cpp').write_text(s)
