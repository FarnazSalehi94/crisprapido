#!/usr/bin/env python3
"""Build a measurement-only native executable; does not edit the prototype."""
import argparse,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out-dir',required=True);p.add_argument('--columba-dir',required=True);p.add_argument('--columba-build',required=True);p.add_argument('--cxx',default='g++');a=p.parse_args();root=Path(__file__).resolve().parents[3];out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True);src=(root/'research/pam_anchored_poc/poc.cpp').read_text()
def replace(old,new):
 global src
 assert src.count(old)==1,old;src=src.replace(old,new)
replace('struct Search {', '''static long long emit_ns=0, replay_ns=0, emit_calls=0, replay_calls=0;
struct ScopeTimer { long long& ns; chrono::steady_clock::time_point start; ScopeTimer(long long& total,long long& calls):ns(total),start(chrono::steady_clock::now()){++calls;} ~ScopeTimer(){ns+=chrono::duration_cast<chrono::nanoseconds>(chrono::steady_clock::now()-start).count();} };
struct Search {''')
replace('   reverse(ops.begin(),ops.end());','   ScopeTimer emit_timer(emit_ns,emit_calls);\n   reverse(ops.begin(),ops.end());')
replace('      // The indexed letter at a reference N is a replacement, so recompute','      ScopeTimer replay_timer(replay_ns,replay_calls);\n      // The indexed letter at a reference N is a replacement, so recompute')
replace(' require_index_binding(argv[2],argv[1]);\n auto ref=fasta(argv[2]), guides=fasta(argv[3]);',' auto verify_start=chrono::steady_clock::now(); require_index_binding(argv[2],argv[1]); cerr<<"PROFILE native_binding_seconds="<<chrono::duration<double>(chrono::steady_clock::now()-verify_start).count()<<"\\n";\n auto ref_start=chrono::steady_clock::now(); auto ref=fasta(argv[2]); cerr<<"PROFILE native_reference_load_seconds="<<chrono::duration<double>(chrono::steady_clock::now()-ref_start).count()<<"\\n"; auto guides=fasta(argv[3]);')
replace(' FMIndex idx(argv[1],0,true,4,false,0); long long states=0; double search_seconds=0;',' auto index_start=chrono::steady_clock::now(); FMIndex idx(argv[1],0,true,4,false,0); cerr<<"PROFILE native_index_load_seconds="<<chrono::duration<double>(chrono::steady_clock::now()-index_start).count()<<"\\n"; long long states=0; double search_seconds=0;')
replace(' cerr<<"states="<<states<<" search_wall_seconds="<<search_seconds', ' cerr<<"PROFILE emit_seconds="<<emit_ns/1e9<<" emit_calls="<<emit_calls<<" replay_seconds="<<replay_ns/1e9<<" replay_calls="<<replay_calls<<"\\n";\n cerr<<"states="<<states<<" search_wall_seconds="<<search_seconds')
cpp=out/'profile_native.cpp';cpp.write_text(src);c=Path(a.columba_dir);b=Path(a.columba_build);helper=root/'research/pam_anchored_poc/compatibility/index_manifest.py'
cmd=[a.cxx,'-std=gnu++14','-O3','-mpopcnt','-DALPHABET=5','-DTHIRTY_TWO','-DHAVE_ZLIB',f'-DINDEX_BINDING_HELPER="{helper}"',f'-DINDEX_BINDING_PYTHON="{Path(sys.executable).resolve()}"','-include','limits','-include','cstdint','-I'+str(c/'src'),'-I'+str(b/'_deps/parallel-hashmap-src'),'-I'+str(b/'_deps/fmtlib-src/include'),str(cpp),str(b/'libAuxiliaryClasses.a'),str(b/'_deps/fmtlib-build/libfmt.a'),'-lz','-pthread','-o',str(out/'profile_native')];subprocess.run(cmd,check=True);print(out/'profile_native')
