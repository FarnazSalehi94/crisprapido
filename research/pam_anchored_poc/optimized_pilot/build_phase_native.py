#!/usr/bin/env python3
"""Build a measurement-only native executable; does not edit the prototype."""
import argparse,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out-dir',required=True);p.add_argument('--columba-dir',required=True);p.add_argument('--columba-build',required=True);p.add_argument('--cxx',default='g++');a=p.parse_args();root=Path(__file__).resolve().parents[3];out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True);src=(root/'research/pam_anchored_poc/poc.cpp').read_text();expected=subprocess.run(['git','show','f7de425:research/pam_anchored_poc/poc.cpp'],cwd=root,text=True,capture_output=True,check=True).stdout
if src!=expected:raise SystemExit('native source differs from f7de425')
def replace(old,new):
 global src
 assert src.count(old)==1,old;src=src.replace(old,new)
replace(' require_index_binding(argv[2],argv[1]);\n auto ref=fasta(argv[2]), guides=fasta(argv[3]);',' auto verify_start=chrono::steady_clock::now(); require_index_binding(argv[2],argv[1]); cerr<<"PROFILE native_binding_seconds="<<chrono::duration<double>(chrono::steady_clock::now()-verify_start).count()<<"\\n";\n auto ref_start=chrono::steady_clock::now(); auto ref=fasta(argv[2]); cerr<<"PROFILE native_reference_load_seconds="<<chrono::duration<double>(chrono::steady_clock::now()-ref_start).count()<<"\\n"; auto guides=fasta(argv[3]);')
replace(' FMIndex idx(argv[1],0,true,4,false,0); long long states=0; double search_seconds=0;',' auto index_start=chrono::steady_clock::now(); FMIndex idx(argv[1],0,true,4,false,0); cerr<<"PROFILE native_index_load_seconds="<<chrono::duration<double>(chrono::steady_clock::now()-index_start).count()<<"\\n"; long long states=0; double search_seconds=0;')
cpp=out/'profile_native.cpp';cpp.write_text(src);c=Path(a.columba_dir);b=Path(a.columba_build);helper=root/'research/pam_anchored_poc/compatibility/index_manifest.py'
cmd=[a.cxx,'-std=gnu++14','-O3','-mpopcnt','-DALPHABET=5','-DTHIRTY_TWO','-DHAVE_ZLIB',f'-DINDEX_BINDING_HELPER="{helper}"',f'-DINDEX_BINDING_PYTHON="{Path(sys.executable).resolve()}"','-include','limits','-include','cstdint','-I'+str(c/'src'),'-I'+str(b/'_deps/parallel-hashmap-src'),'-I'+str(b/'_deps/fmtlib-src/include'),str(cpp),str(b/'libAuxiliaryClasses.a'),str(b/'_deps/fmtlib-build/libfmt.a'),'-lz','-pthread','-o',str(out/'profile_native')];subprocess.run(cmd,check=True);print(out/'profile_native')
