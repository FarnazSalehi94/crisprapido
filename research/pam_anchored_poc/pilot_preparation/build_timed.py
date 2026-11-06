#!/usr/bin/env python3
"""Build a separate native prototype with an index-constructor timer."""
import argparse,os,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--scratch',required=True);p.add_argument('--columba-dir',required=True);p.add_argument('--columba-build',required=True);p.add_argument('--cxx',default='g++');a=p.parse_args()
root=Path(__file__).resolve().parents[3];scratch=Path(a.scratch);scratch.mkdir(parents=True,exist_ok=True);src=(root/'research/pam_anchored_poc/poc.cpp').read_text()
old=' FMIndex idx(argv[1],0,true,4,false,0); long long states=0; double search_seconds=0;'
new=' auto index_start=chrono::steady_clock::now(); FMIndex idx(argv[1],0,true,4,false,0); cerr << "index_load_seconds=" << chrono::duration<double>(chrono::steady_clock::now()-index_start).count() << endl; long long states=0; double search_seconds=0;'
assert src.count(old)==1;src=src.replace(old,new);(scratch/'poc_timed.cpp').write_text(src)
helper=root/'research/pam_anchored_poc/compatibility/index_manifest.py';python=Path(sys.executable).resolve();build=Path(a.columba_build);source=Path(a.columba_dir)
cmd=[a.cxx,'-std=gnu++14','-O3','-mpopcnt','-DALPHABET=5','-DTHIRTY_TWO','-DHAVE_ZLIB',f'-DINDEX_BINDING_HELPER="{helper}"',f'-DINDEX_BINDING_PYTHON="{python}"','-include','limits','-include','cstdint','-I'+str(source/'src'),'-I'+str(build/'_deps/parallel-hashmap-src'),'-I'+str(build/'_deps/fmtlib-src/include'),str(scratch/'poc_timed.cpp'),str(build/'libAuxiliaryClasses.a'),str(build/'_deps/fmtlib-build/libfmt.a'),'-lz','-pthread','-o',str(scratch/'poc_timed')]
subprocess.run(cmd,check=True);print(scratch/'poc_timed')
