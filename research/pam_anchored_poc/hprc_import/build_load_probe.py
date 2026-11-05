#!/usr/bin/env python3
"""Build a no-query native load-only probe using the prototype's Columba flags."""
import argparse,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--scratch',required=True);p.add_argument('--columba-dir',required=True);p.add_argument('--columba-build',required=True);p.add_argument('--cxx',default='g++');a=p.parse_args();root=Path(__file__).resolve().parents[3];out=Path(a.scratch);out.mkdir(parents=True,exist_ok=True)
prototype=(root/'research/pam_anchored_poc/poc.cpp').read_text();source=prototype[:prototype.index('struct Seq {')]+'''int main(int argc,char** argv) { try {
 if(argc!=3) {cerr<<"usage: load_probe INDEX REF.fa\\n";return 2;}
 require_index_binding(argv[2],argv[1]);
 auto start=chrono::steady_clock::now();FMIndex idx(argv[1],0,true,4,false,0);
 cerr<<"load_only_seconds="<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<" text_size="<<idx.getTextSize()<<"\\n";
 return 0;
 } catch(const exception& e) {cerr<<e.what()<<"\\n";return 1;} }
''';src=out/'load_probe.cpp';src.write_text(source);c=Path(a.columba_dir);b=Path(a.columba_build);helper=root/'research/pam_anchored_poc/compatibility/index_manifest.py'
cmd=[a.cxx,'-std=gnu++14','-O3','-mpopcnt','-DALPHABET=5','-DTHIRTY_TWO','-DHAVE_ZLIB',f'-DINDEX_BINDING_HELPER="{helper}"',f'-DINDEX_BINDING_PYTHON="{Path(sys.executable).resolve()}"','-include','limits','-include','cstdint','-I'+str(c/'src'),'-I'+str(b/'_deps/parallel-hashmap-src'),'-I'+str(b/'_deps/fmtlib-src/include'),str(src),str(b/'libAuxiliaryClasses.a'),str(b/'_deps/fmtlib-build/libfmt.a'),'-lz','-pthread','-o',str(out/'load_probe')];subprocess.run(cmd,check=True);print(out/'load_probe')
