#!/usr/bin/env python3
"""Build only the isolated research executable into a NEW external directory."""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--out', type=Path, required=True)
p.add_argument('--columba-dir', type=Path, required=True)
p.add_argument('--columba-build', type=Path, required=True)
p.add_argument('--cxx', required=True)
a = p.parse_args()
root = Path(__file__).resolve().parents[3]
source = root/'research/pam_anchored_poc/poc.cpp'
expected = subprocess.check_output(['git','show','252a6bb:research/pam_anchored_poc/poc.cpp'], cwd=root)
assert source.read_bytes() == expected, 'review changes to unpruned source before rebuilding'
a.out.mkdir(exist_ok=False, parents=True)
ldlib = '/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib'
gcclib = '/gnu/store/gmv6n5vy5qcsn71pkapg2hnknyn1p7g3-gcc-12.3.0-lib/lib'
cmd = [a.cxx, '-std=gnu++14', '-O3', '-mpopcnt', '-fno-fast-math', '-ffp-contract=off',
       '-DALPHABET=5', '-DTHIRTY_TWO', '-DHAVE_ZLIB',
       f'-DINDEX_BINDING_HELPER="{root}/research/pam_anchored_poc/compatibility/index_manifest.py"',
       f'-DINDEX_BINDING_PYTHON="{Path(sys.executable).resolve()}"',
       '-include', 'limits', '-include', 'cstdint', '-I'+str(a.columba_dir/'src'),
       '-I'+str(a.columba_build/'_deps/parallel-hashmap-src'),
       '-I'+str(a.columba_build/'_deps/fmtlib-src/include'),
       str(Path(__file__).with_name('native.cpp')), str(a.columba_build/'libAuxiliaryClasses.a'),
       str(a.columba_build/'_deps/fmtlib-build/libfmt.a'), '-lz', '-pthread',
       '-Wl,--dynamic-linker='+ldlib+'/ld-linux-x86-64.so.2', '-Wl,-rpath,'+ldlib+':'+gcclib,
       '-o', str(a.out/'native')]
subprocess.run(cmd, check=True)
hashfile = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
files = [source, Path(__file__).with_name('native.cpp'), a.out/'native',
         a.columba_build/'libAuxiliaryClasses.a', a.columba_build/'_deps/fmtlib-build/libfmt.a']
(a.out/'build.json').write_text(json.dumps(dict(command=cmd, sha256={str(p):hashfile(p) for p in files},
    compiler=subprocess.check_output([a.cxx,'--version'],text=True).splitlines()[0]), indent=2)+'\n')
print(a.out/'native')
