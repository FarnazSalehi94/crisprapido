#!/usr/bin/env python3
"""Instrument a separate f7de425 source snapshot's runner with phase timers."""
import argparse,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--source',required=True);a=p.parse_args();root=Path(__file__).resolve().parents[3];snapshot=Path(a.source);path=snapshot/'research/pam_anchored_poc/integration/run_snapshot.py';original=subprocess.run(['git','show','f7de425:research/pam_anchored_poc/integration/run_snapshot.py'],cwd=root,text=True,capture_output=True,check=True).stdout
if path.read_text()!=original:raise SystemExit('archived runner differs from f7de425; refusing modification')
s=original
def rep(old,new):
 global s
 assert s.count(old)==1,old;s=s.replace(old,new)
rep('import tempfile\n','import tempfile\nimport time\n')
rep('def main():\n', 'def main():\n    runner_start = time.monotonic()\n')
rep('        shutil.copyfile(a.reference, ref)','        phase_start = time.monotonic()\n        shutil.copyfile(a.reference, ref)')
rep('        guides.write_text(sys.stdin.read(), encoding=\'ascii\')','        print("PROFILE snapshot_copy_seconds=%.6f" % (time.monotonic()-phase_start), file=sys.stderr)\n        guides.write_text(sys.stdin.read(), encoding=\'ascii\')')
rep('        verify(ref, idx)\n        for path in root.iterdir():','        phase_start = time.monotonic()\n        verify(ref, idx)\n        print("PROFILE snapshot_verify_pre_seconds=%.6f" % (time.monotonic()-phase_start), file=sys.stderr)\n        for path in root.iterdir():')
rep('        run = subprocess.run([a.binary,','        phase_start = time.monotonic()\n        run = subprocess.run([a.binary,')
rep('        sys.stderr.write(run.stderr)','        print("PROFILE native_process_seconds=%.6f" % (time.monotonic()-phase_start), file=sys.stderr)\n        sys.stderr.write(run.stderr)')
rep('        verify(ref, idx)\n        # Extract targets', '        phase_start = time.monotonic()\n        verify(ref, idx)\n        print("PROFILE snapshot_verify_post_seconds=%.6f" % (time.monotonic()-phase_start), file=sys.stderr)\n        # Extract targets')
rep('        seqs = {}\n        name = None','        phase_start = time.monotonic()\n        seqs = {}\n        name = None')
rep("        lines = run.stdout.splitlines()", "        print(\"PROFILE python_reference_load_seconds=%.6f\" % (time.monotonic()-phase_start), file=sys.stderr)\n        lines = run.stdout.splitlines()")
rep("            print(line + '\\t' + target + '\\t' + str(len(seq)))\n", "            print(line + '\\t' + target + '\\t' + str(len(seq)))\n        print(\"PROFILE runner_total_seconds=%.6f\" % (time.monotonic()-runner_start), file=sys.stderr)\n")
backup=path.with_name('run_snapshot.uninstrumented.py');backup.write_text(original);path.write_text(s);print(path)
