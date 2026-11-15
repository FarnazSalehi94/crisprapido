import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).with_name('wave_guard.py')
spec = importlib.util.spec_from_file_location('wave_guard', SCRIPT)
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)


class WaveGuardTests(unittest.TestCase):
    def run_guard(self, root, code):
        return subprocess.run([sys.executable, str(SCRIPT), '--wave-dir', str(root),
            '--scratch', str(root), '--storage', str(root), '--scratch-floor', '0',
            '--storage-floor', '0', '--memory-floor', '0', '--interval', '.01',
            '--', sys.executable, '-c', code], capture_output=True, text=True, timeout=5)

    def test_success(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.assertEqual(self.run_guard(root, 'pass').returncode, 0)
            self.assertFalse((root/'STOP').exists())

    def test_failure_blocks_next_task(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self.assertNotEqual(self.run_guard(root, 'raise SystemExit(7)').returncode, 0)
            self.assertIn('child exit 7', (root/'STOP').read_text())
            marker = root/'should-not-exist'
            self.assertNotEqual(self.run_guard(root, f'open({str(marker)!r}, "w").close()').returncode, 0)
            self.assertFalse(marker.exists())

    def test_running_task_stops_on_latch(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            code = f'import pathlib,time; pathlib.Path({str(root/"STOP")!r}).write_text("other task failed"); time.sleep(30)'
            self.assertNotEqual(self.run_guard(root, code).returncode, 0)

    def test_checksum_failure_latches_stop(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root/'payload').write_text('changed')
            (root/'pins').write_text('0'*64 + '  payload\n')
            code = f'import subprocess; subprocess.run(["sha256sum", "--check", "pins"], cwd={d!r}, check=True)'
            self.assertNotEqual(self.run_guard(root, code).returncode, 0)
            self.assertTrue((root/'STOP').exists())
            self.assertNotEqual(self.run_guard(root, 'pass').returncode, 0)

    def test_subreaper_failure_latches_without_launch(self):
        import signal
        with tempfile.TemporaryDirectory() as d:
            old = {s: signal.getsignal(s) for s in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)}
            try:
                with patch.object(sys, 'argv', ['guard', '--wave-dir', d, '--', 'true']), patch.object(guard, 'become_subreaper', side_effect=RuntimeError('prctl failure')):
                    with patch.object(guard.subprocess, 'Popen') as launch:
                        with self.assertRaises(SystemExit): guard.main()
                        launch.assert_not_called()
                        self.assertIn('prctl failure', (Path(d)/'STOP').read_text())
            finally:
                for sig, handler in old.items(): signal.signal(sig, handler)

    def test_latch_failure_still_cleans_up(self):
        proc = object()
        with patch.object(guard, 'latch', side_effect=OSError('storage unavailable')):
            with patch.object(guard, 'cleanup') as clean:
                guard.fail_closed(Path('/unused'), 'failure', proc)
                clean.assert_called_once_with(proc)

    def test_paused_supervisor_separate_session_worker(self):
        # Run subreaping in an isolated process, not in the test runner.
        code = '''import importlib.util,sys,subprocess,os,signal,time
spec=importlib.util.spec_from_file_location('guard',sys.argv[1])
g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
g.become_subreaper()
worker="import subprocess,sys,time; p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)'],start_new_session=True); print(p.pid,flush=True); time.sleep(60)"
p=subprocess.Popen([sys.executable,'-c',worker],start_new_session=True,stdout=subprocess.PIPE,text=True)
child=int(p.stdout.readline())
os.kill(p.pid,signal.SIGSTOP)
g.cleanup(p,timeout=.05)
assert not g.descendants(),g.descendants()
try: os.kill(child,0)
except ProcessLookupError: pass
else: raise AssertionError('separate-session worker survived')
'''
        result = subprocess.run([sys.executable, '-c', code, str(SCRIPT)], capture_output=True, text=True, timeout=5)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_signal_during_preflight_prevents_launch(self):
        import os
        import signal
        with tempfile.TemporaryDirectory() as d:
            argv = ['guard', '--wave-dir', d, '--', sys.executable, '-c', 'pass']
            old = {s: signal.getsignal(s) for s in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)}
            try:
                def preflight(*args):
                    os.kill(os.getpid(), signal.SIGTERM)
                    return None
                with patch.object(sys, 'argv', argv), patch.object(guard, 'become_subreaper'):
                    with patch.object(guard, 'resource_error', side_effect=preflight), patch.object(guard.subprocess, 'Popen') as launch:
                        with self.assertRaises(SystemExit): guard.main()
                        launch.assert_not_called()
            finally:
                for sig, handler in old.items(): signal.signal(sig, handler)

    def test_signal_at_fast_child_exit_is_not_success(self):
        import os
        import signal
        with tempfile.TemporaryDirectory() as d:
            old = {s: signal.getsignal(s) for s in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)}
            try:
                with patch.object(sys, 'argv', ['guard', '--wave-dir', d, '--', 'true']), patch.object(guard, 'become_subreaper'):
                    with patch.object(guard, 'resource_error', return_value=None), patch.object(guard.subprocess, 'Popen') as launch, patch.object(guard, 'cleanup'):
                        def exited():
                            os.kill(os.getpid(), signal.SIGTERM)
                            return 0
                        launch.return_value.poll.side_effect = exited
                        launch.return_value.returncode = 0
                        with self.assertRaises(SystemExit): guard.main()
                        self.assertTrue((Path(d)/'STOP').exists())
            finally:
                for sig, handler in old.items(): signal.signal(sig, handler)

    def test_resource_floors(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            with patch.object(guard.shutil, 'disk_usage') as disk:
                disk.return_value.free = 10
                self.assertIn('free-space', guard.resource_error(root/'STOP', root, root, 11, 0, 0))
                with patch.object(guard, 'available_bytes', return_value=0):
                    self.assertIn('memory', guard.resource_error(root/'STOP', root, root, 0, 0, 1))


if __name__ == '__main__':
    unittest.main()
