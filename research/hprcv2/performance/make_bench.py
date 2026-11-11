#!/usr/bin/env python3
"""Freeze the bounded performance-benchmark plan. No job submission, no full library."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from inventory import digest, library

OLD_BINARY_SHA256 = 'ee543e772ac6e956e8fc3eb88e6ea52517e5192eff9da255902c3bad6073faed'
OLD_COLUMBA_SHA256 = '7e44974c81dad4943292fa64a71e54073908b4909ce327480464d769322c40d1'
LIBRARY_SHA256 = 'a4bac2dc4c8f5ad0c37713751835b416d1669f899de1cd05ca1a0e828be48971'
REFERENCE_SHA256 = 'bae8df3c5cf8ca96d3316b59eaa63fa37b7d3b58551c8a0700f2c530dc1e4f8b'
SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--library', required=True, type=Path)
    a = p.parse_args()
    a.root = a.root.resolve()
    if a.root.name != 'performance_investigation_v1':
        raise ValueError('unexpected benchmark root name')
    if digest(a.library) != LIBRARY_SHA256:
        raise ValueError('frozen library checksum mismatch')
    repo = Path(__file__).resolve().parents[3]
    runtime = Path('/export/local/home/salehi/crisprapido-research-workflow-v1')
    bound = Path('/export/local/home/salehi/crisprapido-hprc-import-hg00735-20261105-v1')
    old_repo = Path('/moosefs/raid5/salehi/crispr2/crispr-progress/columba-crispr-benchmark')
    plan_dir = a.root/'00_plan'
    tooling = plan_dir/'bench-tooling-v1'
    tooling.mkdir(parents=True)
    for name in ['bench_worker.py', 'validate_old.py', 'test_performance.py']:
        shutil.copyfile(Path(__file__).parent/name, tooling/name)
    for name in ['compare_sbds.py', 'collect.py', 'inventory.py']:
        shutil.copyfile(Path(__file__).parent.parent/name, tooling/name)
    shutil.copyfile(repo/'research/pam_anchored_poc/monitored_pilot/supervise.py', tooling/'supervise.py')
    records = library(a.library)
    if len(records) != 77441:
        raise ValueError('not the frozen Brunello library')
    guides = {}
    for size in [1, 100, 1000]:
        path = plan_dir/('bench-guides-N%d.fa' % size)
        with path.open('x') as f:
            for name, sequence in records[:size]:
                f.write('>'+name+'\n'+sequence+'\n')
        guides[str(size)] = str(path)
    import_json = json.loads((bound/'import.json').read_text())
    components = {s: import_json['components'][s]['sha256'] for s in SUFFIXES}
    old = dict(binary=str(old_repo/'crisprapido/target/release/crisprapido_newCFD'),
               binary_sha256=OLD_BINARY_SHA256,
               columba=str(old_repo/'columba/build_Vanilla/columba'),
               columba_sha256=OLD_COLUMBA_SHA256,
               loader='/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib/ld-linux-x86-64.so.2',
               library_path='/gnu/store/ln6hxqjvz6m9gdd9s97pivlqck7hzs99-glibc-2.35/lib:/gnu/store/gmv6n5vy5qcsn71pkapg2hnknyn1p7g3-gcc-12.3.0-lib/lib:/gnu/store/ygm2rzg08w0w96h9587k3n2yvhndcqkn-zlib-1.2.13/lib',
               time_bin='/home/salehi/miniconda3/bin/time',
               index_prefix='/moosefs/salehi/HPRCv2_CRISPRapido/indices/HG00735_pat/HG00735_pat',
               components=components, threads=16, cpus=16, mem_request='128G')
    new = dict(cli=str(runtime/'target-final/release/crisprapido'),
               native=str(runtime/'pam-search'),
               index=str(bound/'index'),
               binding=str(bound/'index.crispr-binding.json'),
               samtools='/home/salehi/miniconda3/envs/assembly_eval/bin/samtools',
               threads=1, cpus=8, mem_request='64G')
    plan = dict(format='hprcv2-performance-bench-v1', root=str(a.root), node='octopus02',
                assembly='HG00735_pat',
                reference=str(bound/'reference.fa'), reference_sha256=REFERENCE_SHA256,
                fai='/moosefs/pangenomes/HPRCv2/HG00735_pat_hprc_r2_v1.0.1.fa.gz.fai',
                fai_sha256='7f6187771a76e8af51b983e2df4e00970dd72722dea0b9bd5c425ef817e3c26d',
                settings=dict(m=2, b=1, z=1, f=0.75, pam='GG'), sizes=[1, 100, 1000],
                library=dict(path=str(a.library), sha256=LIBRARY_SHA256, records=len(records)),
                new=new, old=old, guides=guides,
                matrices=dict(mismatch=str(repo/'mismatch_scores.txt'), pam=str(repo/'pam_scores.txt')),
                release_scripts=str(repo/'scripts/release_validation'),
                supervisor=str(tooling/'supervise.py'), worker=str(tooling/'bench_worker.py'),
                validator=str(tooling/'validate_old.py'),
                limits=dict(stage_timeout_seconds=3600, min_start_available_kib=67108864,
                            min_running_available_kib=8388608, min_scratch_bytes=60000000000,
                            min_running_disk_bytes=10000000000, file_limit_bytes=10737418240,
                            new_rss_kib=41943040, old_rss_kib=115343360),
                provenance=dict(git_commit=subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], universal_newlines=True).strip(),
                    scope='Bounded 1/100/1000-guide matched-settings benchmark on one verified assembly. '
                          'Historical binary output is diagnostic only (known legacy PAM/geometry semantics); '
                          'it is never a valid search result or a recall standard. No full-library launch.',
                    new_engine='Experimental PAM-anchored backend, single-threaded, strict NGG, CFD default (no neural stage).',
                    old_engine='Preserved production binary crisprapido_newCFD (Columba candidates, WFA2 verification), '
                               'production -t 16 configuration, production Guix runtime, index staged to local scratch exactly as in job 2856949.'))
    pin_paths = [Path(new[k]) for k in ['cli', 'native', 'binding', 'samtools']]
    pin_paths += [Path(old[k]) for k in ['binary', 'columba', 'loader', 'time_bin']]
    pin_paths += [Path(plan['reference']), Path(plan['fai']), Path(plan['supervisor']),
                  Path(plan['worker']), Path(plan['validator']), Path(tooling/'compare_sbds.py'),
                  Path(tooling/'collect.py'), Path(tooling/'inventory.py'),
                  Path(plan['matrices']['mismatch']), Path(plan['matrices']['pam']),
                  Path(plan['release_scripts'])/'audit_emitted.py', a.library]
    pin_paths += [Path(p) for p in guides.values()]
    plan['pins'] = {str(p): digest(p) for p in sorted(set(pin_paths))}
    for model in json.loads((runtime/'source/research/crispr_bulge/model_manifest.json').read_text())['models']:
        pass  # No neural inference in this benchmark; models deliberately not involved.
    with (plan_dir/'bench-plan.json').open('x') as f:
        json.dump(plan, f, indent=2)
        f.write('\n')
    print(json.dumps(dict(plan=str(plan_dir/'bench-plan.json'), sha256=digest(plan_dir/'bench-plan.json'))))


if __name__ == '__main__':
    main()
