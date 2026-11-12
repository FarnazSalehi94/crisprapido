#!/usr/bin/env python3
"""Stage-1 deep verification wrapper for one unverified SBDS assembly.

Runs the repaired prepare_binding in DEEP mode: full component reads produce
new trusted checksums; nothing is compared to (absent) historical evidence.
Never rebuilds or overwrites an original index: a failed verification records
FAILED evidence (EIO labelled separately), preserves the partial stage, and
leaves originals untouched. Committed bindings are re-checked for reusability.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path('/moosefs/salehi/HPRCv2_CRISPRapido/runs/sbds_ngg_nrr_v1')
SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()
PREPARE = ROOT/'00_manifest/prepare_binding-v1.py'


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def main():
    assembly = sys.argv[1]
    cohort = {l.split('\t')[0]: l.split('\t')[1] for l in
              Path('/moosefs/salehi/HPRCv2_CRISPRapido/analysis_genomes.tsv').read_text().splitlines() if l.strip()}
    fasta_gz = cohort[assembly]
    index_prefix = '/moosefs/salehi/HPRCv2_CRISPRapido/indices/{}/{}'.format(assembly, assembly)
    records = ROOT/'00_manifest/verification-batch-1'
    records.mkdir(exist_ok=True)
    record = records/(assembly+'.json')
    if record.exists():
        previous = json.loads(record.read_text())
        if previous['status'] == 'BOUND':
            # Reusability check: all committed binding files unchanged.
            stage = ROOT/'00_manifest/bindings'/assembly
            binding = json.loads((stage/'binding.json').read_text())
            ok = (digest(stage/'reference.fa') == binding['reference_sha256']
                  and (stage/'INCOMPLETE').exists() is False)
            print(json.dumps(dict(assembly=assembly, status='RECHECKED', reusable=ok)))
            return
    stage = ROOT/'00_manifest/bindings'/assembly
    component_bytes = sum(os.path.getsize(index_prefix+'.'+s) for s in SUFFIXES)
    read_bytes = component_bytes + os.path.getsize(fasta_gz) + os.path.getsize(fasta_gz+'.fai')
    started = time.monotonic()
    result = subprocess.run([sys.executable, str(PREPARE), '--assembly', assembly,
        '--index-prefix', index_prefix, '--fasta-gz', fasta_gz, '--fai', fasta_gz+'.fai',
        '--out', str(stage),
        '--scoring-python', '/usr/bin/python3',
        '--reference-index-helper', '/export/local/home/salehi/crisprapido-pam-validation/research/crispr_bulge/reference_index.py'],
        capture_output=True, text=True)
    elapsed = time.monotonic()-started
    if result.returncode != 0:
        error = result.stderr.strip().splitlines()[-1] if result.stderr.strip() else 'unknown'
        kind = 'EIO' if ('Input/output error' in result.stderr or 'Errno 5' in result.stderr) else 'OTHER'
        data = dict(assembly=assembly, status='FAILED', failure_kind=kind, error=error,
                    attempted_stage=str(stage), elapsed_seconds=elapsed, read_bytes_planned=read_bytes,
                    note='Originals untouched; partial stage preserved; no rebuild attempted.')
        with record.open('x') as f: json.dump(data, f, indent=2); f.write('\n')
        print(json.dumps(data))
        sys.exit(1)
    binding = json.loads((stage/'binding.json').read_text())
    data = dict(assembly=assembly, status='BOUND', elapsed_seconds=elapsed,
                read_bytes=read_bytes, written_bytes=sum(f.stat().st_size for f in stage.rglob('*') if f.is_file()),
                timings=binding['timings'], reference_sha256=binding['reference_sha256'],
                components_sha256_count=len(binding.get('components_sha256', {})),
                method='deep full-content verification; new checksums produced this run',
                reusable_check='digest(reference.fa) matches binding record')
    with record.open('x') as f: json.dump(data, f, indent=2); f.write('\n')
    print(json.dumps({k: data[k] for k in ('assembly', 'status', 'elapsed_seconds', 'read_bytes')}))


if __name__ == '__main__':
    main()
