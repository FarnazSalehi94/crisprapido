#!/usr/bin/env python3
"""Read-only historical audit. Metadata presence is NEVER content verification."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil

SBDS = 'ACTGAAATCTGTAAGCAGGC'
LIBRARY_SHA256 = 'a4bac2dc4c8f5ad0c37713751835b416d1669f899de1cd05ca1a0e828be48971'
SUFFIXES = 'brt bwt cct fsid headerSN.bin meta pos rev.brt sa.4 sa.bv.4 sna txt.bin'.split()
REFERENCES = {'GRCh38', 'CHM13'}


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def tsv(path, fields, rows):
    with Path(path).open('x') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def pairs(path):
    if Path(path).stat().st_size > 1024*1024:
        raise ValueError('not a small metadata table: '+str(path))
    rows = [line.split('\t', 1) for line in Path(path).read_text().splitlines() if line]
    if any(len(r) != 2 for r in rows) or len({r[0] for r in rows}) != len(rows):
        raise ValueError('invalid key/value table: '+str(path))
    return dict(rows)


def library(path):
    rows, name, seq = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if name is not None:
                rows.append((name, ''.join(seq)))
            name, seq = line[1:].strip(), []
        elif line.strip():
            if name is None:
                raise ValueError('sequence before FASTA header')
            seq.append(line.strip())
    if name is not None:
        rows.append((name, ''.join(seq)))
    if not rows or len({n for n, _ in rows}) != len(rows):
        raise ValueError('empty library or duplicate guide ID')
    if any(not n or re.search(r'\s', n) or len(s) != 20 or set(s)-set('ACGT') for n, s in rows):
        raise ValueError('invalid guide identifier/sequence; no silent normalization')
    return rows


def choose_pilot(rows):
    # Fixed sequence-derived strata, not selected for favourable measured timing.
    # A small safety pilot, NOT a statistically representative library sample.
    ranked = sorted(rows, key=lambda x: (sum(c in 'GC' for c in x[1]), digest_text(x[0])))
    selected = []
    for label, fraction in [('GC_q10', .10), ('GC_q50', .50), ('GC_q90', .90)]:
        n, s = ranked[int((len(ranked)-1)*fraction)]
        selected.append(dict(guide_id=n, sequence=s, stratum=label, gc_bases=sum(c in 'GC' for c in s)))
    if len({r['guide_id'] for r in selected}) != 3:
        raise ValueError('pilot strata overlap')
    return selected


def digest_text(text):
    return hashlib.sha256(text.encode()).hexdigest()


def historical_sbds(path):
    raw = path.read_bytes()
    stats = collections.Counter()
    anchors, perfect = set(), []
    for line in raw.decode().splitlines():
        if not line:
            continue
        f = line.split('\t')
        if len(f) < 12:
            raise ValueError('short historical PAF row: '+str(path))
        tags = {}
        for tag in f[12:]:
            k, ty, val = tag.split(':', 2)
            if k in tags:
                raise ValueError('duplicate historical tag')
            tags[k] = (ty, val)
        stats['rows'] += 1
        key = (f[0], f[5], f[4], f[8] if f[4] == '+' else f[7])
        stats['duplicate_reported_anchors'] += key in anchors
        anchors.add(key)
        if 'cf' in tags:
            score = float(tags['cf'][1])
            for threshold in (.02, .05, .10, .20):
                stats['cfd_ge_'+format(threshold, '.2f')] += score >= threshold
        if tags.get('cg') == ('Z', '20='):
            perfect.append(dict(contig=f[5], strand=f[4], start0=int(f[7]), end0=int(f[8])))
    return dict(stats, perfect_20eq=len(perfect), perfect_intervals=perfect,
                sha256=hashlib.sha256(raw).hexdigest(), bytes=len(raw))


def historical_id(fasta):
    name = Path(fasta).name
    if not name.endswith('.fa.gz'):
        raise ValueError('unexpected frozen FASTA suffix')
    return name[:-6].split('_hprc_', 1)[0]


def audit(base, out):
    base, out = Path(base), Path(out)
    out.mkdir(exist_ok=False)
    source_files = [base/'manifests/HPRCv2_all_466_fastas.txt', base/'manifests/HPRCv2_SBDS_fastas.txt',
                    base/'analysis_genomes.tsv', base/'guides/brunello_77441.fa', base/'guides/brunello_clean.tsv',
                    base/'guides/Brunello.txt', base/'guides/Brunello_unix.txt']
    sources = {str(p): dict(sha256=digest(p), bytes=p.stat().st_size) for p in source_files}
    all_fastas = source_files[0].read_text().splitlines()
    cohort_fastas = source_files[1].read_text().splitlines()
    assembly_map = pairs(base/'analysis_genomes.tsv')
    if (len(all_fastas) != 466 or len(set(all_fastas)) != 466 or len(cohort_fastas) != 464
            or len(set(cohort_fastas)) != 464 or len(assembly_map) != 466
            or set(assembly_map.values()) != set(all_fastas)
            or set(cohort_fastas) != {v for k, v in assembly_map.items() if k not in REFERENCES}):
        raise ValueError('466/464 cohort/reference identity mismatch')
    recovery = base/'analysis/PAM_Index_Recovery_2026-10-29_v1/verified_assemblies.tsv'
    recovered = {r['assembly_id']: r for r in csv.DictReader(recovery.open(), delimiter='\t')}
    sources[str(recovery)] = dict(sha256=digest(recovery), bytes=recovery.stat().st_size)
    inputs, sbds_rows = [], []
    historical = {}
    aggregate = collections.Counter()
    perfect_distribution = collections.Counter()
    for assembly, fasta in assembly_map.items():
        old_id = historical_id(fasta)
        old = base/'results/SBDS_all_HPRCv2'/old_id
        try:
            summary = pairs(old/'SBDS.summary.tsv')
            summary_error = ''
        except (OSError, ValueError) as exc:
            summary, summary_error = {}, str(exc)
        try:
            result = historical_sbds(old/'SBDS.paf')
            if 'paf_records' in summary and result['rows'] != int(summary['paf_records']):
                raise ValueError('SBDS summary/PAF mismatch')
            result.update(status='READ_OK', error='')
            aggregate.update({k: v for k, v in result.items() if k in ['rows', 'duplicate_reported_anchors'] or k.startswith('cfd_ge_')})
            perfect_distribution[result['perfect_20eq']] += 1
        except (OSError, ValueError) as exc:
            result = dict(status='FAILED_HISTORICAL_READ', error=str(old/'SBDS.paf')+': '+str(exc),
                          rows=None, perfect_20eq=None, sha256=None, bytes=None)
        result['summary_error'] = summary_error
        historical[assembly] = result
        sbds_rows.append(dict(assembly=assembly, **{k: result[k] for k in ['status', 'error', 'rows', 'perfect_20eq', 'sha256', 'bytes']}))
        resolved = recovered.get(assembly)
        original_prefix = summary.get('index_prefix_used', str(base/'indices'/old_id/old_id))
        prefix = resolved['index_prefix'] if resolved else original_prefix
        missing = [s for s in SUFFIXES if not Path(prefix+'.'+s).is_file()]
        fai = Path(fasta+'.fai')
        lengths = [int(s.split('\t')[1]) for s in fai.read_text().splitlines()]
        if not lengths or sum(lengths) >= 2**32:
            raise ValueError('unsupported assembly dimensions: '+assembly)
        inputs.append(dict(assembly_id=assembly, historical_result_id=old_id, cohort='reference' if assembly in REFERENCES else 'HPRCv2',
                           fasta=fasta, fasta_bytes=Path(fasta).stat().st_size,
                           fai_sha256=digest(fai), bases=sum(lengths), contigs=len(lengths),
                           original_index=original_prefix, resolved_index=prefix,
                           missing_component_metadata=','.join(missing),
                           content_evidence='PRIOR_12_SET_VERIFICATION_NOT_REHASHED' if resolved else 'NOT_CONTENT_VERIFIED',
                           source_sha256=resolved['fasta_sha256'] if resolved else 'PENDING_FULL_READ'))
    tsv(out/'assemblies.tsv', list(inputs[0]), inputs)
    tsv(out/'historical_sbds.tsv', list(sbds_rows[0]), sbds_rows)
    (out/'historical_sbds.json').write_text(json.dumps(dict(assemblies=historical, aggregate=dict(aggregate),
                perfect_match_distribution=dict(perfect_distribution),
                caveat='Original reported records only; not strict-PAM/PAF eligibility, homology or biological duplication validation.'), indent=2)+'\n')
    rows = library(base/'guides/brunello_77441.fa')
    clean = list(csv.DictReader((base/'guides/brunello_clean.tsv').open(), delimiter='\t'))
    original = collections.Counter((r['gene_symbol'].replace(' ', '_')+'_'+r['guide'], r['guide']) for r in clean)
    raw_table = list(csv.DictReader((base/'guides/Brunello_unix.txt').open(), delimiter='\t'))
    raw_metadata = collections.Counter((r['Target Gene ID'], r['Target Gene Symbol'], r['sgRNA Target Sequence'], r['PAM Sequence']) for r in raw_table)
    clean_metadata = collections.Counter((r['gene_id'], r['gene_symbol'], r['guide'], r['pam']) for r in clean)
    if raw_metadata != clean_metadata or (base/'guides/Brunello.txt').read_text() != (base/'guides/Brunello_unix.txt').read_text():
        raise ValueError('raw/Unix/clean Brunello metadata mismatch')
    if len(rows) != 77441 or sources[str(base/'guides/brunello_77441.fa')]['sha256'] != LIBRARY_SHA256 or collections.Counter(rows) != original:
        raise ValueError('Brunello FASTA differs from frozen run or clean source table')
    seq_counts = collections.Counter(s for _, s in rows)
    pilot = choose_pilot(rows)
    tsv(out/'brunello_pilot.tsv', list(pilot[0]), pilot)
    (out/'brunello_pilot.fa').write_text(''.join('>'+r['guide_id']+'\n'+r['sequence']+'\n' for r in pilot))
    (out/'sbds.fa').write_text('>SBDS\n'+SBDS+'\n')
    run_rows, totals = [], {}
    for name in ['Brunello_all_HPRCv2', 'Brunello_all_HPRCv2_scored_unscored_v2']:
        root = base/'results'/name
        total, statuses = collections.Counter(), collections.Counter()
        progress = (out/(name+'-progress.jsonl')).open('x', buffering=1)
        for assembly in assembly_map:
            directory = root/historical_id(assembly_map[assembly])
            status = 'NO_DONE_MARKER'
            if (directory/'brunello.done').is_file():
                status = 'SMALL_ARTIFACTS_VERIFIED_LARGE_PAF_NOT_REREAD'
                try:
                    s = pairs(directory/'brunello.summary.tsv')
                    meta = pairs(directory/'brunello.metadata.tsv')
                    checksums = {line.split()[1]: line.split()[0] for line in (directory/'brunello.checksums.sha256').read_text().splitlines()}
                    for file in ['brunello.summary.tsv', 'brunello.metadata.tsv', 'brunello.time']:
                        if (directory/file).stat().st_size > 2*1024*1024:
                            raise ValueError('deferred large historical artifact: '+str(directory/file))
                        if digest(directory/file) != checksums[file]:
                            raise ValueError('historical small-artifact checksum mismatch: '+str(directory/file))
                    if meta['guides_sha256'] != LIBRARY_SHA256:
                        raise ValueError('historical library mismatch')
                    if 'scored_candidates' in s:
                        if int(s['total_candidates']) != int(s['scored_candidates'])+int(s['unscored_candidates']):
                            raise ValueError('historical arithmetic mismatch')
                    for file in ['brunello_CFD020.paf.zst', 'brunello.unscored.paf.zst']:
                        if file in checksums and not (directory/file).is_file():
                            raise ValueError('missing historical saved stream')
                    total.update({k: int(v) for k, v in s.items() if k != 'metric'})
                except (OSError, ValueError, KeyError) as exc:
                    status = 'FAILED_SMALL_ARTIFACT_VALIDATION: '+str(exc)
            statuses[status] += 1
            run_rows.append(dict(run=name, assembly=assembly, status=status))
            progress.write(json.dumps(run_rows[-1])+'\n')
        progress.close()
        totals[name] = dict(status_counts=dict(statuses), summary_totals=dict(total))
    tsv(out/'historical_brunello_status.tsv', list(run_rows[0]), run_rows)
    report = dict(status='INITIAL_AUDIT_COMPLETE_NOT_A_RUN_AUTHORIZATION',
                  cohort=dict(hprcv2=464, references=sorted(REFERENCES), sbds_total=466, brunello_total=464),
                  source_files=sources, historical_sbds=dict(status_counts=dict(collections.Counter(r['status'] for r in sbds_rows)),
                      failures=[r for r in sbds_rows if r['status'] != 'READ_OK'],
                      aggregate_readable_only=aggregate, perfect_match_distribution_readable_only=dict(perfect_distribution)), 
                  historical_brunello=totals,
                  library=dict(records=len(rows), unique_ids=len(rows), unique_sequences=len(seq_counts),
                               duplicate_sequence_groups=sum(n > 1 for n in seq_counts.values()),
                               gene_symbols=len({r['gene_symbol'] for r in clean}),
                               non_targeting_controls=sum(r['gene_symbol'] == 'Non-Targeting Control' for r in clean),
                               raw_unix_clean_gene_sequence_pam_multiset_match=True,
                               raw_unix_universal_newline_text_identical=True,
                               clean_table_id_sequence_match_with_documented_label_mapping=True,
                               existing_fasta_label_mapping='Spaces in source gene symbol become underscores; frozen FASTA IDs are preserved unchanged.',
                               published_library_version_attestation='NOT_INDEPENDENTLY_DOWNLOADED; local raw/clean/FASTA digests retained'),
                  index_content_status=dict(prior_verified_subset=len(recovered), pending=len(inputs)-len(recovered),
                                            warning='Stat/is_file is not a full read; damaged original sets must not be reused.'),
                  storage={str(p): dict(zip(['total', 'used', 'free'], shutil.disk_usage(p))) for p in [base, Path('/export/local/home/salehi')]},
                  pilot=dict(assembly='HG00735_pat', guides=pilot, sampling='GC quantiles; three-guide safety/feasibility sample, no library-wide statistical precision'))
    (out/'audit.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    try:
        result = audit(a.base, a.out)
    except Exception as exc:
        if a.out.is_dir():
            (a.out/'INCOMPLETE').write_text(str(exc)+'\n')
        raise
    print(json.dumps({k: result[k] for k in ['status', 'cohort', 'library', 'index_content_status', 'historical_sbds', 'historical_brunello', 'pilot']}, indent=2))


if __name__ == '__main__':
    main()
