#!/usr/bin/env python3
"""Bounded search-mode concordance. No arbitrary genome input, no algorithm edits.

Three current complete pipelines are comparable. Two read-only historical
executables are explicitly schema/policy-limited diagnostics, not a five-way
fair performance ranking. Raw output is never repaired or overwritten.
"""
import argparse
from collections import Counter
import csv
import itertools
import json
from pathlib import Path
import re
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'research/crispr_bulge'), str(ROOT/'research/pam_anchored_poc/cfd_feasibility')]
from adapter import digest
from study import fixtures, load_tables, oracle, rc, representatives, independent_cfd, templates, alignment

CURRENT = ('window_wfa2', 'columba_wfa2', 'pam_columba')
ALL = (*CURRENT, 'historical_wfa2', 'sassy_diagnostic')


def run(cmd, **kw):
    try:
        return subprocess.run(list(map(str, cmd)), check=True, cwd=ROOT, timeout=120, **kw)
    except subprocess.CalledProcessError as error:
        if error.stderr:
            print(error.stderr, file=sys.stderr)
        raise


def rows(text, guide, refs, m, z, tables):
    """Replay raw exported geometry. Schema adaptation is labelled, never repair.

    Sassy's eleven-base-field export omits block length. Reading its known field
    positions is diagnostic only; no missing field or corrected span is invented.
    """
    result = []
    for n, raw in enumerate(text.splitlines()):
        f = raw.split('\t')
        if not raw or not f[0] == 'Guide':
            raise ValueError('unexpected non-PAF stdout')
        schema = 'paf12' if len(f) > 11 and f[11] == '255' else 'legacy11_missing_block_length'
        start_tags = 12 if schema == 'paf12' else 11
        tags = {x.split(':', 2)[0]: x.split(':', 2)[2] for x in f[start_tags:]}
        lo, hi = int(f[7]), int(f[8]); seq = refs[f[5]]
        key = (f[5], f[4], hi if f[4] == '+' else lo)
        errors = []; export = []
        if schema != 'paf12':
            export.append('missing_mandatory_block_length')
        if list(map(int, f[1:4])) != [len(guide), 0, len(guide)]:
            export.append('query_coordinates_not_full_guide')
        if int(f[6]) != len(seq) or not 0 <= lo <= hi <= len(seq):
            errors.append('reference_coordinates')
        cigar = tags.get('cg', '')
        ops = [(int(n), op) for n, op in re.findall(r'(\d+)([=XMID])', cigar)]
        if not ops or ''.join(str(n)+op for n, op in ops) != cigar:
            errors.append('malformed_cigar')
        consumed_g = sum(n for n, op in ops if op != 'D')
        consumed_t = sum(n for n, op in ops if op != 'I')
        if consumed_g != len(guide):
            errors.append('cigar_guide_consumption')
        if consumed_t != hi-lo:
            errors.append('cigar_reference_span_disagrees_with_export')
        target = seq[lo:hi] if f[4] == '+' else rc(seq[lo:hi])
        pam = seq[hi:hi+3] if f[4] == '+' else rc(seq[max(0, lo-3):lo])
        mm = matches = ins = dele = events = gaplen = maxgap = 0
        ag = []; at = []; expanded = []; qi = ti = 0; previous = False
        if not errors:
            for count, op in ops:
                gap = op in 'ID'
                if gap:
                    events += not previous; gaplen = gaplen+count if previous else count
                    maxgap = max(maxgap, gaplen)
                else:
                    gaplen = 0
                previous = gap
                for _ in range(count):
                    a = '-' if op == 'D' else guide[qi]
                    b = '-' if op == 'I' else target[ti]
                    qi += op != 'D'; ti += op != 'I'
                    ins += op == 'I'; dele += op == 'D'
                    if not gap:
                        if (op == '=' and a != b) or (op == 'X' and a == b):
                            errors.append('cigar_base_assertion')
                        mm += a != b; matches += a == b
                        op2 = '=' if a == b else 'X'
                    else:
                        op2 = op
                    ag.append(a); at.append(b); expanded.append(op2)
        if len(pam) != 3 or set(pam)-set('ACGT') or pam[1:] != 'GG':
            errors.append('not_requested_complete_NGG')
        if not errors and (mm > m or events > 1 or maxgap > z or matches < 15):
            errors.append('outside_common_mismatch_bulge_fraction_limits')
        score = None
        if not errors:
            canonical = [(len(list(xs)), op) for op, xs in itertools.groupby(expanded)]
            score = independent_cfd(''.join(ag), ''.join(at), canonical, pam, *tables)
        reported = float(tags['cf']) if 'cf' in tags else None
        score_agrees = (abs(reported-score) <= .000051 if reported is not None and score is not None else
                        None if reported is None else False)
        if 'nm' in tags and not errors and int(tags['nm']) != mm:
            export.append('reported_mismatch_count_disagrees')
        if 'ng' in tags and not errors and int(tags['ng']) != events:
            export.append('reported_gap_events_disagree')
        result.append(dict(row=n, key=list(key), start=lo, end=hi, cigar=cigar, schema=schema,
            export_errors=sorted(set(export)), alignment_errors=sorted(set(errors)),
            eligible=not errors, comparable_record=not errors and not export,
            mismatches=mm, inserted=ins, deleted=dele, events=events, max_gap=maxgap,
            pam=pam, independent_cfd=score, reported_cfd=reported, cfd_agrees=score_agrees,
            aligned_guide=''.join(ag), aligned_target=''.join(at)))
    return result


def site_map(records, strict=True):
    out = {}
    for row in records:
        if row['comparable_record' if strict else 'eligible']:
            out.setdefault(tuple(row['key']), []).append(row)
    return out


def witness_at(refs, guide, key, m, z, tables):
    """Anchor-local exhaustive witness search, separate from any emitted CIGAR."""
    contig, strand, anchor = key; seq = refs[contig]; oriented = seq if strand == '+' else rc(seq)
    end = anchor if strand == '+' else len(seq)-anchor
    found = []
    if not 0 <= end <= len(oriented)-3:
        return found
    pam = oriented[end:end+3]
    if len(pam) != 3 or set(pam)-set('ACGT') or pam[1:] != 'GG':
        return found
    for ops in templates(len(guide), z):
        span = len(ops)-ops.count('I')
        if span <= 0 or span > end:
            continue
        lo, hi = (end-span, end) if strand == '+' else (len(seq)-end, len(seq)-end+span)
        hit = alignment(guide, oriented[end-span:end], ops, key, lo, hi, pam, tables, m, .75)
        if hit:
            found.append(hit)
    return found


def witness_sam(key, hit):
    ops = re.findall(r'(\d+)([=XID])', hit.cigar)
    cigar = ''.join(n+op for n, op in (reversed(ops) if key[1]=='-' else ops))
    distance = hit.mm+hit.ins+hit.dele
    return f'guide\t{16 if key[1]=="-" else 0}\t{key[0]}\t{hit.start+1}\t255\t{cigar}\t*\t0\t0\t*\t*\tNM:i:{distance}\tAS:i:{distance}\n'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('cli', 'native', 'columba', 'builder', 'historical', 'sassy', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args(); a.out.mkdir(exist_ok=False, parents=True)
    guide, allrefs, _ = fixtures(); tables = load_tables()
    names = ['suffix_TCAC', 'suffix_ACAC', 'suffix_AAAC', 'suffix_AAAA', 'suffix_GGGG',
             'D5_1', 'D5_2', 'I5_1', 'I5_2', 'D0_1', 'I0_1', 'D20_1', 'mixed5']
    refs = {name+strand: allrefs[name+strand] for name in names for strand in '+-'}
    refs['dense_plus'] = allrefs['suffix_TCAC+']+'TTTT'+allrefs['suffix_TCAC+']
    refs['dense_minus'] = rc(refs['dense_plus'])
    mixed = guide[:5]+'C'+guide[6:]+'AGG'
    datasets = [('common_z1', refs, 3, 1), ('z2_radius_diagnostic', refs, 3, 2),
                ('mixed_only', {'mixed+': mixed, 'mixed-': rc(mixed)}, 0, 2)]
    sources = [a.cli, a.native, a.columba, a.builder, a.historical, a.sassy, Path(__file__),
               ROOT/'src/main.rs', ROOT/'src/verification.rs', ROOT/'src/reporting.rs',
               ROOT/'src/columba.rs', ROOT/'research/pam_anchored_poc/poc.cpp',
               ROOT/'mismatch_scores.txt', ROOT/'pam_scores.txt']
    pins = {str(x): digest(x) for x in sources}
    summary = dict(baseline='ce52c15', current_applicable_modes=list(CURRENT), historical_modes='diagnostics only; externally restricted PAMs, nonstandard/invalid exports, non-equivalent z2 Sassy budget, debug Sassy build', panels=[], sha256=pins)
    for label, sequences, m, z in datasets:
        folder = a.out/label; folder.mkdir(); ref = folder/'reference.fa'; guides = folder/'guide.fa'; idx = folder/'index'
        ref.write_text(''.join(f'>{c}\n{s}\n' for c, s in sequences.items())); guides.write_text('>guide\n'+guide+'\n')
        with (folder/'index-build.log').open('x') as log:
            run([sys.executable, ROOT/'research/pam_anchored_poc/compatibility/index_manifest.py',
                 'build', '--builder', a.builder, '--reference', ref, '--index', idx], stdout=log, stderr=subprocess.STDOUT)
        truth = representatives(oracle(sequences, guide, 'GG', m, z, .75, tables))
        shared = ['-r', ref, '-g', guide, '-m', m, '-b', '1', '-z', z, '-t', '1']
        current = [a.cli, *shared, '-p', 'GG', '-f', '.75']
        commands = dict(window_wfa2=current,
            columba_wfa2=[*current, '--columba-bin', a.columba, '--columba-index', idx],
            pam_columba=[*current, '--search-backend', 'pam-anchored', '--pam-search-bin', a.native, '--columba-index', idx],
            historical_wfa2=[a.historical, *shared],
            sassy_diagnostic=[a.sassy, *shared, '-p', 'GG', '-f', '.75', '--allow-indels',
                              '--mismatch-scores', ROOT/'mismatch_scores.txt', '--pam-scores', ROOT/'pam_scores.txt'])
        parsed = {}; measures = []; duplicates = {}; maps = {}; diagnostics = {}
        for repetition in range(3 if label == 'common_z1' else 1):
            order = list(ALL)
            order = order[repetition:]+order[:repetition]
            for mode in order:
                prefix = folder/f'{mode}-{repetition}'
                with Path(str(prefix)+'.monitor.stdout').open('x') as log:
                    run([sys.executable, ROOT/'research/pam_anchored_poc/monitored_pilot/supervise.py',
                         '--limit-kib', '2097152', '--interval', '.02', '--log', str(prefix)+'.rss',
                         '--summary', str(prefix)+'.tree.json', '--', sys.executable, ROOT/'research/crispr_bulge/measure.py',
                         '--metrics', str(prefix)+'.time.json', '--stdout', str(prefix)+'.paf', '--stderr', str(prefix)+'.stderr',
                         '--timeout', '120', '--', *commands[mode]], stdout=log)
                output = Path(str(prefix)+'.paf'); records = rows(output.read_text(), guide, sequences, m, z, tables)
                for record in records:
                    if mode in CURRENT and record['comparable_record']:
                        expected_score, reported_score = record['independent_cfd'], record['reported_cfd']
                        if (expected_score is None) != (reported_score is None) or (expected_score is not None and not record['cfd_agrees']):
                            raise ValueError('current pipeline CFD disagrees with independent replay')
                if repetition == 0:
                    parsed[mode] = records
                    maps[mode] = site_map(records)
                    diagnostics[mode] = site_map(records, strict=False)
                    duplicates[mode] = sum(len(v)-1 for v in maps[mode].values())
                elif digest(output) != digest(folder/f'{mode}-0.paf'):
                    # Historical threaded exporters may reorder rows: do not assert byte identity.
                    signature = lambda rs: sorted(json.dumps(r, sort_keys=True) for r in rs)
                    clean = lambda rs: [{k:v for k,v in r.items() if k != 'row'} for r in rs]
                    if signature(clean(records)) != signature(clean(parsed[mode])):
                        raise ValueError('nondeterministic records')
                measures.append(dict(mode=mode, repetition=repetition,
                    measurement=json.loads(Path(str(prefix)+'.time.json').read_text()),
                    tree=json.loads(Path(str(prefix)+'.tree.json').read_text()), output_sha256=digest(output)))
        # Diagnostic only: expose the window-preferred alignment even when it fails
        # eligibility/PAM. These rows NEVER enter the primary comparison or timing.
        unfiltered = run([*current, '--no-filter'], capture_output=True, text=True).stdout
        (folder/'window-unfiltered-diagnostic.paf').write_text(unfiltered)
        preferred_window = rows(unfiltered, guide, sequences, m, z, tables)
        # Raw candidate controls: same upper-bound and exact searches as automatic CLI.
        sams = []
        for distance in sorted({0, m+z}):
            sam = folder/f'candidates-e{distance}.sam'
            with (folder/f'candidates-e{distance}.log').open('x') as log:
                run([a.columba, '-r', idx, '-f', guides, '-a', 'all', '-m', 'edit', '-e', distance, '-t', '1', '-o', sam], stdout=log, stderr=subprocess.STDOUT)
            sams += [line for line in sam.read_text().splitlines() if not line.startswith('@')]
        merged = folder/'candidates-union.sam'; merged.write_text('\n'.join(dict.fromkeys(sams))+'\n')
        imported = run([*current, '--columba-sam', merged], capture_output=True, text=True).stdout
        if imported != (folder/'columba_wfa2-0.paf').read_text():
            raise ValueError('raw-candidate control differs from automatic pipeline')
        missing_columba = set(truth)-set(maps['columba_wfa2'])
        injected = folder/'oracle-witnesses.sam'
        with injected.open('x') as target:
            for k in sorted(missing_columba):
                target.write(witness_sam(k, truth[k]))
        injected_out = run([*current, '--columba-sam', injected], capture_output=True, text=True).stdout if missing_columba else ''
        injected_records = rows(injected_out, guide, sequences, m, z, tables) if injected_out else []
        anchored_unfiltered = run([*current, '--columba-sam', injected, '--no-filter'], capture_output=True, text=True).stdout if missing_columba else ''
        (folder/'anchored-unfiltered-diagnostic.paf').write_text(anchored_unfiltered)
        preferred_anchored = rows(anchored_unfiltered, guide, sequences, m, z, tables) if anchored_unfiltered else []
        rescued = set(site_map(injected_records))
        sam_intervals = []
        for line in sams:
            f = line.split('\t'); flag = int(f[1])
            if flag & 4 or f[2] == '*': continue
            lo = int(f[3])-1; span = sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', f[5]) if op in 'MDN=X')
            sam_intervals.append((f[2], '-' if flag & 16 else '+', lo, lo+span))
        universe = set(truth).union(*(set(v) for v in diagnostics.values()),
                                   *(set(tuple(r['key']) for r in parsed[mode]) for mode in ALL))
        ledger = []; pair_counts = {}
        for key in sorted(universe):
            present = [mode for mode in ALL if key in (maps if mode in CURRENT else diagnostics)[mode]]
            witness = witness_at(sequences, guide, key, m, z, tables)
            if not witness and any(key in maps[mode] for mode in CURRENT):
                raise ValueError('modern eligible hit has no independent witness')
            entry = dict(key=list(key), independent_eligible=bool(witness), witness_alignments=len(witness), present=present,
                         missing=[mode for mode in ALL if mode not in present], rows={mode:[r for r in parsed[mode] if tuple(r['key'])==key] for mode in ALL})
            if key in truth:
                hit = truth[key]; entry['oracle_representative'] = dict(start=hit.start, end=hit.end, cigar=hit.cigar, mm=hit.mm, ins=hit.ins, dele=hit.dele, cfd=hit.score)
            if key in truth and key not in maps['window_wfa2']:
                entry['window_missing_investigation'] = dict(
                    preferred_unfiltered_same_contig_strand=[r for r in preferred_window if r['key'][:2]==list(key[:2])],
                    interpretation='single ends-free optimum per window/strand then constraints/PAM; eligible alternative witness exists; see controlled preferred-alignment rows and source')
            if key in missing_columba:
                hit = truth[key]
                entry['columba_missing_investigation'] = dict(
                    identical_candidate_interval=any((c, s, lo, hi)==(key[0], key[1], hit.start, hit.end) for c,s,lo,hi in sam_intervals),
                    candidate_at_same_anchor=any(c==key[0] and s==key[1] and (hi if s=='+' else lo)==key[2] for c,s,lo,hi in sam_intervals),
                    any_candidate_same_contig_strand=any(c==key[0] and s==key[1] for c,s,_,_ in sam_intervals),
                    oracle_candidate_injection_recovers_anchor=key in rescued,
                    preferred_anchored_unfiltered_same_contig_strand=[{k:v for k,v in r.items() if k not in ('aligned_guide','aligned_target')} for r in preferred_anchored if r['key'][:2]==list(key[:2])],
                    interpretation='eligible witness supplied: recovery localizes loss upstream of/sensitive to candidate representation' if key in rescued else 'eligible witness supplied but pipeline still omits anchor: verifier/fallback/report selection limitation')
            if (len(present) != len(ALL) or len({(r['start'],r['end'],r['cigar']) for mode in ALL for r in entry['rows'][mode] if r['eligible']})>1
                    or any(r['export_errors'] or r['alignment_errors'] or r['cfd_agrees'] is False for mode in ALL for r in entry['rows'][mode])):
                ledger.append(entry)
        for left, right in itertools.combinations(CURRENT, 2):
            common = set(maps[left]) & set(maps[right]); diff = Counter(); examples=[]
            for key in sorted(common):
                x, y = maps[left][key][0], maps[right][key][0]
                flags = {name:x[name]!=y[name] for name in ('start','end','cigar','mismatches','inserted','deleted','events','max_gap')}
                score_diff = (x['independent_cfd'] is None) != (y['independent_cfd'] is None) or (x['independent_cfd'] is not None and y['independent_cfd'] is not None and abs(x['independent_cfd']-y['independent_cfd'])>1e-9)
                for name, flag in flags.items(): diff[name] += flag
                diff['cfd'] += score_diff
                if any(flags.values()): examples.append(dict(key=list(key),left=x,right=y,cfd_differs=score_diff))
            pair_counts[left+'__'+right] = dict(shared=len(common), left_only=len(set(maps[left])-set(maps[right])),
                right_only=len(set(maps[right])-set(maps[left])), differences=dict(diff), differing_pairs=examples)
        mode_summary = {}
        for mode in ALL:
            eligible = maps[mode] if mode in CURRENT else diagnostics[mode]
            if set(eligible)-set(truth):
                raise ValueError('independently eligible exported anchor absent from exhaustive oracle')
            mode_summary[mode] = dict(raw_rows=len(parsed[mode]), standard_valid_sites=len(maps[mode]),
                diagnostic_eligible_sites=len(diagnostics[mode]), missing_expected=len(set(truth)-set(eligible)),
                unexpected_eligible=0, duplicate_standard_anchors=duplicates[mode],
                duplicate_diagnostic_anchors=sum(len(v)-1 for v in diagnostics[mode].values()),
                alignment_rejection_reasons=dict(Counter(e for r in parsed[mode] for e in r['alignment_errors'])),
                export_errors=dict(Counter(e for r in parsed[mode] for e in r['export_errors'])),
                reported_cfd_agreement=dict(Counter(str(r['cfd_agrees']) for r in parsed[mode] if r['eligible'])))
        panel = dict(name=label, bases=sum(map(len,sequences.values())), contigs=len(sequences), m=m,b=1,z=z,f=.75,pam='NGG',
            expected=len(truth), modern_all_shared=len(set.intersection(*(set(maps[k]) for k in CURRENT))),
            five_diagnostic_all_shared=len(set.intersection(*(set(diagnostics[k]) for k in ALL))),
            modes=mode_summary,pairs=pair_counts,ledger=ledger,commands={k:list(map(str,v)) for k,v in commands.items()},
            raw_candidate_control_matches_automatic=True, runs=measures)
        if label=='common_z1':
            panel['performance']={mode:dict(repetitions=3, median_wall_seconds=statistics.median(r['measurement']['wall_seconds'] for r in measures if r['mode']==mode),
                max_individual_rss_kib=max(r['measurement']['max_individual_child_rss_kib'] for r in measures if r['mode']==mode),
                max_tree_rss_kib=max(r['tree']['sampled_peak_rss_kib'] for r in measures if r['mode']==mode),
                ranking_applicable=mode in CURRENT) for mode in ALL}
        (folder/'audited-rows.json').write_text(json.dumps(parsed,indent=2)+'\n')
        summary['panels'].append(panel)
        print(label, panel['expected'], panel['modern_all_shared'], {k:v['diagnostic_eligible_sites'] for k,v in mode_summary.items()}, flush=True)
    for path, sha in pins.items():
        if digest(path)!=sha: raise ValueError('changed executable/source')
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')


if __name__=='__main__': main()
