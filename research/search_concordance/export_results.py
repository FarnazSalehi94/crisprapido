#!/usr/bin/env python3
"""Publish small complete case/record ledgers; keep raw outputs externally."""
import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'crispr_bulge'))
from adapter import digest

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--run', type=Path, required=True)
p.add_argument('--out', type=Path, required=True, help='NEW output directory')
a = p.parse_args(); a.out.mkdir(parents=True, exist_ok=False)
summary = json.loads((a.run/'summary.json').read_text())
summary['external_summary_sha256'] = digest(a.run/'summary.json')
summary['audited_rows_sha256'] = {}
cases = []
columns = ['panel','guide_sequence','mode','row','contig','strand','anchor','start','end','cigar','schema',
    'eligible','comparable_record','export_errors','alignment_errors','mismatches',
    'inserted','deleted','events','max_gap','pam','independent_cfd_status','independent_cfd','reported_cfd','cfd_agrees']
with (a.out/'records.tsv').open('x', newline='') as output:
    writer = csv.writer(output, delimiter='\t', lineterminator='\n'); writer.writerow(columns)
    for panel in summary['panels']:
        command = panel['commands']['window_wfa2']
        panel['guide_sequence'] = command[command.index('-g')+1]
        path = a.run/panel['name']/'audited-rows.json'
        summary['audited_rows_sha256'][panel['name']] = digest(path)
        parsed = json.loads(path.read_text())
        for mode, rows in parsed.items():
            for row in rows:
                status = 'not_auditable' if not row['eligible'] else 'unavailable' if row['independent_cfd'] is None else 'defined'
                data = dict(panel=panel['name'], guide_sequence=panel['guide_sequence'], mode=mode, **row, contig=row['key'][0], strand=row['key'][1], anchor=row['key'][2],
                            independent_cfd_status=status)
                writer.writerow([';'.join(map(str, data[k])) if isinstance(data[k], list) else
                    'N/A' if data[k] is None else data[k] for k in columns])
        for case in panel.pop('ledger'):
            case['panel'] = panel['name']
            case['record_ids'] = {mode:[r['row'] for r in rows] for mode, rows in case.pop('rows').items()}
            if 'window_missing_investigation' in case:
                # Keep the controlled no-filter preferred rows: they are separate
                # diagnostics, not production predictions or performance input.
                for row in case['window_missing_investigation']['preferred_unfiltered_same_contig_strand']:
                    row.pop('aligned_guide'); row.pop('aligned_target')
            cases.append(case)
        for pair in panel['pairs'].values():
            examples = pair.pop('differing_pairs')
            pair['examples'] = [{k:v for k,v in example.items() if k not in ('left','right')} |
                {side:{k:v for k,v in example[side].items() if k not in ('aligned_guide','aligned_target')} for side in ('left','right')}
                for example in examples]
(a.out/'cases.json').write_text(json.dumps(cases, indent=2, sort_keys=True)+'\n')
summary['case_ledger_sha256'] = digest(a.out/'cases.json')
summary['record_ledger_sha256'] = digest(a.out/'records.tsv')
(a.out/'summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True)+'\n')
print('Published', len(cases), 'fully linked anchor cases and all raw-record audit outcomes')
