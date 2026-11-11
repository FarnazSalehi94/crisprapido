#!/usr/bin/env python3
"""Data-derived safety-pilot diagnostics, NOT the six final cohort figures."""
import argparse
import csv
import json
from pathlib import Path

from collect import only
from inventory import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--historical-comparison', type=Path)
    a = p.parse_args()
    report = json.loads(a.report.read_text())
    if report['status'] != 'COMPLETE_SAFETY_PILOT':
        raise ValueError('cannot plot a missing/incomplete result as zero')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    a.out.mkdir()
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    configs = report['configurations']
    names = list(configs)
    colors = ['#1b7c83', '#be6841']
    axes = axes.ravel()
    if report['project'] == 'brunello':
        guides = list(configs['A']['qc']['per_guide'])
        x = np.arange(len(guides))
        for i, name in enumerate(names):
            y = [configs[name]['qc']['per_guide'][g]['eligible_hits'] for g in guides]
            axes[0].bar(x+(i-.5)*.36, y, width=.36, label=name, color=colors[i])
        axes[0].set_xticks(x, [g.rsplit('_', 1)[0] for g in guides])
        axes[0].set_yscale('log')
        axes[0].set_ylabel('Eligible anchors (log scale)')
        axes[0].set_title('A: m4/b1/z1; B: m5/b1/z2')
        axes[0].legend()
        for i, name in enumerate(names):
            stats = configs[name]['qc']
            total, supported = stats['eligible_hits'], stats['cb_statuses'].get('scored', 0)
            fraction = supported/total
            axes[1].bar(i, fraction, color=colors[i])
            axes[1].bar(i, 1-fraction, bottom=fraction, color='#d7d7d7')
            axes[1].text(i, fraction/2, str(supported)+'\n'+format(100*fraction, '.1f')+'%', ha='center', va='center', color='white')
        axes[1].set_xticks(range(2), names)
        axes[1].set_ylim(0, 1.1)
        axes[1].set_ylabel('Fraction of eligible hits')
        axes[1].set_title('Neural defined (colour) / unavailable (grey)')
        labels = ['none', '0.05', '0.20']
        for i, name in enumerate(names):
            stats = configs[name]['qc']['per_guide']
            axes[2].plot(labels, [sum(g['view_'+v+'_rows'] for g in stats.values()) for v in labels],
                         'o-', color=colors[i], label=name)
        axes[2].set_yscale('log')
        axes[2].set_xlabel('Analytical full-precision CFD threshold')
        axes[2].set_ylabel('Rows in downstream view (log)')
        axes[2].set_title('No score filtering of eligible hits; all retained')
        axes[2].legend()
        timing_colors = ['#375b8c', '#c6d2e4', '#a58cb8', '#d8c4a8']
        bottom = np.zeros(2)
        series = [([configs[n]['native_search_seconds_including_replay'] for n in names], 'Native traversal/replay'),
                  ([configs[n]['other_search_wall_residual_seconds'] for n in names], 'Other search work (not isolated setup)'),
                  ([configs[n]['scoring_metrics']['wall_seconds'] for n in names], 'Post-search both-score annotation'),
                  ([configs[n]['monitor']['elapsed_seconds']-configs[n]['search_metrics']['wall_seconds']-configs[n]['scoring_metrics']['wall_seconds'] for n in names], 'Other checks/startup')]
        for (values, label), color in zip(series, timing_colors):
            axes[3].bar(range(2), values, bottom=bottom, color=color, label=label)
            bottom += values
        axes[3].set_xticks(range(2), names)
        axes[3].set_ylabel('Wall seconds')
        axes[3].set_title('One fresh call/config; uncontrolled caches')
        axes[3].legend(fontsize=7)
        title = 'Brunello SAFETY PILOT: three fixed guides, one HG00735_pat assembly'
    else:
        if not a.historical_comparison:
            raise ValueError('SBDS diagnostics need historical comparison')
        historical = json.loads(a.historical_comparison.read_text())
        counts = historical['counts']
        values = [counts['old_raw_rows'], counts['old_strict_eligible_witnesses'],
                  configs['historical-domain-control']['qc']['eligible_hits'], configs['primary']['qc']['eligible_hits']]
        bars = axes[0].bar(range(4), values, color=['#c0c0c0', '#a3b8bb', '#1b7c83', '#be6841'])
        bars[0].set_hatch('///')
        axes[0].set_xticks(range(4), ['Old raw\n(not comparable)', 'Old eligible\nwitnesses', 'Corrected\nm2/b1/z1', 'Primary\nm5/b1/z2'])
        axes[0].set_yscale('log')
        axes[0].set_ylabel('Reported rows / anchors (log)')
        axes[0].set_title('Eligibility and setting effects; not a recall ratio')
        for i, name in enumerate(['historical-domain-control', 'primary']):
            stats = configs[name]['qc']
            defined = stats['cb_statuses'].get('scored', 0)
            axes[1].bar(i, defined, color='#1b7c83')
            axes[1].bar(i, stats['eligible_hits']-defined, bottom=defined, color='#d7d7d7')
            axes[1].text(i, stats['eligible_hits']+500, str(defined)+' / '+str(stats['eligible_hits']), ha='center')
        axes[1].set_xticks(range(2), ['m2 control', 'm5 primary'])
        axes[1].set_ylim(0, 1.16*max(v['qc']['eligible_hits'] for v in configs.values()))
        axes[1].set_ylabel('Hits (defined / all labelled above)')
        axes[1].set_title('Neural defined (colour); unavailable retained (grey)')
        marker_path = Path(configs['primary']['completion_marker'])
        if digest(marker_path) != configs['primary']['completion_sha256']:
            raise ValueError('completion marker changed')
        marker = json.loads(marker_path.read_text())
        annotation = only(marker['outputs'], 'annotations.tsv')
        if digest(annotation) != marker['outputs'][str(annotation)]['sha256']:
            raise ValueError('annotation changed')
        with annotation.open() as f:
            paired = [(float(r['cfd']), float(r['crispr_bulge'])) for r in csv.DictReader(f, delimiter='\t')
                      if r['cfd'] != 'N/A' and r['crispr_bulge'] != 'N/A']
        axes[2].scatter([p[0] for p in paired], [p[1] for p in paired], s=3, alpha=.3, color='#375b8c', rasterized=True)
        axes[2].set_xscale('symlog', linthresh=.001)
        axes[2].set_xlabel('Stored CFD (four decimals; symlog retains true zero)')
        axes[2].set_ylabel('Experimental CRISPR-Bulge score')
        axes[2].set_title('Paired supported subset only (n='+str(len(paired))+')')
        axes[3].axis('off')
        records_path = a.historical_comparison.parent/'records.tsv'
        if digest(records_path) != historical['records_sha256']:
            raise ValueError('historical comparison records changed')
        with records_path.open() as f:
            records = list(csv.DictReader(f, delimiter='\t'))
        examples = []
        for label, cigar in [('chr7 pattern', '20='), ('chr12 pattern', '6=1I4=1X6=1X1='), ('chr17 pattern', '2=1I1X6=1X9=')]:
            found = [r for r in records if r['old_cigar'] == cigar]
            if not found:
                raise ValueError('missing historical pattern; do not invent a table row')
            for r in found:
                status = 'ineligible PAM' if r['pam'][1:] != 'GG' else 'retained' if r['primary_same_anchor'] == 'True' else 'investigate'
                examples.append([label, r['pam'], r['old_cfd'], status])
        table = axes[3].table(cellText=examples,
                    colLabels=['Prior context', 'Observed PAM', 'Old CFD', 'Corrected status'], loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.6)
        axes[3].set_title('Historical patterns: chromosome context from prior study')
        title = 'SBDS SAFETY PILOT: one HG00735_pat assembly, NOT a 466-assembly result'
    fig.suptitle(title, fontsize=12)
    fig.text(.5, .01, 'Strict NGG, f=0.75; no neural/CFD pruning. Scores are not calibrated cleavage probabilities or interchangeable. No genome-wide recall claim.', ha='center', fontsize=8)
    fig.tight_layout(rect=[0, .04, 1, .95])
    outputs = []
    for extension in ['pdf', 'png']:
        path = a.out/('pilot-diagnostics.'+extension)
        fig.savefig(path, dpi=180)
        outputs.append(path)
    plt.close(fig)
    with (a.out/'provenance.json').open('x') as f:
        json.dump(dict(scope='Safety-pilot diagnostic, not final SBDS six figures or population analysis',
                       source_report=str(a.report), report_sha256=digest(a.report), matplotlib=matplotlib.__version__, numpy=np.__version__,
                       renderer_sha256=digest(__file__), outputs={str(p): digest(p) for p in outputs}), f, indent=2)
        f.write('\n')
    print('\n'.join(str(p) for p in outputs))


if __name__ == '__main__':
    main()
