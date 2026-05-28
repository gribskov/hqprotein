#!/usr/bin/env python3
#
# 05_parse_metrics.py - Load AF3 metrics and run per-metric analysis
# Alex Baena, Mar 2026
#
# Wrapper that loads extracted metrics (from 04_extract_af3_metrics.py)
# and runs selected metric parsers. Each parser gets its own subdirectory
# under --outdir and writes its feature table, plots, etc. there.
#
# Usage:
#   python3 05_parse_metrics.py --indir <metrics_dir> --outdir <results_dir>
#   python3 05_parse_metrics.py --indir <metrics_dir> --outdir <results_dir> --parsers plddt

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from metric_parsers.loader import load_dataset
from metric_parsers import plddt
from metric_parsers import sequence_composition


AVAILABLE_PARSERS = {
    'plddt': plddt,
    'sequence_composition': sequence_composition,
    # 'ptm': ptm,       # future
    # 'pae': pae,        # future
}


def main(metrics_dir, outdir, parsers='all', plots=False, seqdir=None):
    dataset, problems = load_dataset(metrics_dir, seqdir=seqdir)

    if problems:
        for pid, issues in problems:
            print(f'  warning: {pid}: {", ".join(issues)}', file=sys.stderr)

    # select parsers
    if parsers == 'all':
        active = AVAILABLE_PARSERS
    else:
        names = [p.strip() for p in parsers.split(',')]
        unknown = [n for n in names if n not in AVAILABLE_PARSERS]
        if unknown:
            print(f'Unknown parsers: {", ".join(unknown)}. '
                  f'Available: {", ".join(AVAILABLE_PARSERS)}', file=sys.stderr)
            sys.exit(1)
        active = {n: AVAILABLE_PARSERS[n] for n in names}

    # run each parser — each gets its own subdirectory
    outdir = Path(outdir)
    for name, module in active.items():
        parser_dir = outdir / name
        parser_dir.mkdir(parents=True, exist_ok=True)
        module.analyze(dataset, parser_dir, plots=plots)

    # summary
    cats = {}
    for p in dataset:
        cats[p['category']] = cats.get(p['category'], 0) + 1
    print(f'{len(dataset)} proteins, ' +
          ', '.join(f'{n} {c}' for c, n in sorted(cats.items())) +
          f' | parsers: {", ".join(active)} -> {outdir}/')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Load AF3 metrics, run per-metric analysis')
    parser.add_argument('--indir', required=True,
                        help='Directory with extracted metrics (from step 04)')
    parser.add_argument('--outdir', required=True,
                        help='Results directory (subdirs created per parser)')
    parser.add_argument('--parsers', default='all',
                        help='Comma-separated list of parsers to run (default: all)')
    parser.add_argument('--plots', action='store_true',
                        help='Generate plots in each parser output directory')
    parser.add_argument('--seqdir', default=None,
                        help='Directory with AF3 input JSONs (needed for sequence_composition)')
    args = parser.parse_args()

    main(args.indir, args.outdir, args.parsers, args.plots, args.seqdir)
