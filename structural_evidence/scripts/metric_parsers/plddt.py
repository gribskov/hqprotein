#
# plddt.py - pLDDT analysis (Layers 1-3)
#
# Each layer adds columns to the per-protein feature table.
# analyze() runs all layers and writes results to outdir/

import csv
import numpy as np
import ruptures as rpt
from pathlib import Path
from scipy import stats as sp_stats


# --- Layer 1: basic statistics -------------------------------------------

def basic_stats(scores):
    """Summary statistics and threshold fractions from a pLDDT array."""
    n = len(scores)
    return {
        'plddt_mean':       np.mean(scores),
        'plddt_median':     np.median(scores),
        'plddt_std':        np.std(scores),
        'plddt_min':        np.min(scores),
        'plddt_max':        np.max(scores),
        'plddt_q10':        np.percentile(scores, 10),
        'plddt_q25':        np.percentile(scores, 25),
        'plddt_q75':        np.percentile(scores, 75),
        'plddt_q90':        np.percentile(scores, 90),
        'plddt_iqr':        np.percentile(scores, 75) - np.percentile(scores, 25),
        # standard AlphaFold confidence bins
        'frac_above90':     np.sum(scores >= 90) / n,
        'frac_above70':     np.sum(scores >= 70) / n,
        'frac_below50':     np.sum(scores < 50) / n,
        'frac_below30':     np.sum(scores < 30) / n,
    }


# --- Layer 2: distribution shape -----------------------------------------

def distribution_shape(scores):
    """Beta fit, entropy, moments."""
    features = {}

    # higher moments
    features['plddt_skewness'] = sp_stats.skew(scores)
    features['plddt_kurtosis'] = sp_stats.kurtosis(scores)

    # beta distribution fit (pLDDT is bounded [0, 100])
    scaled = np.clip(scores / 100.0, 0.001, 0.999)
    a, b, _, _ = sp_stats.beta.fit(scaled, floc=0, fscale=1)
    features['beta_a'] = a
    features['beta_b'] = b
    features['beta_mean'] = a / (a + b)
    features['beta_concentration'] = a + b

    # shannon entropy of binned distribution
    hist, _ = np.histogram(scores, bins=20, range=(0, 100))
    hist = hist / hist.sum()
    hist = hist[hist > 0]  # avoid log(0)
    features['plddt_entropy'] = -np.sum(hist * np.log2(hist))

    return features


# --- Layer 3: signal processing on per-residue trace ---------------------

def signal_features(scores):
    """Treat pLDDT as a 1D signal: runs, autocorrelation, dips."""
    features = {}
    n = len(scores)

    # -- run-length encoding at threshold 70 --
    binary = (scores >= 70).astype(int)
    transitions = np.diff(binary)
    features['n_transitions'] = int(np.sum(transitions != 0))

    run_starts = np.concatenate(([0], np.where(transitions != 0)[0] + 1))
    run_lengths = np.diff(np.concatenate((run_starts, [n])))
    run_values = binary[run_starts]

    high_runs = run_lengths[run_values == 1]
    low_runs = run_lengths[run_values == 0]

    features['max_high_run'] = int(np.max(high_runs)) if len(high_runs) > 0 else 0
    features['max_low_run'] = int(np.max(low_runs)) if len(low_runs) > 0 else 0
    features['mean_high_run'] = float(np.mean(high_runs)) if len(high_runs) > 0 else 0.0
    features['mean_low_run'] = float(np.mean(low_runs)) if len(low_runs) > 0 else 0.0
    features['frac_high_length'] = float(np.sum(high_runs)) / n if len(high_runs) > 0 else 0.0

    # -- autocorrelation (lags 1, 5, 10) --
    centered = scores - np.mean(scores)
    var = np.var(scores)
    for lag in [1, 5, 10]:
        if var > 0 and lag < n:
            features[f'autocorr_lag{lag}'] = np.sum(centered[:n-lag] * centered[lag:]) / (n * var)
        else:
            features[f'autocorr_lag{lag}'] = 0.0

    # -- PELT changepoint segmentation --
    # segments the trace into regions of similar pLDDT, then identifies
    # dips: segments whose mean is lower than both neighbors by >= 15 points
    if n >= 20:
        algo = rpt.Pelt(model='l2', min_size=10).fit(scores.reshape(-1, 1))
        bkps = algo.predict(pen=20)
        segments = np.split(scores, bkps[:-1])
        seg_means = [s.mean() for s in segments]

        features['n_segments'] = len(segments)

        # dips: segments lower than both neighbors by a meaningful amount
        n_dips = 0
        max_dip_drop = 0.0
        for i in range(1, len(seg_means) - 1):
            if seg_means[i] < seg_means[i-1] and seg_means[i] < seg_means[i+1]:
                drop = min(seg_means[i-1], seg_means[i+1]) - seg_means[i]
                if drop >= 15:
                    n_dips += 1
                    max_dip_drop = max(max_dip_drop, drop)

        features['n_dips'] = n_dips
        features['max_dip_drop'] = max_dip_drop

        # segment length stats
        seg_lengths = [len(s) for s in segments]
        features['mean_segment_len'] = float(np.mean(seg_lengths))
        features['std_segment_means'] = float(np.std(seg_means))
    else:
        features['n_segments'] = 1
        features['n_dips'] = 0
        features['max_dip_drop'] = 0.0
        features['mean_segment_len'] = float(n)
        features['std_segment_means'] = 0.0

    return features


# --- output ---------------------------------------------------------------

def write_feature_table(rows, outdir):
    """Write feature_table.tsv to outdir."""
    path = Path(outdir) / 'feature_table.tsv'
    fields = list(rows[0].keys())
    with open(path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


# --- plots ----------------------------------------------------------------

def make_plots(rows, dataset, outdir):
    """Generate L1-L3 plots in outdir/plots/."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns

    plotdir = Path(outdir) / 'plots'
    plotdir.mkdir(exist_ok=True)

    cats = ['Verified', 'Uncharacterized', 'Dubious']
    palette = {'Verified': '#2ecc71', 'Uncharacterized': '#f39c12', 'Dubious': '#e74c3c'}

    def _violin_strip(ax, rows, col, title):
        data = {cat: [float(r[col]) for r in rows if r['category'] == cat] for cat in cats}
        parts = ax.violinplot([data[c] for c in cats], positions=range(len(cats)),
                              showmedians=False, showextrema=False)
        for i, (pc, cat) in enumerate(zip(parts['bodies'], cats)):
            pc.set_facecolor(palette[cat])
            pc.set_alpha(0.4)
        for i, cat in enumerate(cats):
            jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(data[cat]))
            ax.scatter(i + jitter, data[cat], s=6, alpha=0.6,
                       color=palette[cat], edgecolors='none')
        ax.set_xticks(range(len(cats)))
        ax.set_xticklabels(cats)
        ax.set_title(title, fontsize=10)
        ax.tick_params(axis='x', rotation=15)

    # -- L1: basic stats (violin + strip) --
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    fig.suptitle('Layer 1 — Basic Statistics', fontsize=13, fontweight='bold')
    for ax, (col, title) in zip(axes.flat, [
        ('plddt_mean', 'Mean pLDDT'), ('plddt_median', 'Median pLDDT'),
        ('plddt_std', 'Std Dev'), ('frac_above90', 'Fraction >= 90'),
        ('frac_above70', 'Fraction >= 70'), ('frac_below50', 'Fraction < 50'),
    ]):
        _violin_strip(ax, rows, col, title)
    plt.tight_layout()
    plt.savefig(plotdir / 'L1_basic_stats.png', dpi=150)
    plt.close()

    # -- L2: distribution shape (violin + strip) --
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    fig.suptitle('Layer 2 — Distribution Shape', fontsize=13, fontweight='bold')
    for ax, (col, title) in zip(axes.flat, [
        ('plddt_skewness', 'Skewness'), ('plddt_kurtosis', 'Kurtosis'),
        ('beta_a', 'Beta a'), ('beta_b', 'Beta b'),
        ('beta_mean', 'Beta Mean — a/(a+b)'), ('plddt_entropy', 'Shannon Entropy'),
    ]):
        _violin_strip(ax, rows, col, title)
    plt.tight_layout()
    plt.savefig(plotdir / 'L2_distribution_shape.png', dpi=150)
    plt.close()

    # -- L2: pooled per-residue density by category --
    residues_by_cat = {c: [] for c in cats}
    for p in dataset:
        if p['plddt'] is not None:
            residues_by_cat[p['category']].extend(p['plddt'].tolist())

    fig, ax = plt.subplots(figsize=(8, 5))
    bins = np.linspace(0, 100, 100)
    for cat in cats:
        vals = np.array(residues_by_cat[cat])
        hist, edges = np.histogram(vals, bins=bins, density=True)
        centers = (edges[:-1] + edges[1:]) / 2
        ax.fill_between(centers, hist, alpha=0.25, color=palette[cat])
        ax.plot(centers, hist, linewidth=1.5, color=palette[cat],
                label=f'{cat} ({len(vals):,} residues)')
    ax.set_xlabel('pLDDT', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title('Per-Residue pLDDT Distribution by Category', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_xlim(0, 100)
    plt.tight_layout()
    plt.savefig(plotdir / 'L2_pooled_density.png', dpi=150)
    plt.close()

    # -- L3: signal features (3 panels) --
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(15, 5))
    gs = GridSpec(1, 3, width_ratios=[1, 1.3, 1], wspace=0.35)

    # panel 1: % proteins with dips
    ax1 = fig.add_subplot(gs[0])
    pcts = []
    for cat in cats:
        subset = [r for r in rows if r['category'] == cat]
        with_dips = sum(1 for r in subset if int(r['n_dips']) > 0)
        pcts.append(100 * with_dips / len(subset))
    bars = ax1.bar(cats, pcts, color=[palette[c] for c in cats], alpha=0.8, edgecolor='white')
    for bar, pct in zip(bars, pcts):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
                 f'{pct:.0f}%', ha='center', fontsize=10, fontweight='bold')
    ax1.set_ylabel('% of proteins', fontsize=10)
    ax1.set_title('Proteins with structural dips', fontsize=11, fontweight='bold')
    ax1.set_ylim(0, 60)
    ax1.tick_params(axis='x', rotation=15)

    # panel 2: scatter — std_segment_means vs plddt_mean
    ax2 = fig.add_subplot(gs[1])
    for cat in cats:
        subset = [r for r in rows if r['category'] == cat]
        x = [float(r['plddt_mean']) for r in subset]
        y = [float(r['std_segment_means']) for r in subset]
        ax2.scatter(x, y, c=palette[cat], alpha=0.5, s=20, label=cat, edgecolors='none')
    ax2.set_xlabel('Mean pLDDT', fontsize=10)
    ax2.set_ylabel('Std of segment means', fontsize=10)
    ax2.set_title('Structural complexity vs confidence', fontsize=11, fontweight='bold')
    ax2.legend(fontsize=8, loc='upper left')

    # panel 3: run lengths
    ax3 = fig.add_subplot(gs[2])
    x_pos = np.arange(len(cats))
    width = 0.35
    high_means = []
    low_means = []
    for cat in cats:
        subset = [r for r in rows if r['category'] == cat]
        high_means.append(np.mean([float(r['max_high_run']) for r in subset]))
        low_means.append(np.mean([float(r['max_low_run']) for r in subset]))
    ax3.bar(x_pos - width/2, high_means, width, color='#2ecc71', alpha=0.7, label='Longest confident run')
    ax3.bar(x_pos + width/2, low_means, width, color='#e74c3c', alpha=0.7, label='Longest low run')
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(cats, rotation=15)
    ax3.set_ylabel('Residues', fontsize=10)
    ax3.set_title('Run lengths by category', fontsize=11, fontweight='bold')
    ax3.legend(fontsize=8)

    plt.savefig(plotdir / 'L3_signal_features.png', dpi=150)
    plt.close()

    # -- Example traces: 2 per category --
    examples = [
        ('YAL062W', 'Verified — GDH3, clean high confidence'),
        ('YDL154W', 'Verified — 5 dips, complex domain architecture'),
        ('YML079W', 'Uncharacterized — high confidence, looks real'),
        ('YER158C', 'Uncharacterized — low confidence, looks dubious'),
        ('YAR073W', 'Dubious — pseudogene, fools pLDDT'),
        ('YAL037C-B', 'Dubious — flat, no structure'),
    ]
    traces = {p['gene_id']: p['plddt'] for p in dataset if p['plddt'] is not None}

    fig, axes = plt.subplots(len(examples), 1, figsize=(12, 2.5 * len(examples)))
    fig.suptitle('Example pLDDT Traces', fontsize=13, fontweight='bold')
    for ax, (pid, label) in zip(axes, examples):
        scores = traces[pid]
        cat = label.split(' — ')[0]
        mean_val = np.mean(scores)
        median_val = np.median(scores)
        ax.plot(scores, linewidth=0.8, color='#2c3e50')
        ax.axhline(y=mean_val, color='#3498db', linestyle='-', alpha=0.6, linewidth=1.0)
        ax.axhline(y=median_val, color='#e67e22', linestyle='--', alpha=0.6, linewidth=1.0)
        ax.fill_between(range(len(scores)), scores, alpha=0.2, color=palette.get(cat, '#999'))
        ax.text(len(scores)*0.98, mean_val+2, f'mean={mean_val:.1f}', fontsize=7,
                color='#3498db', ha='right', va='bottom')
        ax.text(len(scores)*0.98, median_val-2, f'med={median_val:.1f}', fontsize=7,
                color='#e67e22', ha='right', va='top')
        ax.set_ylim(0, 100)
        ax.set_ylabel('pLDDT')
        ax.set_title(f'{pid} — {label} ({len(scores)} aa)', fontsize=9)
    axes[-1].set_xlabel('Residue number')
    plt.tight_layout()
    plt.savefig(plotdir / 'example_traces.png', dpi=150)
    plt.close()


# --- main entry point ----------------------------------------------------

def analyze(dataset, outdir, plots=False):
    """Run all pLDDT layers and write results to outdir/."""
    rows = []
    for protein in dataset:
        scores = protein['plddt']
        if scores is None:
            continue

        row = {
            'gene_id': protein['gene_id'],
            'category': protein['category'],
            'length_aa': protein['length_aa'],
            'n_residues': len(scores),
        }
        row.update(basic_stats(scores))          # L1
        row.update(distribution_shape(scores))    # L2
        row.update(signal_features(scores))       # L3
        rows.append(row)

    write_feature_table(rows, outdir)

    if plots:
        make_plots(rows, dataset, outdir)
