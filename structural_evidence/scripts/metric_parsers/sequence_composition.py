#
# sequence_composition.py - Amino acid composition features
#
# Sequence-level features independent of structure prediction:
# charge, pI, hydrophobicity, composition bias. Dubious ORFs often
# have skewed composition because they weren't selected for foldability.

import csv
import numpy as np
from pathlib import Path


# --- constants -----------------------------------------------------------

# Kyte-Doolittle hydrophobicity scale
KD = {
    'A':  1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C':  2.5,
    'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I':  4.5,
    'L':  3.8, 'K': -3.9, 'M':  1.9, 'F':  2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V':  4.2,
}

# residue pKa values for pI calculation (Bjellqvist, EMBOSS-style)
PKA_SIDE = {'C': 8.5, 'D': 3.86, 'E': 4.25, 'H': 6.0,
            'K': 10.53, 'R': 12.48, 'Y': 10.07}
PKA_NTERM = 9.69
PKA_CTERM = 2.34

POSITIVE = set('KR')         # always positive at pH 7
NEGATIVE = set('DE')         # always negative at pH 7
CHARGED = POSITIVE | NEGATIVE | {'H'}
HYDROPHOBIC = set('AILMFVWY')
POLAR = set('STNQ')
AROMATIC = set('FWY')


# --- charge / pI ---------------------------------------------------------

def net_charge_at_ph(seq, pH=7.0):
    """Approximate net charge at given pH using Henderson-Hasselbalch."""
    charge = 0.0
    # N-terminus (+)
    charge += 1.0 / (1.0 + 10 ** (pH - PKA_NTERM))
    # C-terminus (-)
    charge -= 1.0 / (1.0 + 10 ** (PKA_CTERM - pH))
    for aa in seq:
        pKa = PKA_SIDE.get(aa)
        if pKa is None:
            continue
        if aa in ('K', 'R', 'H'):
            charge += 1.0 / (1.0 + 10 ** (pH - pKa))
        else:  # D, E, C, Y
            charge -= 1.0 / (1.0 + 10 ** (pKa - pH))
    return charge


def isoelectric_point(seq):
    """Bisection search for pH where net charge = 0."""
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo + hi) / 2
        c = net_charge_at_ph(seq, mid)
        if c > 0:
            lo = mid
        else:
            hi = mid
        if abs(c) < 1e-4:
            break
    return mid


# --- feature layers ------------------------------------------------------

def composition_features(seq):
    """Amino acid fractions and simple class summaries."""
    n = len(seq)
    counts = {aa: seq.count(aa) for aa in 'ACDEFGHIKLMNPQRSTVWY'}
    features = {f'frac_{aa}': counts[aa] / n for aa in counts}

    features['frac_hydrophobic'] = sum(counts[a] for a in HYDROPHOBIC) / n
    features['frac_polar'] = sum(counts[a] for a in POLAR) / n
    features['frac_charged'] = sum(counts[a] for a in CHARGED) / n
    features['frac_aromatic'] = sum(counts[a] for a in AROMATIC) / n
    features['frac_positive'] = sum(counts[a] for a in POSITIVE) / n
    features['frac_negative'] = sum(counts[a] for a in NEGATIVE) / n

    # Shannon entropy of composition — low entropy means skewed toward few AAs
    freqs = np.array([counts[a] / n for a in counts if counts[a] > 0])
    features['composition_entropy'] = float(-np.sum(freqs * np.log2(freqs)))

    return features


def charge_features(seq):
    """Charge and pI features."""
    n = len(seq)
    n_pos = sum(1 for a in seq if a in POSITIVE)
    n_neg = sum(1 for a in seq if a in NEGATIVE)

    charge7 = net_charge_at_ph(seq, 7.0)
    return {
        'net_charge_ph7':        charge7,
        'net_charge_per_res':    charge7 / n,
        'abs_charge_per_res':    abs(charge7) / n,
        'charge_asymmetry':      (n_pos - n_neg) / n,
        'pI':                    isoelectric_point(seq),
    }


def hydrophobicity_features(seq):
    """Kyte-Doolittle hydrophobicity summaries + TM-helix candidate runs.

    A sliding window >= 19 residues with mean KD >= 1.6 is a candidate
    transmembrane segment — used to flag membrane proteins so the model
    can learn they are real despite high pLDDT variance.
    """
    n = len(seq)
    kd = np.array([KD.get(a, 0.0) for a in seq])

    features = {
        'kd_mean': float(np.mean(kd)),
        'kd_std':  float(np.std(kd)),
        'kd_max':  float(np.max(kd)),
        'kd_min':  float(np.min(kd)),
    }

    # sliding window 19 residues (standard TM helix length)
    window = 19
    if n >= window:
        means = np.convolve(kd, np.ones(window) / window, mode='valid')
        features['kd_window_max'] = float(np.max(means))
        features['n_tm_candidate'] = int(np.sum(means >= 1.6))
        # collapse overlapping windows into segments
        above = means >= 1.6
        if above.any():
            transitions = np.diff(above.astype(int))
            starts = int(np.sum(transitions == 1)) + (1 if above[0] else 0)
            features['n_tm_segments'] = starts
        else:
            features['n_tm_segments'] = 0
    else:
        features['kd_window_max'] = float(np.mean(kd))
        features['n_tm_candidate'] = 0
        features['n_tm_segments'] = 0

    return features


# --- output --------------------------------------------------------------

def write_feature_table(rows, outdir):
    path = Path(outdir) / 'feature_table.tsv'
    fields = list(rows[0].keys())
    with open(path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


# --- plots ---------------------------------------------------------------

def make_plots(rows, outdir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plotdir = Path(outdir) / 'plots'
    plotdir.mkdir(exist_ok=True)

    cats = ['Verified', 'Uncharacterized', 'Dubious']
    palette = {'Verified': '#2ecc71', 'Uncharacterized': '#f39c12', 'Dubious': '#e74c3c'}

    def _violin_strip(ax, rows, col, title):
        data = {cat: [float(r[col]) for r in rows if r['category'] == cat] for cat in cats}
        parts = ax.violinplot([data[c] for c in cats], positions=range(len(cats)),
                              showmedians=False, showextrema=False)
        for pc, cat in zip(parts['bodies'], cats):
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

    # charge panel
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    fig.suptitle('Sequence Composition — Charge & pI', fontsize=13, fontweight='bold')
    for ax, (col, title) in zip(axes.flat, [
        ('net_charge_per_res', 'Net charge / residue'),
        ('abs_charge_per_res', '|Net charge| / residue'),
        ('charge_asymmetry',   '(K+R − D−E) / length'),
        ('pI',                 'Isoelectric point'),
        ('frac_charged',       'Fraction charged'),
        ('composition_entropy','Composition entropy'),
    ]):
        _violin_strip(ax, rows, col, title)
    plt.tight_layout()
    plt.savefig(plotdir / 'charge_pi.png', dpi=150)
    plt.close()

    # hydrophobicity panel
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    fig.suptitle('Sequence Composition — Hydrophobicity', fontsize=13, fontweight='bold')
    for ax, (col, title) in zip(axes.flat, [
        ('kd_mean',        'Mean Kyte-Doolittle'),
        ('kd_std',         'Std Kyte-Doolittle'),
        ('kd_window_max',  'Max 19-residue window'),
        ('n_tm_segments',  'Candidate TM segments'),
        ('frac_hydrophobic','Fraction hydrophobic'),
        ('frac_aromatic',  'Fraction aromatic'),
    ]):
        _violin_strip(ax, rows, col, title)
    plt.tight_layout()
    plt.savefig(plotdir / 'hydrophobicity.png', dpi=150)
    plt.close()


# --- main entry point ----------------------------------------------------

def analyze(dataset, outdir, plots=False):
    """Run sequence composition analysis and write results to outdir/."""
    rows = []
    skipped = 0
    for protein in dataset:
        seq = protein.get('sequence')
        if not seq:
            skipped += 1
            continue

        row = {
            'gene_id': protein['gene_id'],
            'category': protein['category'],
            'length_aa': protein['length_aa'],
            'n_residues': len(seq),
        }
        row.update(composition_features(seq))
        row.update(charge_features(seq))
        row.update(hydrophobicity_features(seq))
        rows.append(row)

    if not rows:
        print('  sequence_composition: no sequences available (pass --seqdir)')
        return

    write_feature_table(rows, outdir)
    if plots:
        make_plots(rows, outdir)
    if skipped:
        print(f'  sequence_composition: skipped {skipped} proteins with no sequence')
