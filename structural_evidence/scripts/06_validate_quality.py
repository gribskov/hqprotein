#!/usr/bin/env python3
#
# 06_validate_quality.py - Classify gene models using extracted metric features
# Alex Baena, Apr 2026
#
# Reads feature tables from 05_parse_metrics.py, trains a logistic regression
# on Verified vs Dubious, evaluates with cross-validation, ranks features,
# and scores Uncharacterized proteins.
#
# Usage:
#   python3 06_validate_quality.py --indir parsed_metrics/ --outdir validation/
#   python3 06_validate_quality.py --indir parsed_metrics/ --outdir validation/ --plots

import sys
import csv
import numpy as np
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve, classification_report
from sklearn.calibration import CalibratedClassifierCV
from sklearn.decomposition import PCA


def discover_metrics(indir):
    """Find available metric subdirs (each containing a feature_table.tsv)."""
    indir = Path(indir)
    metrics = {}
    for subdir in sorted(indir.iterdir()):
        ft = subdir / 'feature_table.tsv'
        if subdir.is_dir() and ft.exists():
            metrics[subdir.name] = ft
    return metrics


def load_features(path):
    """Read a feature_table.tsv into rows (list of dicts)."""
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            rows.append(r)
    return rows


def extract_xy(rows):
    """Split rows into feature matrix X, labels y, and gene_ids.
    Only Verified and Dubious — Uncharacterized returned separately."""
    # identify feature columns (everything that's not metadata)
    skip = {'gene_id', 'category', 'length_aa', 'n_residues'}
    feature_cols = [k for k in rows[0].keys() if k not in skip]

    train_rows = [r for r in rows if r['category'] in ('Verified', 'Dubious')]
    unchar_rows = [r for r in rows if r['category'] == 'Uncharacterized']

    def to_matrix(subset):
        X = np.array([[float(r[c]) for c in feature_cols] for r in subset])
        ids = [r['gene_id'] for r in subset]
        return X, ids

    X_train, ids_train = to_matrix(train_rows)
    y_train = np.array([1 if r['category'] == 'Verified' else 0 for r in train_rows])

    X_unchar, ids_unchar = to_matrix(unchar_rows)

    return X_train, y_train, ids_train, X_unchar, ids_unchar, feature_cols


def classify(X, y, ids, X_unchar, ids_unchar, feature_cols):
    """Train logistic regression with cross-validation. Returns results dict."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    X_unchar_scaled = scaler.transform(X_unchar)

    # cross-validated predictions (each sample predicted when held out)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    model = LogisticRegression(max_iter=1000, random_state=42)
    cv_probs = cross_val_predict(model, X_scaled, y, cv=cv, method='predict_proba')[:, 1]
    cv_auc = roc_auc_score(y, cv_probs)

    # train final model on all data
    model.fit(X_scaled, y)
    coefs = model.coef_[0]

    # feature importance (absolute coefficient, standardized features)
    importance = sorted(zip(feature_cols, coefs), key=lambda x: -abs(x[1]))

    # score uncharacterized
    unchar_probs = model.predict_proba(X_unchar_scaled)[:, 1]

    # ROC curve data
    fpr, tpr, thresholds = roc_curve(y, cv_probs)

    # bootstrap 95% CI for AUC
    rng = np.random.default_rng(42)
    n_boot = 2000
    boot_aucs = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.choice(len(y), size=len(y), replace=True)
        if len(np.unique(y[idx])) < 2:
            boot_aucs[i] = np.nan
            continue
        boot_aucs[i] = roc_auc_score(y[idx], cv_probs[idx])
    boot_aucs = boot_aucs[~np.isnan(boot_aucs)]
    auc_ci_lo, auc_ci_hi = np.percentile(boot_aucs, [2.5, 97.5])

    # PCA on scaled training + uncharacterized data
    X_all = np.vstack([X_scaled, X_unchar_scaled])
    cats_all = (['Verified' if yi == 1 else 'Dubious' for yi in y]
                + ['Uncharacterized'] * len(ids_unchar))
    pca = PCA(n_components=min(5, X_all.shape[1]))
    X_pca = pca.fit_transform(X_all)
    loadings = pca.components_  # (n_components, n_features)

    return {
        'cv_auc': cv_auc,
        'auc_ci': (auc_ci_lo, auc_ci_hi),
        'cv_probs': cv_probs,
        'y': y,
        'ids_train': ids,
        'importance': importance,
        'fpr': fpr,
        'tpr': tpr,
        'thresholds': thresholds,
        'unchar_probs': unchar_probs,
        'ids_unchar': ids_unchar,
        'X_pca': X_pca,
        'pca_cats': cats_all,
        'pca_loadings': loadings,
        'pca_variance': pca.explained_variance_ratio_,
        'feature_cols': feature_cols,
    }


def write_results(results, feature_cols, outdir):
    """Write classification outputs to outdir/."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # feature importance
    with open(outdir / 'feature_importance.tsv', 'w') as f:
        f.write('feature\tcoefficient\tabs_coefficient\n')
        for feat, coef in results['importance']:
            f.write(f'{feat}\t{coef:.4f}\t{abs(coef):.4f}\n')

    # cross-validation predictions (Verified + Dubious)
    with open(outdir / 'cv_predictions.tsv', 'w') as f:
        f.write('gene_id\ttrue_label\tpredicted_prob\n')
        labels = {1: 'Verified', 0: 'Dubious'}
        for gid, true, prob in zip(results['ids_train'], results['y'], results['cv_probs']):
            f.write(f'{gid}\t{labels[true]}\t{prob:.4f}\n')

    # uncharacterized predictions
    with open(outdir / 'uncharacterized_predictions.tsv', 'w') as f:
        f.write('gene_id\tpredicted_prob\n')
        for gid, prob in sorted(zip(results['ids_unchar'], results['unchar_probs']),
                                 key=lambda x: -x[1]):
            f.write(f'{gid}\t{prob:.4f}\n')

    # summary
    with open(outdir / 'summary.txt', 'w') as f:
        ci_lo, ci_hi = results['auc_ci']
        f.write(f'Cross-validated AUC: {results["cv_auc"]:.4f}  (95% CI: {ci_lo:.4f} - {ci_hi:.4f})\n')
        f.write(f'Training set: {sum(results["y"])} Verified, {len(results["y"]) - sum(results["y"])} Dubious\n')
        f.write(f'Uncharacterized scored: {len(results["ids_unchar"])}\n\n')
        f.write('Top 10 features:\n')
        for feat, coef in results['importance'][:10]:
            direction = 'Verified' if coef > 0 else 'Dubious'
            f.write(f'  {feat:<30} {coef:>8.4f}  ({direction})\n')

        # uncharacterized breakdown
        probs = results['unchar_probs']
        f.write(f'\nUncharacterized predictions:\n')
        f.write(f'  P > 0.9 (likely real):     {sum(probs > 0.9)}\n')
        f.write(f'  P 0.5-0.9:                 {sum((probs > 0.5) & (probs <= 0.9))}\n')
        f.write(f'  P 0.1-0.5:                 {sum((probs >= 0.1) & (probs <= 0.5))}\n')
        f.write(f'  P < 0.1 (likely dubious):  {sum(probs < 0.1)}\n')


def make_plots(results, outdir):
    """Generate ROC curve, feature importance, and prediction distribution plots."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plotdir = Path(outdir) / 'plots'
    plotdir.mkdir(exist_ok=True)

    # -- ROC curve with optimal threshold --
    fig, ax = plt.subplots(figsize=(6, 6))
    fpr, tpr = results['fpr'], results['tpr']
    ci_lo, ci_hi = results['auc_ci']
    ax.plot(fpr, tpr, color='#2c3e50', linewidth=2,
            label=f'AUC = {results["cv_auc"]:.3f} ({ci_lo:.3f}-{ci_hi:.3f})')
    ax.plot([0, 1], [0, 1], color='#bdc3c7', linestyle='--', linewidth=1)

    # optimal point: maximizes TPR - FPR (Youden's J)
    j_scores = tpr - fpr
    best_idx = np.argmax(j_scores)
    best_fpr, best_tpr = fpr[best_idx], tpr[best_idx]
    best_threshold = results['thresholds'][best_idx]
    ax.scatter(best_fpr, best_tpr, color='#e74c3c', s=80, zorder=5)
    ax.annotate(f'threshold={best_threshold:.2f}\nTPR={best_tpr:.2f}, FPR={best_fpr:.2f}',
                xy=(best_fpr, best_tpr), xytext=(best_fpr + 0.15, best_tpr - 0.15),
                fontsize=9, arrowprops=dict(arrowstyle='->', color='#e74c3c'),
                color='#e74c3c')

    ax.set_xlabel('False Positive Rate', fontsize=11)
    ax.set_ylabel('True Positive Rate', fontsize=11)
    ax.set_title('ROC Curve — Verified vs Dubious (5-fold CV)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=11, loc='lower right')
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    plt.tight_layout()
    plt.savefig(plotdir / 'roc_curve.png', dpi=150)
    plt.close()

    # -- Feature importance (top 15) --
    top = results['importance'][:15]
    feats = [f for f, _ in top][::-1]
    coefs = [c for _, c in top][::-1]
    colors = ['#2ecc71' if c > 0 else '#e74c3c' for c in coefs]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(feats, coefs, color=colors, alpha=0.8)
    ax.set_xlabel('Coefficient (positive = Verified, negative = Dubious)', fontsize=10)
    ax.set_title('Feature Importance — Top 15', fontsize=12, fontweight='bold')
    ax.axvline(x=0, color='#2c3e50', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(plotdir / 'feature_importance.png', dpi=150)
    plt.close()

    # -- Prediction distribution for Uncharacterized --
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(results['unchar_probs'], bins=20, range=(0, 1), color='#f39c12',
            alpha=0.7, edgecolor='white')
    ax.axvline(x=0.5, color='#2c3e50', linestyle='--', linewidth=1)
    ax.set_xlabel('P(Verified)', fontsize=11)
    ax.set_ylabel('Number of proteins', fontsize=11)
    ax.set_title('Predicted Probabilities — Uncharacterized Proteins', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(plotdir / 'uncharacterized_distribution.png', dpi=150)
    plt.close()

    # -- PCA biplot --
    palette = {'Verified': '#2ecc71', 'Uncharacterized': '#f39c12', 'Dubious': '#e74c3c'}
    X_pca = results['X_pca']
    cats = results['pca_cats']
    loadings = results['pca_loadings']
    var_exp = results['pca_variance']
    feature_names = results['feature_cols']

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # left panel: scatter PC1 vs PC2
    ax = axes[0]
    for cat in ['Verified', 'Uncharacterized', 'Dubious']:
        mask = [c == cat for c in cats]
        ax.scatter(X_pca[mask, 0], X_pca[mask, 1], c=palette[cat],
                   alpha=0.5, s=25, label=cat, edgecolors='none')
    ax.set_xlabel(f'PC1 ({var_exp[0]*100:.1f}% variance)', fontsize=11)
    ax.set_ylabel(f'PC2 ({var_exp[1]*100:.1f}% variance)', fontsize=11)
    ax.set_title('PCA — Proteins in Feature Space', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)

    # right panel: loading arrows (biplot)
    ax = axes[1]
    n_arrows = min(12, len(feature_names))
    # rank by total loading magnitude across PC1+PC2
    mag = np.sqrt(loadings[0]**2 + loadings[1]**2)
    top_idx = np.argsort(mag)[::-1][:n_arrows]

    for i in top_idx:
        ax.annotate('', xy=(loadings[0][i], loadings[1][i]), xytext=(0, 0),
                     arrowprops=dict(arrowstyle='->', color='#2c3e50', lw=1.3))
        ax.text(loadings[0][i] * 1.08, loadings[1][i] * 1.08,
                feature_names[i], fontsize=7, ha='center', va='center',
                color='#2c3e50')
    ax.set_xlabel(f'PC1 loading ({var_exp[0]*100:.1f}%)', fontsize=11)
    ax.set_ylabel(f'PC2 loading ({var_exp[1]*100:.1f}%)', fontsize=11)
    ax.set_title('Feature Loadings — Top 12', fontsize=12, fontweight='bold')
    ax.axhline(0, color='#bdc3c7', linewidth=0.5)
    ax.axvline(0, color='#bdc3c7', linewidth=0.5)
    lim = max(abs(loadings[0].max()), abs(loadings[1].max())) * 1.3
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(plotdir / 'pca_biplot.png', dpi=150)
    plt.close()

    # -- PCA variance explained bar chart --
    fig, ax = plt.subplots(figsize=(6, 4))
    n_pcs = len(var_exp)
    cumulative = np.cumsum(var_exp)
    ax.bar(range(1, n_pcs+1), var_exp * 100, color='#3498db', alpha=0.7, label='Individual')
    ax.step(range(1, n_pcs+1), cumulative * 100, where='mid', color='#2c3e50',
            linewidth=2, label='Cumulative')
    ax.set_xlabel('Principal Component', fontsize=11)
    ax.set_ylabel('Variance Explained (%)', fontsize=11)
    ax.set_title('PCA — Variance Explained', fontsize=12, fontweight='bold')
    ax.set_xticks(range(1, n_pcs+1))
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(plotdir / 'pca_variance.png', dpi=150)
    plt.close()


def main(indir, outdir, plots=False):
    metrics = discover_metrics(indir)
    if not metrics:
        print('No feature tables found', file=sys.stderr)
        sys.exit(1)

    # run classification per metric
    for name, path in metrics.items():
        rows = load_features(path)
        X, y, ids, X_unc, ids_unc, feature_cols = extract_xy(rows)
        results = classify(X, y, ids, X_unc, ids_unc, feature_cols)

        metric_outdir = Path(outdir) / name
        write_results(results, feature_cols, metric_outdir)

        if plots:
            make_plots(results, metric_outdir)

        ci = results['auc_ci']
        print(f'{name}: AUC={results["cv_auc"]:.3f} (95% CI: {ci[0]:.3f}-{ci[1]:.3f}) | '
              f'{sum(y)} Verified vs {len(y)-sum(y)} Dubious | '
              f'{len(ids_unc)} Uncharacterized scored')

    # combined classification if multiple metrics exist
    if len(metrics) > 1:
        all_rows = {}
        for name, path in metrics.items():
            for r in load_features(path):
                gid = r['gene_id']
                if gid not in all_rows:
                    all_rows[gid] = {'gene_id': gid, 'category': r['category'],
                                     'length_aa': r['length_aa'], 'n_residues': r.get('n_residues', '')}
                skip = {'gene_id', 'category', 'length_aa', 'n_residues'}
                for k, v in r.items():
                    if k not in skip:
                        all_rows[gid][f'{name}_{k}'] = v

        combined = list(all_rows.values())
        X, y, ids, X_unc, ids_unc, feature_cols = extract_xy(combined)
        results = classify(X, y, ids, X_unc, ids_unc, feature_cols)

        combined_outdir = Path(outdir) / 'combined'
        write_results(results, feature_cols, combined_outdir)
        if plots:
            make_plots(results, combined_outdir)
        ci = results['auc_ci']
        print(f'combined: AUC={results["cv_auc"]:.3f} (95% CI: {ci[0]:.3f}-{ci[1]:.3f})')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Classify gene models using extracted metric features')
    parser.add_argument('--indir', required=True,
                        help='Directory with parsed metrics (from step 05)')
    parser.add_argument('--outdir', required=True,
                        help='Output directory for classification results')
    parser.add_argument('--plots', action='store_true',
                        help='Generate ROC curve, feature importance, and prediction plots')
    args = parser.parse_args()

    main(args.indir, args.outdir, args.plots)
