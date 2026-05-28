#
# loader.py - Load extracted AF3 metrics into a unified dataset
#
# Reads per-protein TSV files from 04_extract_af3_metrics.py output.
# Metadata is read from each protein's own metadata.tsv (written by step 04).

import csv
import json
import numpy as np
from pathlib import Path


def load_plddt(protein_dir):
    """Read plddt.tsv → numpy array of per-residue scores."""
    path = protein_dir / 'plddt.tsv'
    if not path.exists():
        return None
    vals = []
    with open(path) as f:
        next(f)
        for line in f:
            _, score = line.strip().split('\t')
            vals.append(float(score))
    return np.array(vals)


def load_scalar(protein_dir, filename):
    """Read a single-value TSV (ptm.tsv, fraction_disordered.tsv)."""
    path = protein_dir / filename
    if not path.exists():
        return None
    with open(path) as f:
        next(f)
        val = f.readline().strip()
    return float(val) if val != 'NA' else None


def load_sequence(seqdir, protein_id):
    """Read amino acid sequence from AF3 input JSON."""
    jpath = seqdir / protein_id / f'{protein_id}.json'
    if not jpath.exists():
        return None
    with open(jpath) as f:
        data = json.load(f)
    for entry in data.get('sequences', []):
        if 'protein' in entry:
            return entry['protein']['sequence']
    return None


def load_protein_metadata(protein_dir):
    """Read per-protein metadata.tsv (written by step 04)."""
    path = protein_dir / 'metadata.tsv'
    if not path.exists():
        return None
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        return next(reader)


def validate_plddt(scores):
    """Sanity checks on a pLDDT array. Returns list of issues (empty = clean)."""
    issues = []
    if len(scores) == 0:
        issues.append('empty')
    if np.any(np.isnan(scores)):
        issues.append(f'{np.sum(np.isnan(scores))} NaNs')
    if np.any(scores < 0) or np.any(scores > 100):
        issues.append('values outside [0, 100]')
    return issues


def load_dataset(metrics_dir, seqdir=None):
    """Load all proteins into a list of dicts.

    Returns (dataset, problems) where problems is a list of (gene_id, issues).
    If seqdir is given, also loads amino acid sequences from AF3 input JSONs.
    """
    metrics_dir = Path(metrics_dir)
    seqdir = Path(seqdir) if seqdir else None
    protein_dirs = sorted([d for d in metrics_dir.iterdir() if d.is_dir()])

    dataset = []
    problems = []

    for pdir in protein_dirs:
        pid = pdir.name

        # metadata from step 04
        m = load_protein_metadata(pdir)
        if m is None:
            problems.append((pid, ['no metadata.tsv']))
            continue

        plddt = load_plddt(pdir)
        ptm = load_scalar(pdir, 'ptm.tsv')
        frac_dis = load_scalar(pdir, 'fraction_disordered.tsv')
        has_pae = (pdir / 'pae.tsv').exists()

        if plddt is not None:
            issues = validate_plddt(plddt)
            if issues:
                problems.append((pid, issues))

        seq = load_sequence(seqdir, pid) if seqdir else None

        dataset.append({
            'gene_id': pid,
            'category': m['category'],
            'length_aa': int(m['length_aa']),
            'n_exons': int(m['n_exons']),
            'has_introns': m['has_introns'],
            'length_bin': m['length_bin'],
            'chromosome': m['chromosome'],
            'gene_name': m.get('gene_name', ''),
            'plddt': plddt,
            'ptm': ptm,
            'fraction_disordered': frac_dis,
            'has_pae': has_pae,
            'sequence': seq,
        })

    return dataset, problems
