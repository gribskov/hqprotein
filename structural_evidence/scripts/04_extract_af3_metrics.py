#!/usr/bin/env python3
#
# 04_extract_af3_metrics.py - Extract raw metrics from AlphaFold 3 output
# Alex Baena, Mar 2026
#
# Processes an AF3 output directory and writes per-protein metric files:
#   extracted_metrics_<tag>/
#     Q0075/
#       plddt.tsv                (residue_number, plddt)
#       ptm.tsv                  (ptm value)
#       fraction_disordered.tsv  (fraction value)
#       pae.tsv                  (NxN matrix, tab-delimited)
#       metadata.tsv             (gene_id, category, length, exons, etc.)
#     YAL046C/
#       ...
#

import sys
import csv
import json
import numpy as np
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def extract_id(filename):
    """Protein ID from AF3 filename. Q0075_model.cif -> Q0075"""
    name = Path(filename).stem
    if name.endswith('.cif'):
        name = Path(name).stem
    for suffix in ('_model', '_confidences', '_summary_confidences', '_data'):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


def find_proteins(af3_dir):
    """Locate top-ranked model directories in AF3 output.
    Returns dict: protein_id -> path to directory containing the output files.
    Skips seed-N_sample-N subdirs (individual samples)."""
    proteins = {}
    for cif in Path(af3_dir).glob('**/*_model.cif'):
        if 'seed-' in cif.parent.name:
            continue
        pid = extract_id(cif.name)
        proteins[pid] = cif.parent
    return proteins


def extract_plddt(cif_path):
    """Per-residue pLDDT from CIF B-factor column.
    Returns list of (residue_number, plddt)."""
    cif = MMCIF2Dict(str(cif_path))
    res_nums = cif['_atom_site.auth_seq_id']
    bfactors = cif['_atom_site.B_iso_or_equiv']

    # collapse atoms to residues — one pLDDT per residue
    seen = {}
    for rn, bf in zip(res_nums, bfactors):
        if rn not in seen:
            seen[rn] = float(bf)
    return sorted(seen.items(), key=lambda x: int(x[0]))


def extract_summary(json_path):
    """pTM and fraction_disordered from summary_confidences.json."""
    with open(json_path) as f:
        d = json.load(f)
    return d.get('ptm'), d.get('fraction_disordered')


def extract_pae(json_path):
    """PAE matrix from confidences.json. Returns NxN numpy array."""
    with open(json_path) as f:
        d = json.load(f)
    return np.array(d['pae'])


def load_metadata(metadata_path):
    """Load sample_metadata.tsv into dict keyed by gene_id."""
    meta = {}
    with open(metadata_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            meta[row['gene_id']] = row
    return meta


def process_protein(pid, protein_dir, out_dir, meta_row):
    """Extract all metrics for one protein and write to out_dir/pid/."""
    pdir = Path(out_dir) / pid
    pdir.mkdir(parents=True, exist_ok=True)

    # pLDDT
    cif_path = protein_dir / f'{pid}_model.cif'
    residues = extract_plddt(cif_path)
    with open(pdir / 'plddt.tsv', 'w') as f:
        f.write("residue_number\tplddt\n")
        for rn, val in residues:
            f.write(f"{rn}\t{val:.2f}\n")

    # pTM and fraction_disordered
    summary_path = protein_dir / f'{pid}_summary_confidences.json'
    if summary_path.exists():
        ptm, frac_dis = extract_summary(summary_path)

        with open(pdir / 'ptm.tsv', 'w') as f:
            f.write("ptm\n")
            f.write(f"{ptm if ptm is not None else 'NA'}\n")

        with open(pdir / 'fraction_disordered.tsv', 'w') as f:
            f.write("fraction_disordered\n")
            f.write(f"{frac_dis if frac_dis is not None else 'NA'}\n")

    # PAE matrix
    conf_path = protein_dir / f'{pid}_confidences.json'
    if conf_path.exists():
        pae = extract_pae(conf_path)
        np.savetxt(pdir / 'pae.tsv', pae, fmt='%.2f', delimiter='\t')

    # metadata
    cols = list(meta_row.keys())
    with open(pdir / 'metadata.tsv', 'w') as f:
        f.write('\t'.join(cols) + '\n')
        f.write('\t'.join(str(meta_row[c]) for c in cols) + '\n')


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: 04_extract_af3_metrics.py <af3_output_dir> <output_dir> "
              "<sample_metadata.tsv>", file=sys.stderr)
        sys.exit(1)

    af3_dir = sys.argv[1]
    out_dir = sys.argv[2]
    meta = load_metadata(sys.argv[3])

    proteins = find_proteins(af3_dir)
    n = 0
    for pid in sorted(proteins):
        try:
            process_protein(pid, proteins[pid], out_dir, meta[pid])
            n += 1
        except Exception as e:
            print(f"skip {pid}: {e}", file=sys.stderr)

    print(f"{n} / {len(proteins)} -> {out_dir}/", file=sys.stderr)
