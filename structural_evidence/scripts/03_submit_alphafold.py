#!/usr/bin/env python3
#
# 03_submit_alphafold.py - Sample proteins and generate AlphaFold 3 input JSONs
# Alex Baena, Feb 2026
#
# Builds a stratified dev set of 330 proteins from yeast (S. cerevisiae):
#   - 110 Verified, 110 Uncharacterized, 110 Dubious
#   - 20 per length bin per category
#   - 3 intron-containing genes forced per category
#   - Random seed for reproducibility
#
# Reads:
#   1. Protein FASTA (orf_trans_all) — sequences + basic metadata from headers
#   2. GFF annotation — exon counts per gene
#
# Writes:
#   1. One AF3 input JSON per protein in output_dir/
#   2. A sample_metadata.tsv with all gene info for downstream analysis
#

import sys
import gzip
import json
import random
from pathlib import Path
from collections import defaultdict


# Length bins: (label, min_aa, max_aa)
BINS = [
    ('<=100',   0,   100),
    ('101-200', 101, 200),
    ('201-300', 201, 300),
    ('301-500', 301, 500),
    ('501-800', 501, 800),
    ('>800',    801, 99999),
]

# samples per bin per category
BIN_QUOTA = {
    '<=100':   20,
    '101-200': 20,
    '201-300': 20,
    '301-500': 20,
    '501-800': 20,
    '>800':    20,
}

# Force this many intron-containing genes per category.
# Dubious only has 3 intron genes, so 3 is the max.
FORCED_INTRON = 3

SEED = 42


def assign_bin(length):
    """Map a protein length to its bin label."""
    for label, lo, hi in BINS:
        if lo <= length <= hi:
            return label
    return '>800'


def parse_fasta(fasta_path):
    """Read protein sequences and header metadata from SGD FASTA.

    SGD FASTA headers look like:
    >YAL068C PAU8 SGDID:S000002142, Chr I from 2169-1807, ... Verified ORF, "..."

    The gene_id is the first field, gene_name is the second (same as gene_id
    if the gene has no common name), and the category is embedded as
    'Verified ORF', 'Dubious ORF', or 'Uncharacterized ORF'.
    """
    genes = {}
    current_id = None
    current_seq = []

    # handle both gzipped and plain FASTA
    opener = gzip.open if str(fasta_path).endswith('.gz') else open
    mode = 'rt' if str(fasta_path).endswith('.gz') else 'r'

    with opener(fasta_path, mode) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # save previous gene
                if current_id is not None:
                    seq = ''.join(current_seq).rstrip('*')
                    genes[current_id]['sequence'] = seq
                    genes[current_id]['length_aa'] = len(seq)

                # parse header
                header = line[1:]  # drop '>'
                fields = header.split()
                gene_id = fields[0]

                # gene_name: second field, but if it equals gene_id there's
                # no common name. Also skip if it looks like 'SGDID:...'
                gene_name = fields[1] if len(fields) > 1 else gene_id
                if gene_name.startswith('SGDID:') or gene_name == gene_id:
                    gene_name = ''

                # category from header text
                if 'Verified ORF' in header:
                    category = 'Verified'
                elif 'Uncharacterized ORF' in header:
                    category = 'Uncharacterized'
                elif 'Dubious ORF' in header:
                    category = 'Dubious'
                else:
                    category = 'Unknown'

                # chromosome from "Chr X" in header
                chrom = ''
                if ', Chr ' in header:
                    chrom_part = header.split(', Chr ')[1]
                    chrom = 'chr' + chrom_part.split()[0]

                genes[gene_id] = {
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'category': category,
                    'chromosome': chrom,
                }
                current_id = gene_id
                current_seq = []
            else:
                current_seq.append(line)

    # don't forget the last gene
    if current_id is not None:
        seq = ''.join(current_seq).rstrip('*')
        genes[current_id]['sequence'] = seq
        genes[current_id]['length_aa'] = len(seq)

    return genes


def count_exons_from_gff(gff_path):
    """Count CDS features per gene from GFF — each CDS is one exon.

    Returns dict: gene_id -> number of exons.
    Single-exon genes have count=1, multi-exon genes have count>1.
    """
    exon_counts = defaultdict(int)

    opener = gzip.open if str(gff_path).endswith('.gz') else open
    mode = 'rt' if str(gff_path).endswith('.gz') else 'r'

    with opener(gff_path, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'CDS':
                continue

            # CDS lines have Parent=YAL069W_mRNA — strip _mRNA to get gene_id
            attrs = parts[8]
            for attr in attrs.split(';'):
                if attr.startswith('Parent='):
                    parent = attr[7:]
                    # some parents have _mRNA suffix, some have _id001 etc
                    # normalize to gene_id by stripping _mRNA
                    if '_mRNA' in parent:
                        gene_id = parent.split('_mRNA')[0]
                    elif '_id' in parent:
                        # genes with multiple transcript isoforms:
                        # Parent=YPL198W_id003,YPL198W_id001
                        # take first one, strip _id part
                        gene_id = parent.split(',')[0]
                        gene_id = gene_id.split('_id')[0]
                    else:
                        gene_id = parent
                    exon_counts[gene_id] += 1

    return dict(exon_counts)


def sample_proteins(genes, seed=SEED):
    """Stratified sampling: equal per length bin, force intron genes.

    Strategy:
    1. For each category, find intron-containing genes
    2. Randomly pick FORCED_INTRON of them, place in their natural bins
    3. Fill remaining bin slots with random non-intron genes
    """
    random.seed(seed)

    # organize genes by category and bin
    by_cat_bin = defaultdict(lambda: defaultdict(list))
    intron_genes = defaultdict(list)

    for gid, info in genes.items():
        cat = info['category']
        if cat not in ('Verified', 'Uncharacterized', 'Dubious'):
            continue
        length_bin = assign_bin(info['length_aa'])
        info['length_bin'] = length_bin
        by_cat_bin[cat][length_bin].append(gid)

        if info.get('n_exons', 1) > 1:
            intron_genes[cat].append(gid)

    selected = []

    for cat in ('Verified', 'Uncharacterized', 'Dubious'):
        # step 1: pick forced intron genes
        available_intron = intron_genes[cat][:]
        random.shuffle(available_intron)
        forced = available_intron[:FORCED_INTRON]

        # fill each bin
        for label, _, _ in BINS:
            quota = BIN_QUOTA[label]
            pool = by_cat_bin[cat][label]

            # forced intron genes that fall in this bin
            forced_in_bin = [g for g in forced if genes[g]['length_bin'] == label]

            # remaining slots to fill randomly (excluding forced ones)
            remaining = quota - len(forced_in_bin)
            non_forced = [g for g in pool if g not in forced]
            random.shuffle(non_forced)

            if remaining > len(non_forced):
                pick = non_forced
            else:
                pick = non_forced[:remaining]

            selected.extend(forced_in_bin + pick)

    return selected


def write_af3_json(gene_info, output_dir):
    """Write one AlphaFold 3 input JSON for a single protein.

    Format follows AF3 v4 spec:
    - name: gene_id (used for output directory naming)
    - modelSeeds: [1] (single prediction, reproducible)
    - sequences: one protein chain with id "A"
    - description: metadata string (AF3 ignores it, useful for tracking)
    """
    gid = gene_info['gene_id']
    desc = (f"{gene_info['category']} | {gene_info['n_exons']} exons | "
            f"{gene_info['length_aa']} AA | {gene_info['chromosome']}")

    af3_input = {
        'name': gid,
        'modelSeeds': [1],
        'sequences': [
            {
                'protein': {
                    'id': 'A',
                    'sequence': gene_info['sequence'],
                    'description': desc,
                }
            }
        ],
        'dialect': 'alphafold3',
        'version': 4,
    }

    out_path = Path(output_dir) / f'{gid}.json'
    with open(out_path, 'w') as f:
        json.dump(af3_input, f, indent=2)


def write_metadata(genes, selected_ids, output_path):
    """Write the master metadata TSV for all sampled proteins.

    This file links gene_id to category, length, exon info, etc.
    We join this with AF3 results during analysis.
    """
    with open(output_path, 'w') as f:
        f.write('gene_id\tcategory\tlength_aa\tn_exons\t'
                'has_introns\tlength_bin\tchromosome\tgene_name\n')
        for gid in sorted(selected_ids):
            g = genes[gid]
            has_introns = 'yes' if g.get('n_exons', 1) > 1 else 'no'
            f.write(f"{g['gene_id']}\t{g['category']}\t{g['length_aa']}\t"
                    f"{g.get('n_exons', 1)}\t{has_introns}\t"
                    f"{g.get('length_bin', '')}\t{g['chromosome']}\t"
                    f"{g['gene_name']}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("Usage: 03_submit_alphafold.py <proteins.fasta[.gz]> "
              "<annotation.gff[.gz]> <output_dir> <metadata.tsv>",
              file=sys.stderr)
        sys.exit(1)

    fasta_path = sys.argv[1]
    gff_path = sys.argv[2]
    output_dir = sys.argv[3]
    metadata_path = sys.argv[4]

    # make output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    genes = parse_fasta(fasta_path)

    exon_counts = count_exons_from_gff(gff_path)
    for gid in genes:
        genes[gid]['n_exons'] = exon_counts.get(gid, 1)

    selected = sample_proteins(genes)

    for gid in selected:
        write_af3_json(genes[gid], output_dir)

    write_metadata(genes, selected, metadata_path)

    # one-line summary: how many sampled out of how many total
    print(f"{len(selected)} / {len(genes)} -> {output_dir}/", file=sys.stderr)
