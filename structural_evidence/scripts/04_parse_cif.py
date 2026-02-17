#!/usr/bin/env python3
#
# 04_parse_cif.py - Extract per-residue pLDDT from AlphaFold CIF files
# Alex Baena, Feb 2026
#

import sys
import tarfile       # read CIF files directly from AlphaFold tar archives
import gzip          # CIF files inside the tar are gzipped (.cif.gz)
import io            # convert binary tar streams to text for BioPython
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict  # BioPython CIF reader


def parse_single_cif(cif_path_or_fileobj):
    """Extract per-residue pLDDT from one CIF file.
    Returns list of (residue_number, plddt) tuples."""

    # MMCIF2Dict parses CIF into a dict: key = column name, value = list of strings.
    # A typical AlphaFold CIF has ~236 keys. We only need 2 of them.
    # It accepts either a file path (string) or a file object (open stream),
    # so we check which one we got — file object when reading from tar,
    # string when reading a single .cif from disk.
    if hasattr(cif_path_or_fileobj, 'read'):
        cif = MMCIF2Dict(cif_path_or_fileobj)
    else:
        cif = MMCIF2Dict(str(cif_path_or_fileobj))

    # auth_seq_id = residue number (1, 2, 3, ... up to protein length).
    # Each residue appears multiple times because each amino acid is made
    # of ~5-8 atoms (C, N, O, etc), so a 145-residue protein has ~1122 rows.
    residue_numbers = cif['_atom_site.auth_seq_id']

    # B_iso_or_equiv = the B-factor column. In experimental crystallography
    # this measures atomic vibration, but AlphaFold repurposes it to store
    # pLDDT (0-100 confidence). This is universal — works on CIF from AF2,
    # AF3, Boltz, ColabFold, any structure predictor.
    bfactors = cif['_atom_site.B_iso_or_equiv']

    # Collapse atoms to residues: all atoms in the same residue share the
    # same pLDDT (AlphaFold scores residues, not individual atoms), so we
    # just keep the first atom we see for each residue number.
    plddt = {}
    for res_num, bfactor in zip(residue_numbers, bfactors):
        if res_num not in plddt:
            plddt[res_num] = float(bfactor)

    # Return sorted by residue number (1, 2, 3, ...) for sequential output
    return sorted(plddt.items(), key=lambda x: int(x[0]))


def extract_uniprot_from_filename(filename):
    """Pull UniProt accession from AlphaFold filename.
    e.g. AF-A0A023PYF4-F1-model_v6.cif.gz -> A0A023PYF4
    We need this to link pLDDT scores back to gene annotations."""
    name = Path(filename).stem       # strip .cif or .gz
    if name.endswith('.cif'):         # handle double extension (.cif.gz)
        name = Path(name).stem
    # split on dashes: ['AF', 'A0A023PYF4', 'F1', 'model_v6']
    # UniProt ID is always the second element
    return name.split('-')[1]


def parse_tar(tar_path, output_path):
    """Parse all CIF files inside an AlphaFold DB tar archive.
    Reads directly from tar — no need to extract 6000+ files to disk.
    The yeast tar has 12,110 files (6,055 CIF + 6,055 PDB)."""

    n_proteins = 0
    with tarfile.open(tar_path, 'r') as tar, open(output_path, 'w') as out:
        out.write("uniprot_id\tresidue_number\tplddt\n")

        for member in tar:
            # The tar contains both .cif.gz and .pdb.gz — skip everything
            # that isn't a CIF file
            if not (member.name.endswith('.cif') or member.name.endswith('.cif.gz')):
                continue

            # Extract file into memory (not to disk)
            fileobj = tar.extractfile(member)
            if fileobj is None:
                continue

            # tarfile gives us raw bytes; BioPython needs text.
            # If gzipped, gzip.open handles both decompression and text decoding.
            # If not gzipped, TextIOWrapper converts binary to text.
            if member.name.endswith('.gz'):
                fileobj = gzip.open(fileobj, 'rt', encoding='utf-8')
            else:
                fileobj = io.TextIOWrapper(fileobj, encoding='utf-8')

            uniprot_id = extract_uniprot_from_filename(member.name)

            # If one CIF is corrupted, skip it instead of crashing the whole run
            try:
                residues = parse_single_cif(fileobj)
            except Exception as e:
                print(f"skip {member.name}: {e}", file=sys.stderr)
                continue

            # Write one TSV row per residue
            for res_num, plddt in residues:
                out.write(f"{uniprot_id}\t{res_num}\t{plddt:.2f}\n")

            n_proteins += 1

    print(f"{n_proteins} proteins -> {output_path}", file=sys.stderr)


def parse_single_file(cif_path, output_path):
    """Parse one CIF file from disk — for testing or one-off checks."""

    uniprot_id = extract_uniprot_from_filename(cif_path)
    residues = parse_single_cif(cif_path)

    with open(output_path, 'w') as out:
        out.write("uniprot_id\tresidue_number\tplddt\n")
        for res_num, plddt in residues:
            out.write(f"{uniprot_id}\t{res_num}\t{plddt:.2f}\n")

    print(f"1 protein, {len(residues)} residues -> {output_path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI: accepts a .tar (whole proteome) or a single .cif file
# Output is always a 3-column TSV: uniprot_id, residue_number, plddt
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: 04_parse_cif.py <input.tar|input.cif> <output.tsv>",
              file=sys.stderr)
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    if input_path.endswith('.tar'):
        parse_tar(input_path, output_path)
    elif input_path.endswith('.cif'):
        parse_single_file(input_path, output_path)
    else:
        print(f"Error: expected .tar or .cif, got {input_path}", file=sys.stderr)
        sys.exit(1)
