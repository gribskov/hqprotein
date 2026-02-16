#!/usr/bin/env python3
#
# 02_extract_proteins.py - Extract CDS from genome and translate to protein
# Alex Baena, Feb 2026
#
# Not needed yet. S. cerevisiae already has pre-translated proteins
# (orf_trans_all_R64-5-1_20240529.fasta.gz from SGD).
#
# This becomes necessary when a user provides a new genome + GFF
# with no pre-extracted proteins. Takes CDS coordinates from 01,
# pulls DNA from the genome FASTA, and translates to protein sequences.
#
