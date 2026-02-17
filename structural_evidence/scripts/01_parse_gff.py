#!/usr/bin/env python3
#
# parse_gff.py - GFF3 parser for structural evidence pipeline
# Alex Baena, Feb 2026
#
# Not needed yet. For now we have pre-extracted proteins (FASTA) and
# pre-computed AlphaFold structures (CIF), so we can skip GFF parsing
# and just read those files directly.
#
# This becomes necessary when:
#   - A user provides a new GFF + genome with no pre-computed structures
#     (tool needs to extract CDS, translate, and run AlphaFold)
#   - Segment-level analysis requires exon coordinates from the GFF
#   - Corruption experiments need to modify CDS boundaries
#
# Until then, any GFF data we need (e.g. gene-to-UniProt mapping) can
# be grabbed inline with a few lines of Python or grep.
#
