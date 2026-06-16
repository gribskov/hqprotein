"""=====================================================================================================================
select_isoforms.py

Construct a table of multiple predicted isofroms based on Orthofinder WorkingDirectory/SequenceIDs
Format:
gene_gid: transcript_id
0_0: g1.t1
0_1: g2.t1
...
0_13: g14.t1
0_14: g14.t2
0_15: g15.t1

row number gives the ID in sequential numbering (used in some orthofinder outputs)
gene_id: species_gid - sequential numbering withing species
transcript.id: gene.transcript given in original GFF input (in this case braker3)

Usage:
select_isoforms.py SequenceIDs.txt isoform.list.txt

2026-06-16 gribskov
====================================================================================================================="""
import sys

# ======================================================================================================================
# Main
# ======================================================================================================================
if __name__ == '__main__':
    exit(0)
