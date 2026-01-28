"""=====================================================================================================================
# codon_model.py
#
# Michael Gribskov 1/28/2026
====================================================================================================================="""
from include.fasta import Fasta
from codon import Codon
from collections import defaultdict

dnafile = 'data/z.tritici.IP0323.reannot.cds.fasta'
cds = Fasta(filename=dnafile)

# ref holds the codon information for the three reading frames
rf = [Codon(), Codon(), Codon()]

seq_n = 0
while cds.next():
    seq_n += 1
    print(f'{seq_n}\t{cds.id}')
    for frame in range(3):
        rf[frame].add_from_dna(cds.seq, frame)

    if seq_n > 10:
        break

for frame in range(3):
    print(f'frame: {frame}\tcodons: {rf[frame].n}')
    rf[frame].update_frequencies()

aacount = []
frame = 0
aacount[frame] = Codon()
aacount[frame].add_from_codon(rf[0])



# for aa in sorted(aa2codon):
#     print(f'{aa}\t{aa2codon[aa]}')
#     for codon in sorted(aa2codon[aa]):
#         print(f'\t{codon}\t{rf[0].frequency[codon]:.4f}')

exit(0)
