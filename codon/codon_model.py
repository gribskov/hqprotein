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

aa2codon = defaultdict(list)
aa_count = {}
for codon in Codon.codon2aa:
    aa = Codon.codon2aa[codon]
    aa2codon[aa].append(codon)

# aa_count = { aa:[0,0,0] for aa in aa2codon}
# for aa in aa_count:
#     for codon in aa2codon[aa]:
#         for frame in range(3):
#             aa_count[aa] +=


for aa in sorted(aa2codon):
    print(f'{aa}\t{aa2codon[aa]}')
    for codon in sorted(aa2codon[aa]):
        print(f'\t{codon}\t{rf[0].frequency[codon]:.4f}')



exit(0)
