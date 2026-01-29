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

# rf holds the codon information for the three reading frames
rf = [Codon(), Codon(), Codon()]

seq_n = 0
while cds.next():
    # sum codon usage over all CDS sequences for all three reading frames
    seq_n += 1
    print(f'{seq_n}\t{cds.id}')
    for frame in range(3):
        # codon counts for all three frames, includes stop codons
        rf[frame].add_from_dna(cds.seq, frame)

    if seq_n > 10:
        break

for frame in range(3):
    # calculate codon frequencies from counts, all three reading frames
    print(f'frame: {frame}\tcodons: {rf[frame].n}')
    rf[frame].update_frequencies()

familycount = [Codon(), Codon(), Codon()]
for frame in range(3):
    # calculate frequenccy of codon families in all three frames
    familycount[frame].add_from_codon(rf[frame])

preference = []
for frame in range(3):
    # codon preference for all three frames
    preference.append(rf[frame] / familycount[frame])

exit(0)
