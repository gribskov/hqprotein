"""=====================================================================================================================
# codon_model.py
#
# Michael Gribskov 1/28/2026
====================================================================================================================="""
from include.fasta import Fasta
from codon import Codon
from collections import defaultdict

def calculate_lr( rf ):
    """-----------------------------------------------------------------------------------------------------------------
    frequencies in rf are the P(codon|frame=coding), if the priors for the three reading
    frames are equal, they factor out and
    P(coding|codon) = P(codon|coding) / (P(codon|coding) + P(codon|f1) + P(codon|f2))

    :param rf: list of Codon        counts and frequencies of codons in the three reading frames
    :return: Codon                  likelihood ration for coding vs noncoding
    -----------------------------------------------------------------------------------------------------------------"""
    result = Codon()
    a = Codon.codon2aa
    coding = rf[0] / rf[0].n
    psum = coding + 0
    psum += rf[1] / rf[1].n
    psum += rf[2] / rf[2].n

    result = coding / psum
    psum.update_frequencies()

    return result



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

    if seq_n > 100:
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

lr = calculate_lr(rf)

exit(0)
