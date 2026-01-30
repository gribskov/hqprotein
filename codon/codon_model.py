"""=====================================================================================================================
# codon_model.py
#
# Michael Gribskov 1/28/2026
====================================================================================================================="""
from include.fasta import Fasta
from codon import Codon
import matplotlib.pyplot as plt
import numpy as np


def calculate_lr(rf):
    """-----------------------------------------------------------------------------------------------------------------
    frequencies in rf are the P(codon|frame=coding), if the priors for the three reading
    frames are equal, they factor out and
    P(coding|codon) = P(codon|coding) / (P(codon|coding) + P(codon|f1) + P(codon|f2))

    :param rf: list of Codon        counts and frequencies of codons in the three reading frames
    :return: Codon                  likelihood ratio for coding vs noncoding
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


def window_average(dna, lr, window=5):
    """-----------------------------------------------------------------------------------------------------------------
    sliding window average of the codon score (lr)

    :param dna: string          DNA sequence
    :param lr: Codon            statistic to average (lr.count)
    :param window: integer      window size in codons
    :return: list               elements are of [pos,ave_lr.count]
    -----------------------------------------------------------------------------------------------------------------"""
    result = []
    stat = 0.0
    for pos in range(0, window * 3, 3):
        # fill first window
        codon = dna[pos:pos + 3]
        stat += lr.count[codon]

    # fill successive windows to end
    half = (window - 1) / 2
    for pos in range(window * 3, len(dna) - window * 3, 3):
        start = pos - window * 3
        result.append([start/3 + half, stat / window])
        stat -= lr.count[dna[start:start + 3]]
        stat += lr.count[dna[pos:pos + 3]]
        start += 3

    result.append([start/3+half, stat / window])

    return result


def plot(xydata, window, line):
    """-----------------------------------------------------------------------------------------------------------------
    plot with matplotlib

    :param xydata: list     [ [x0,y0], [x1,y1], ... ]
    :return: None
    -----------------------------------------------------------------------------------------------------------------"""
    # Data for plotting
    x = [a[0] for a in xydata]
    y = [a[1] for a in xydata]
    fig, ax = plt.subplots()
    ax.plot(x, y)

    ax.set(xlabel='codon', ylabel='ave lr stat',
           title=f'codon usage lr, window = {window}')
    ax.grid()
    plt.axvline(x=line, color='red', label='axvline - full height')

    fig.savefig("test.png")
    plt.show()

    return None


# DNA coding sequence for testing
dna = """GGCTTTGCCTCCGCCGCCAGCTTGTTGCGGAACCAGGCATTCGAGCTTTCAGCTCATAGGACGGAGGAGG
GGGTCGAAGATATCGAGCACAGGAAGAGCGTGAAGAACGGTCTTGCACGGGGGGAGGAGAGCGTCAAAGA
CAAGGTTGCCGGGATCGAGAAAGAACCTGTGAACAAGCGAGCGAGGAGCAGGAAAGCGCAAGAACAGACG
GACGAACATGACACACAAGCACCACAGCCTAGAGTACAGGCTGTCGATGATCTGTCGAGACCTTCCAAGA
CTGCTACCGGCAAGCTGAAGAAGACTGCAAAGTCGAGCAAGAGCGCAAAGTCTGGAGCAGAAAGTCTGTC
TCGCTCCCCTTTGAGCGAGTTCGCGCTGCCTGCTCAAACATCTGCGATTCGCAATGTGAGTCCGGCAGCA
GCCAGCGTTTCTGTCCTTGCCGGACGTCAAGAGACTCCTGACCTCGCAGCTTCCCGCAGAACGAGCATCT
CACTCTCCGAATTTGACTTCAACCCGAACAGCCCAAGCTCGAGACGAACCGTCTCGCCACTTGGCAAAAA
GGCGAGCGTGGCCAAGACTGGAAAGCATGGTTCCTCGCACCTCAACTCCGGGCGGACATCTGCCAGCAAC
AAAGTGTCAAAGTCTGCTCCTCGCAAATCCAAGACTACCGAATCACGACCAAAGTCCGGATCTCGGACGG
CGAGCATCTGTATTTCCGAGTACGCGTTGCAACCTGATAGCCCGCCGAAGTACGCTTCGCCGCCAGCGAT
GGAGACCTTGGATGACGTCCCAAAGCTTCGAGAGAATCAAGAGACCGACGTAGCACCACCCTTCATGGTG
CGATCGGCCACATCTTACGCAAGCGACTCCAACATACTCGGTCGATCACCCTCCCAATTTGCCTCGAGGT
ATCAGTATGGCGAGACACACTCAGCGCTGTCGACCTCTGAGCCGCAGCAGCCAGTCATCAGTCTGCCGGT
GACTACGGATGTTGCGACAGCGGCGTTCAAGGATAAGGCAATGGACACGGTATCTTTGAGTGAGCAGCCC
TCCCTCGCCAAGGCTCCTCGGAAACGCATGCCAAAAGCCAGGCAGTCAAAGGGTTCTGGTCCCGCGAGGT
CGAAGAAGCTTGCCAAGTCGGCCTCCATCATTCTCAATTCTGATGAGCCGGATCTGGACGCGGGCGACCA
CGCTGGCATTGCGCCGAATCCTGGTCCACCAGCGCCGGTCTTGCATGTGCTGACTGAACCAGCGGCGAAA
GCCAGGAAGGCTGTCAAGAAGCCGAGGGTCGCTCGTGGCAATAAGGAGGTTGAGAACGCCACAGCCTCAG
CCTCAGCTAATGCCATCGTGCCTCCGCTTCACGACATGCATGGCGGAGGCGTTCTAAAGGACATAGATGC
ACAGGGCGGATTCGACGCTCATCGCAACGCTTGGAGCCAATCTGCTTACTTCGTCGACAACGTCATCGCC
AATTGCCAGGCTCTACCACTGGCGGCGACAGAGGAATCTCTTCAGTCGGTGTTTGAAGTCGAAAGCATGA
GAAGGAACTCTACTCCACCTCCAGATCTCCCGCTGCCTCACTTTCGCCGACGAAGAAGCTGGACTCCAGC
ACGGAACACTGAATCGGCGACTGTCGACACGCCCCGTACGGTCGATTCGGGACGATTGGAAGAGTCTCCT
CCCTCTGGACCGAGTCTTGCTGCGTTGGTCAGTGGTTTCGGTTACAGTGCGCAGGCCACCAGTACCACTG
CCAACGCCCGCACTGGCACCGGTGAGCCTGCCGTGAAGAAACGGCGTATTGAGCTCGCGCATGAGACAGC
CGGCGGTGTTGTGGTTCGCAAAGCCAAAGATGCAGCTCCCGCGAAGCGAGTCAAGGCAGTCAAGAAAAAG
CCGCAAACAATTACCGATCTTGCCACCAAAGCATATCGCGCCGAACCCGAAATCGTGGCGGCGCAATCGA
CTGTGTCCGAATTCTTTGCTCCACAGAAGGAGTCGCTCACCTTGACCGAGGACGCGGCTGCGGGTGAGGC
ACTAACAGAGAAGCCAGCCAAAGTGAAAAAGACAAGAAAACGCAGGATCAAGGTTGGTGGAGATCCCGAG
GCTTTGGCAGCGATGGCGGAAAGGGCAGCAAAGCTGCAGAAGCCCGCGAAGGTCAAGAAAGTGAAGTTCA
ACGATGCGAATGCGGTACCGCCTCTGCTCTCTCCACAACGAGCTCGGGCGATCGAGACGCAGCAAGGATT
TCTGTTCGGCACCTCAAGTCAACTGGCGGTGGACGATTCGCCGACATTCATTCGGCAAATGCAAACGGCG
CTCCGAGAGTCTGAGTTGATGAGCACCCTGCAGGCACGTTTCAGCCCTCGTCGAAAAAGCTGTGCCAAAG
TTCCCAGCGCTCCTCATGGCACTTCCTTGTCAATTGGACAGGCCGCCCGGGACTTGTGGTGCACTGCCTC
GAGAGACCACGGCGGAGATGTTCTTGCGCTTACAGAACGGCCCAAGGCTCCCGTTGTTATCGAGCAGAAA
GCGCATAGCGAAGAAGACGCACATTCGCAATATGCCGCGGAGACATTGCCCGTTGTCGATGCTCTTCATG
CCATCGAGGAGACCACTGATGTCATTCCAGACGCTCAGCTTGACGGACAGGAAGTTGTGGACCTCTGTGC
CACATCGCCCGCGGGCGAGCGCGCTGTGCAAGACCCGGAACCGGCAGTGTTTGACCTTCCCAAGCCCGAA
GACTGCGATGAAAGTGCTCATGCGCTCCCCTCTGAACCGCACAATTCGGACGACAGCTGGATGCTATTGT
CAAGTGATGGCGCTGAGAAGGCGGAAGGGCATGCCTCCACCTTCCCTCTACCGATGATCCCAAACATCAT
GCGTGCAGCACCGGCCTTGGCCCGAGTAGCGACGTCACCAACGCGATTAAGAGCTGCCCTCCAGACCTTG
GACGGGAACGCAAGTGTCCTTGCTCAAAAATCTTCCGCAAGATTCTTGCAACAGCGAGCTTTTGCGAGCG
CAACCTCACCACCAGCGCAGAAGCGCCCTCGAGGGAGACCAAGAAAAGATGCATCATCTACCGAGAACAT
CCCTCCAGCCCAGAAAAAGAGAGGTCGCCCTCCGAAAAGCGCGACCACGTTCACACTGCCGGAAGTGAAG
GAAGGCAGAATCACGTTACCAGCCTCGCAGCAATCCACTTCCTCCGACTTCGTGAACATTGACGAAATTT
ACGACTCTGATCCTCCAACCCCCTCCCCGCCTCGCCGTCGGGCTGCCTCCACTTCCCCTTCGGTCAGTCC
ATTGGAGCTGCGGCTTTCCGGATCGCCTTCACCACGCGCCAAAGGGCCAGCCGTCGTGACAACAATGCTG
AAGTCCGGCGACCCGCAGTGGGCCTCCATCAGCGCCATCCTCTTCCCGCAGATCACGAAAGCTGTCCGCG
ATATACCTCCCAGCAACGACATGGTGAACCCGTCCTGGCACGAAAAGATCTTACTATTCGATCCGATTGT
GTTGGAAGACTTTACAGGTTGGCTCAACGGGCAGGGCTTGCGGGTCGAGACAAGGAAGGTCATTCCAAAG
GCGAAGAAGAAGGGTAAGAAGAAGAAGCAGGACGATGATACTCCGGAGGTGGCTGACGAGCTGGGGGAGT
ACGAGCTGTGCAGGGAAGAGATCAAGCCTTGGATGGCGCAGAAGTGGTGCGAGTCGAAGAGCATTTGCTG
TTTGTGGAAGGAGGGGATGAGGGGAGGGGTGAAGGTACAGTATTGA""".replace('\n', '')

# ===================================================================================================
# Main
# ===================================================================================================
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
print(f'sequence is {len(dna)} bases long')
# delete one base to create frameshift
half = len(dna) // 2
dna = dna[:half] + dna[half + 1:]
print(f'sequence is {len(dna)} bases long with frameshift at {half} (codon {half//3})')
ave = window_average(dna, lr,20)
plot(ave, 20, half/3)

exit(0)
