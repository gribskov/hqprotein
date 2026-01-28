"""=====================================================================================================================
# codon.py
#
# Michael Gribskov 1/28/2026
====================================================================================================================="""


class Codon:
    """-----------------------------------------------------------------------------------------------------------------
    Codon

    codon2aa - class variable with standard codon definitions stored as dict
    aa2codon - class variable with list of codons for each one-letter amino acid residue
    -----------------------------------------------------------------------------------------------------------------"""
    codon2aa = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
                "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
                "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
                "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",

                "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
                "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
                "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
                "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",

                "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
                "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
                "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
                "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",

                "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "T",
                "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
                "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
                "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}
    aa2codon = {}

    def __init__(self):
        """-------------------------------------------------------------------------------------------------------------
        Constructor for Codon class
        -------------------------------------------------------------------------------------------------------------"""
        self.count = {c: 0 for c in Codon.codon2aa}
        self.frequency = {c: 0 for c in Codon.codon2aa}
        self.n = 0

        # construct aa2codon from codon2aa
        for codon in Codon.codon2aa:
            aa = Codon.codon2aa[codon]
            if aa not in Codon.aa2codon:
                Codon.aa2codon[aa] = []
            Codon.aa2codon[aa].append(codon)

    def add_from_dna(self, dna, frame=0):
        """-------------------------------------------------------------------------------------------------------------
        Add the codon counts from a DNA sequence in fasta format to the current count

        :param dna: Fasta       fasta formatted DNA sequence
        :param frame: int       reading frame to use, 0 indicates start with base 0
        :return: int            number of codons added
        -------------------------------------------------------------------------------------------------------------"""
        start = frame
        dna = dna.upper()
        for start in range(frame, len(dna)-2, 3):
            self.n += 1
            self.count[dna[start:start + 3]] += 1

        return self.n

    def add_from_codon(self, codon):
        """-------------------------------------------------------------------------------------------------------------
        build amino acid (codon family) count from an existing Codon object

        :param codon: Codon     Codon object with codon counts
        :return: int            total counts (should equal number of input codons)
        -------------------------------------------------------------------------------------------------------------"""
        self.count = {aa:0 for aa in Codon.aa2codon}
        n = 0
        for codon in codon.count:
            aa = self.codon2aa[codon]
            self.count[aa] += codon.count[codon]
            n + codon.count[codon]

        return n


    def update_frequencies(self):
        """-------------------------------------------------------------------------------------------------------------
        Using the current count and n, calculate the frequency of each codon, P(codon|dna_seq)
        frequency is stored in self.frequency
        Use the returned sum of the frequencies to check for errors (such as incorrect number  of codons, self.b)

        :return: float      sum of frequency, should be 1.0
        -------------------------------------------------------------------------------------------------------------"""
        count = self.count
        n = self.n
        self.frequency = {c: count[c] / n for c in count}

        return sum(self.frequency.values())


# ######################################################################################################################
# Testing
# ######################################################################################################################
if __name__ == '__main__':
    print(f'Read codons of coding frame from DNA')
    coding = Codon()

    seq = '''ATGAGGTTCCACGTTCATTCGACGCCATTCTACCAACGCATAGCCTGCAACACCACATCGACCATCACTG
             CGCACTCCACAAGTAGCCGGCAGACCATACCGTCAAGTTCGAGAATGTCCGTCCAACAAACAATCAACGA
             GGAGAGGATGGAGAACCGGGCAACATTGTCGCTCCAGTGCGAGTCCAGGCTGATCACCTCAGGCCACACG
             ACGCCTGACTTCATCAATCGGATGAACAGCCTCTTACAAGACTACGTCGACATCGAGTACGGTCTCATCC
             TGCGCATTCGCCAAATTCTGAAGAACCGGGGTCTCGACCGGCAATGTTGGACCCCCGACCGCATCGACGG
             GTACGAGAGATTCACGGAGTCATGGATCTACAGGCTGGGCTATTCAACGAACTTTGAAAGGGACTTCTCC
             CCATCCGACATGATGGCTTTGATCGCAGATCTCGACGCTCTCCTGCTAGACGCTCAAATTGCGAATGATA
             TCGAGTTCGAATGGCAGGATGTCTCGGCTGAGGGCGACCGGTTAGGGAAATGCTGCGAAGGTTCGCATCC
             GACCAAGCGGAAGATCATCATCTTTGCACTGGGTCGTACCTCTGCGGACATCTTAATCGGCACCCTTGTG
             CATGAGATGTGTCATGCCTTCATCGATATATCGCTCATTGATGAGCTCGGCCATGAATCAGCCTACGACC
             CTGAGGGTGGCAACACCGAGCCGTGGACCACAGCAATCGGAAGCACAGGCCATGGCTACTGCTGGCAAAT
             GCTCGCTATGATCATGAACGACTGTTCGGAAGCCATTGGGTACCCTTTTGACCTTGGTTTTGAAGACGCT
             GGAAACAGCGACTTCGACGGACACCTCGGCGATCGCCTCTCCGGTCTAGGGCGATACGAAGTGGAGGAAC
             CCAGCTCGCCATGGGCTACTTTACACGGTCGGCAGCAGCGGCGGGCCGAGGTTGACAGGGCTGAGAGGGA
             GGAAGCCAGGCAGGGGGAAGAGGAATCGTCCGAAGGATCTTCACCCGGAGAAGATACGATGATGGAGCTC
             ACACAGAGCGATATCGTGGCGACTTGA'''
    seq = seq.replace('\n', '').replace(' ', '')
    lastcodon = seq[-3:len(seq)]
    print(f'\tSequence is {len(seq)} bases long and ends in {lastcodon}')

    codon_n = coding.add_from_dna(seq)
    print(f'\texpect {len(seq) // 3} codons. codons read: {codon_n}')

    codon_n += coding.add_from_dna(seq)
    print(f'\tRead second time, expect {len(seq) // 3 * 2} codons. ', end='')
    print(f'codons read: {codon_n}\tcount({lastcodon}): {coding.count[lastcodon]}')

    exit(0)
