"""=====================================================================================================================
# codon.py
#
# Michael Gribskov 1/28/2026
====================================================================================================================="""


class Codon:
    """-----------------------------------------------------------------------------------------------------------------
    Codon

    codon2aa - class variable with standard codon definitions stored as dict
    aa2codon - class variable with list of codons for each one-letter amino acid residue, aa2codon can be populated
               from codon2aa so only codon2aa needs to be defined.
    -----------------------------------------------------------------------------------------------------------------"""
    # class variables
    # noinspection DuplicatedCode
    # Carefully proofread, earlier version had TAT as a T codon instead of Y
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

                "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y",
                "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
                "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
                "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}
    aa2codon = {}

    def __init__(self):
        """-------------------------------------------------------------------------------------------------------------
        Constructor for Codon class.
        Instance variables:
            count: integer count of each defined symbol (codon or AA, for instance)
            n: sum of counts for all symbols
        -------------------------------------------------------------------------------------------------------------"""
        self.count = {c: 0.0 for c in Codon.codon2aa}
        self.frequency = {c: 0.0 for c in Codon.codon2aa}
        self.n = 0

        # construct aa2codon from codon2aa
        if not Codon.aa2codon:
            # only construct aa2codon for first instance
            for codon in Codon.codon2aa:
                aa = Codon.codon2aa[codon]
                if aa not in Codon.aa2codon:
                    Codon.aa2codon[aa] = []
                Codon.aa2codon[aa].append(codon)

    def __truediv__(self, denom):
        """-------------------------------------------------------------------------------------------------------------
        overload / operator.
        currently supports division by another instance of Codon, int, or float

        :param denom: various       denominator for division, Codon, int, and float supported
        :return: Codon              new codon with result of division
        -------------------------------------------------------------------------------------------------------------"""
        if isinstance(denom, Codon):
            result = Codon()
            for codon in self.count:
                result.count[codon] = self.count[codon] / denom.count[codon]
                result.n = sum(result.count.values())
            return result
        elif isinstance(denom, (int, float)):
            result = Codon()
            for codon in self.count:
                result.count[codon] = self.count[codon] / denom
                result.n = sum(result.count.values())
            return result

        # only reach here for unknown denominator type
        return NotImplemented

    def __str__(self):
        """-------------------------------------------------------------------------------------------------------------
        string representation. Prints in the traditional format for the codon table.
        :return: str
        -------------------------------------------------------------------------------------------------------------"""
        out = f'{self.n:.3f} observations\n'
        for b0 in 'TCAG':
            for b1 in 'TCAG':
                for b2 in 'TCAG':
                    codon = f'{b0}{b2}{b1}'
                    out += f'   {codon} {self.codon2aa[codon]} {self.count[codon]:6.3f} '
                out += '\n'
            out += '\n'

        return out.rstrip()


    def add_from_dna(self, dna, frame=0):
        """-------------------------------------------------------------------------------------------------------------
        Add the codon counts from a DNA sequence in fasta format to the current count

        :param dna: Fasta       fasta formatted DNA sequence
        :param frame: int       reading frame to use, 0 indicates start with base 0
        :return: int            number of codons added
        -------------------------------------------------------------------------------------------------------------"""
        dna = dna.upper()
        for start in range(frame, len(dna) - 2, 3):
            self.n += 1
            self.count[dna[start:start + 3]] += 1

        return self.n

    def add_from_codon(self, codon):
        """-------------------------------------------------------------------------------------------------------------
        Build amino acid (codon family) count from an existing Codon object. store the family count in each codon in
        the family. This makes it simple to divide codon counts by family counts to get codon preference

        :param codon: Codon     Codon object with codon counts
        :return: int            total counts (should equal number of input codons)
        -------------------------------------------------------------------------------------------------------------"""
        n = 0
        for thiscodon in codon.count:
            aa = self.codon2aa[thiscodon]
            for member in Codon.aa2codon[aa]:
                self.count[member] += codon.count[thiscodon]

            n += codon.count[thiscodon]

        self.n = n
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
    print(f'\n{"*"*80}\nRead codons of coding frame from DNA\n{"*"*80}')
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
    print(f'Sequence is {len(seq)} bases long and ends in {lastcodon}')

    codon_n = coding.add_from_dna(seq)
    print(f'expect {len(seq) // 3} codons, and no stop codons except {lastcodon}.')
    print(f'\n{coding}')

    #-----------------------------------------------------------------------------------------------
    print(f'\n{"*"*80}\nTest update frequencies; all values should be 1/64 = 0.016\n{"*"*80}')
    t1 = Codon()
    for codon in t1.codon2aa:
        t1.count[codon] = 1
        t1.n += 1
    t1.update_frequencies()
    for codon in t1.codon2aa:
        t1.count[codon] = t1.frequency[codon]
    print(t1)

    #-----------------------------------------------------------------------------------------------
    print(f'\n{"*"*80}\nTest conversion to family counts\n{"*"*80}')
    t2 = Codon()
    for codon in t2.codon2aa:
        t2.count[codon] = 1
        t2.n += 1
    t3 = Codon()
    t3.add_from_codon(t2)
    print(t3)

    #-----------------------------------------------------------------------------------------------
    # TODO test division by codon


    #-----------------------------------------------------------------------------------------------
    print(f'\n{"*"*80}\nTest division by 2; should match original count\n{"*"*80}')
    half = coding / 2
    print(f'{half}')


    exit(0)
