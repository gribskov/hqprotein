"""=====================================================================================================================
codon.py implements the Codon class which is used to manipulate codon usage/frequency tables

# create a codon usage table from DNA sequence, frame 0
# this gets the raw counts. frame is optional, defaults to zero (start at first base in sequence)
codon_usage = Codon()
codon_usage.add_from_dna(dnaseq, frame=0)

# codon usage as fraction of total codon, i.e., P(codon|frame). The frequencies are stored in the count attribute
codon_frequency = codon_usage / codon_usage.n

# add plus 1 prior to codon counts
codon_usage += 1

# calculate codon preference (frequencies of codons within each family).. The frequencies are stored in the count
# attribute
family_count = codon_usage.add_from_codon()
preference = codon_usage / family_count

# make a table with a constant value = 1
prior = Codon()
prior += 1

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

    def __add__(self, addend):
        """-------------------------------------------------------------------------------------------------------------
        Add two codon objects
        overload + operator. codon + 1 or codon1 + codon2
        currently supports addition by another instance of Codon, int, or float
        returns a new Codon object

        :param addend: various       addend for addition, Codon, int, and float supported
        :return: Codon              new codon with result of addition
        -------------------------------------------------------------------------------------------------------------"""
        if isinstance(addend, Codon):
            result = Codon()
            for codon in self.count:
                result.count[codon] = self.count[codon] + addend.count[codon]
            result.n = sum(result.count.values())
            return result
        elif isinstance(addend, (int, float)):
            result = Codon()
            for codon in self.count:
                result.count[codon] = self.count[codon] + addend
            result.n = sum(result.count.values())
            return result

        # only reach here for unknown addend type
        return NotImplemented

    def __radd__(self, addend):
        """-------------------------------------------------------------------------------------------------------------
        Add two codon objects
        overload + operator with the Codon object on the right side, i.e., 1 + codon
        currently supports addition by another instance of Codon, int, or float
        returns a new Codon object

        :param addend: various       addend for addition, Codon, int, and float supported
        :return: Codon              new codon with result of addition
        -------------------------------------------------------------------------------------------------------------"""
        if isinstance(addend, Codon):
            result = Codon()
            for codon in self.count:
                result.count[codon] = self.count[codon] + addend.count[codon]
            result.n = sum(result.count.values())
            return result
        elif isinstance(addend, (int, float)):
            result = Codon()
            for codon in self.count:
                result.count[codon] = self.count[codon] + addend
            result.n = sum(result.count.values())
            return result

        # only reach here for unknown addend type
        return NotImplemented

    def __truediv__(self, denom):
        """-------------------------------------------------------------------------------------------------------------
        divide one codon table by another, or by a constant
        overload / operator.
        currently supports division by another instance of Codon, int, or float
        returns a new Codon object

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

    def family(self):
        """-------------------------------------------------------------------------------------------------------------
        Build amino acid (codon family) count from an existing Codon object. store the family count in each codon in
        the family. This makes it simple to divide codon counts by family counts to get codon preference.

        :return: Codon          new Codon with counts of family stored in each docon
        -------------------------------------------------------------------------------------------------------------"""
        n = 0
        result = Codon()
        for thiscodon in Codon.codon2aa:
            aa = self.codon2aa[thiscodon]
            for member in Codon.aa2codon[aa]:
                # add counts for each member of the codon family to each member
                result.count[thiscodon] += self.count[member]

            n += result.count[thiscodon]

        result.n = n
        return result


# ######################################################################################################################
# Testing
# ######################################################################################################################
if __name__ == '__main__':
    print(f'\n{"*" * 80}\nRead codons of coding frame from DNA\n{"*" * 80}')
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

    # -----------------------------------------------------------------------------------------------
    print(f'\n{"*" * 80}\nTest conversion to frequencies\n{"*" * 80}')
    print('Expect all values to be 1/64 = 0.016')
    # test data is table with one count for each codon
    t1 = Codon()
    t1 += 1

    t2 = t1 / t1.n
    print(t2)

    # -----------------------------------------------------------------------------------------------
    print(f'\n{"*" * 80}\nTest conversion to family counts\n{"*" * 80}')
    print(f'Expect counts to be the number of codons in the corresponding synonymous codon family. Total count = 244')
    t1 = Codon()
    t1 += 1
    t2 = t1.family()
    print(t2)
    # print family sums
    print('\nCumulative counts by codon family (family, family size, count)')
    total_count = 0
    for aa in Codon.aa2codon:
        family = Codon.aa2codon[aa]
        for codon in family:
            total_count += t2.count[codon]
        print(f'\t{aa}\t{len(family)}\t{total_count}')

    # -----------------------------------------------------------------------------------------------
    print(f'\n{"*" * 80}\nTest division of one codon table by another\n{"*" * 80}')
    print(' All values should be 1/synonymous_codon_family_size.')
    t1 = Codon()
    t1 += 1

    t3 = t1 / t2
    print(t3)

    # -----------------------------------------------------------------------------------------------
    print(f'\n{"*" * 80}\nTest addition\n{"*" * 80}')

    t1 = Codon()
    t1 += 1
    print('t1 += 1; All values should be 1, observations=64')
    print(t1)

    t2 = 2 + t1
    print('\nt2 = 2 + t1; All values should be 3, observations=192')
    print(t2)

    t3 = t1 + t2
    print('\nt3 = t1 + t2; All values should be 4, observations = 256')
    print(t3)

    # -----------------------------------------------------------------------------------------------
    print(f'\n{"*" * 80}\nTest division by 2; should match original count\n{"*" * 80}')
    print(f'Expect values to be 2, observations = 128.')
    half = t3 / 2
    print(f'{half}')

    exit(0)
