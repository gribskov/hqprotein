"""=====================================================================================================================
splice_model.py

frequency model of splice donor and acceptor based on genome and annotation (GFF3)
frequencies will be P(base(pos)|acceptor) and P(base(pos)|donor)

Synopsis
# Create empty object
splice = SpliceSite()

# See sample code for reading splice junctions from gff below

# read saved positional probabilities from file
fq_from_file = SpliceSite()

# convert counts to frequency and print in specific format
frequency = splice.frequency()
frequency.fieldwidth = 6
frequency.precision = 3
print(f'\n{frequency}')

# Positional Shannon entropy
h = frequency.information()

# sharpen positional probabilities
fq_from_file = SpliceSite()
n = fq_from_file.read('frequency.dat')
print(f'\n{n} frequencies read from frequency.dat')
fq_from_file.fieldwidth = 6
fq_from_file.precision = 3

print(f'Sharpened frequencies')
fq_from_file.sharpen(2)
print(fq_from_file)

Michael Gribskov 2/6/2026
====================================================================================================================="""
from sys import stderr
from math import log
from include.gff.gff2 import GxfSet
from include.sequence.fasta import Fasta
from collections import defaultdict
import numpy as np
from pssm import PSSM


class SpliceSite:
    """=================================================================================================================
    Frequency matrices for splice donor and acceptor
    ================================================================================================================="""
    base_complement = str.maketrans('ACGT', 'TGCA')
    dna = 'ACGT'
    protein = 'ACDEFGHIKLMNPQRSTVWY'

    def __init__(self, pre=5, post=10, seqtype='dna'):
        """-------------------------------------------------------------------------------------------------------------
        holds position specific counts or frequencies for splice donor and acceptor sites
        pre and post are defined with respect to the exon side with before indicating the exon side

               donor           acceptor
              pre|post         pre|post
        exon-----|-----intron-----|-----exon

        for pre = 2 and post = 5m with e and i indicating exon positions
        donor is eeiiiii and acceptor is iiiiiee

        donor       position specific counts
        acceptor    position specific counts
        pre         bases before site (exon side)
        post        bases after site (intron side)
        -------------------------------------------------------------------------------------------------------------"""
        self.donor = PSSM(title='donor, offset=pre')
        self.acceptor = PSSM(title='acceptor', offset=post)
        self.pre = pre
        self.post = post
        if seqtype == 'dna':
            # convert sequence characters to integer index
            alphabet = getattr(self, seqtype)
            a2i = {alphabet[i]: i for i in range(len(alphabet))}
        self.a2i = a2i

    def add_junction(self, strand, donorpos, acceptorpos, sequence):
        """-------------------------------------------------------------------------------------------------------------
        add a donor/acceptor pair from sequence based on coordinates (from GFF). Input sequences are forced to
        uppercase. Currently, no checking for non A,C,G,T bases.

        :param strand: str          '+' or '-', strand of exons
        :param donorpos: int        end of donor exon
        :param acceptorpos: str     strand, '+' or '-'
        :param sequence: str        DNA sequence
        :return: int, int           number of donor/acceptor sequences
        -------------------------------------------------------------------------------------------------------------"""
        sequence = sequence.upper()
        # self.donor_n += 1
        # self.acceptor_n += 1
        left = {}
        right = {}
        if strand == '+':
            # donor/acceptor site endpoints for + strand
            left['donor'] = donorpos - self.pre
            right['donor'] = donorpos + self.post
            left['acceptor'] = acceptorpos - self.post - 1
            right['acceptor'] = acceptorpos + self.pre - 1
            # print(f'{strand}\t{sequence[left['donor']:right['donor']]}\t{sequence[left['acceptor']:right['acceptor']]}')
        else:
            # donor/acceptor site endpoints for - strand
            left['donor'] = donorpos - self.post - 1
            right['donor'] = donorpos + self.pre - 1
            left['acceptor'] = acceptorpos - self.pre
            right['acceptor'] = acceptorpos + self.post
            # d = sequence[left['donor']:right['donor']].translate(SpliceSite.base_complement)[::-1]
            # a = sequence[left['acceptor']:right['acceptor']].translate(SpliceSite.base_complement)[::-1]
            # print(f'{strand}\t{d}\t{a}')

        sitepos = 0
        for jtype in ('donor', 'acceptor'):
            site = sequence[left[jtype]:right[jtype]]
            pssm = getattr(self, jtype)
            if strand != '+':
                # reverse complement - strand
                site = site.translate(SpliceSite.base_complement)[::-1]

            sarr = np.array(list(site))
            count = np.zeros((4, len(site)), dtype=int)
            a2i = self.a2i
            for base in a2i:
                mask = sarr == base
                count[a2i[base], mask] += 1

            pssm.matrix += count

        return self.donor.n, self.acceptor.n

    # @staticmethod
    # def complement(sequence):
    #     """-------------------------------------------------------------------------------------------------------------
    #     complements but does not reverse the sequence string A->T, C->G, G->C, T->A. No checking for ambiguous or
    #     incorrect bases
    #
    #     :param sequence: str    DNA sequence, 5' to 3'
    #     :return: str            complementary strand in 3' to 5' order
    #     -------------------------------------------------------------------------------------------------------------"""


# ======================================================================================================================
# End of class SpliceSite
# ======================================================================================================================

# ======================================================================================================================
# main
# ======================================================================================================================
if __name__ == '__main__':
    gff_file = 'data/z.tritici.IP0323.reannot.gff3'
    gff = GxfSet(file=gff_file, fmt='gff')
    genome = Fasta(filename='data/z.tritici.IP0323.fasta')

    # read exon features from GFF file. Genes with more than one exon have a splice junction and are stored
    # as a list of exons in junction
    feature_n = gff.feature_get(['exon'])
    junction = defaultdict(list)
    feature_count = 0
    current = ''
    exon_set = []
    exon_set_n = 0
    for feature in gff.features:
        # find exons from multi exon genes and store in junction
        feature_count += 1
        # print(f'\t\t{feature.start}\t{feature.end}\t{feature.strand}\t{feature.attribute['Parent']}')
        if current == feature.attribute['Parent']:
            exon_set.append(feature)
        else:
            if len(exon_set) > 1:
                # multiple exons, get acceptor and donor sites
                print(f'\tmultiple exon gene: {feature.seqid}\t{feature.attribute['Parent']}')
                sequence = exon_set[0].seqid
                junction[sequence].append(exon_set)
            exon_set = [feature]
            exon_set_n += 1

        current = feature.attribute['Parent']

        # TODO remove limit input for debug
        if feature_count > 10: break

    print(f'\nexons: {feature_count}')
    print(f'multiple exon sets: {exon_set_n}')

    # read one sequence at a time from genome and process multiple exon genes in that sequence
    splice = SpliceSite()
    jtotal = 0
    for sequence in genome:
        jcount = 0
        print(f'\nProcessing {sequence.id}')
        # id = sequence.id
        for exon_set in junction[sequence.id]:
            parent = exon_set[0].attribute['Parent']
            strand = exon_set[0].strand

            # print(f'\tgene:{parent} strand: {strand}')
            for i in range(len(exon_set) - 1):
                if strand == '+':
                    donor = exon_set[i].end
                    acceptor = exon_set[i + 1].start
                else:
                    acceptor = exon_set[i].end
                    donor = exon_set[i + 1].start

                splice.add_junction(strand, donor, acceptor, sequence.seq)
                jcount += 1
        print(f'\t{jcount} junctions added from {sequence.id}')
        jtotal += jcount

        # TODO remove limit input for debug
        break

        splice.fieldwidth = 7
        print(f'\n{splice}')

    print(f'{jtotal} splice junctions analyzed')

    frequency = splice.frequency()
    frequency.fieldwidth = 6
    frequency.precision = 3
    print(f'\n{frequency}')

    h = frequency.information()
    print(f'\nPositional Shannon entropy:')
    for jtype in ('donor', 'acceptor'):
        linepos = frequency.pre
        if jtype == 'acceptor':
            linepos = frequency.post

        print(f'\n{jtype:<8s}', end='\t')
        pos = 0
        for val in h[jtype]:
            if pos == linepos:
                print(f'{'|':>6s}', end='')
            print(f'{2 + val:6.3f}', end='')

            pos += 1

    print()

    with open('frequency.dat', 'w') as fq:
        print(f'{frequency}', file=fq)

    fq_from_file = SpliceSite()
    n = fq_from_file.read('frequency.dat')
    print(f'\n{n} frequencies read from frequency.dat')
    fq_from_file.fieldwidth = 6
    fq_from_file.precision = 3

    print(f'Sharpened frequencies')
    fq_from_file.sharpen(2)
    print(fq_from_file)

    exit(0)
