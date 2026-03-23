"""=====================================================================================================================
splice_model.py

frequency model of splice donor and acceptor based on genome and annotation (GFF3)
frequencies will be P(base(pos)|acceptor) and P(base(pos)|donor)

Synopsis
# Create empty SpliceSite object
splice = SpliceSite()

1. Read exons from GFF and select transcripts with at least two exons
2. Extract donor and acceptor sites from the reference genome
    a. Read chromosomes on at a time
    b. extract donor and acceptor sites and add to donor and acceptor PSSMs in splice
3. Convert to frequencies and write out

Michael Gribskov 2/6/2026
====================================================================================================================="""
from include.gff.gff2 import GxfSet
from include.sequence.fasta import Fasta
from collections import defaultdict
import numpy as np
from pssm import PSSM


class SpliceSite:
    """=================================================================================================================
    Frequency matrices for splice donor and acceptor. This class is mostly just a container for the donor and acceptor
    sites which are PSSM objects.

    Synopsis
    # Create empty object
    splice = SpliceSite()

    #add count for a splice junction from its sequence
    splice.add_junction(strand, donor, acceptor, sequence.seq)

    # print donor and acceptor
    PSSM.fieldwidth = len(maxval) + 2
    PSSM.precision = 0
    print(splice.donor)
    print(splice.acceptor)
    ================================================================================================================="""
    base_complement = str.maketrans('ACGT', 'TGCA')
    dna = 'ACGT'
    protein = 'ACDEFGHIKLMNPQRSTVWY'

    def __init__(self, pre=5, post=10):
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
        pre         bases before site (exon side for both donors and acceptors)
        post        bases after site (intron side or both donors and acceptors)
        -------------------------------------------------------------------------------------------------------------"""
        self.donor = PSSM(title='donor', rows='ACGT', offset=pre)
        self.acceptor = PSSM(title='acceptor', rows='ACGT', offset=post)
        self.pre = pre
        self.post = post

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

        for jtype in ('donor', 'acceptor'):
            site = sequence[left[jtype]:right[jtype]]
            pssm = getattr(self, jtype)
            if strand != '+':
                # reverse complement - strand
                site = site.translate(SpliceSite.base_complement)[::-1]

            pssm.add_counts_sequence(site)

        return self.donor.n, self.acceptor.n


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

    ####################################################################################################################
    # read exon features from GFF file. Genes with more than one exon have a splice junction and are stored
    # as a list of exons in junction
    ####################################################################################################################
    feature_n = gff.feature_get(['exon'])
    junction = defaultdict(list)
    feature_count = 0
    current = ''
    exon_set = []
    exon_set_n = 0
    for feature in gff.features:
        # find exons from multi exon genes and store in junction
        feature_count += 1
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

        # limit input for debug
        # if feature_count > 500: break

    print(f'\nexons: {feature_count}')
    print(f'multiple exon sets: {exon_set_n}')

    ####################################################################################################################
    # read one sequence at a time from the reference genome file, and process all exons multiple exon genes in that
    # sequence
    ####################################################################################################################
    splice = SpliceSite()
    jtotal = 0
    for sequence in genome:
        jcount = 0
        print(f'\nProcessing {sequence.id}')
        # id = sequence.id
        for exon_set in junction[sequence.id]:
            parent = exon_set[0].attribute['Parent']
            strand = exon_set[0].strand

            # main processing loop that adds counts to donor and acceptor PSSMs
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

        # uncommenting this limits input to just the first sequence for debugging
        # break

    ####################################################################################################################
    # Reports
    ####################################################################################################################
    print(f'{jtotal} splice junctions analyzed')

    # raw counts
    print(f'\n{'#' * 80}\n* counts\n{'#' * 80}')
    maxval = f'{int(np.max(splice.donor.matrix))}'
    PSSM.fieldwidth = len(maxval) + 2
    PSSM.precision = 0
    print(splice.donor)
    print(splice.acceptor)

    # frequencies
    print(f'\n{'#' * 80}\n* frequency\n{'#' * 80}')
    PSSM.fieldwidth = 7
    PSSM.precision = 3
    donor_freq = splice.donor.frequency()
    acceptor_freq = splice.acceptor.frequency()
    print(donor_freq)
    print(acceptor_freq)

    # positional information content
    print(f'\n{'#' * 80}\n# Positional Shannon entropy\n{'#' * 80}')
    I = donor_freq.information()
    J = acceptor_freq.information()
    # print(I)
    fmt = f'{PSSM.fieldwidth}.{PSSM.precision}f'
    print('\ndonor')
    fI = '\t '
    for i in I:
        fI += f'{i:{fmt}}'
    print(fI)
    print('\nacceptor')
    fJ = '\t '
    for j in J:
        fJ += f'{j:{fmt}}'
    print(fJ)

    # write donor and acceptor frequencies
    with open('donor.dat', 'w') as fq:
        print(f'{donor_freq}', file=fq)
    with open('acceptor.dat', 'w') as fq:
        print(f'{acceptor_freq}', file=fq)

    # print(f'Sharpened frequencies')
    # fq_from_file.sharpen(2)
    # print(fq_from_file)

    exit(0)
