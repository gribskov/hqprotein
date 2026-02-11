"""=====================================================================================================================
splice_model.py

frequency model of splice donor and acceptor based on genome and annotation (GFF3)
frequencies will be P(base(pos)|acceptor) and P(base(pos)|donor)

Michael Gribskov 2/6/2026
====================================================================================================================="""
from include.gff.gff2 import GxfSet
from include.sequence.fasta import Fasta
from collections import defaultdict


class SpliceSite:
    """=================================================================================================================
    Frequency matrix for splice donor and acceptor
    ================================================================================================================="""

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
        donor_n     count of donor examples; use to calculate frequencies
        acceptor_n  count of acceptor examples; use to calculate frequencies
        pre         bases before site (exon side)
        post        bases after site (intron side)
        -------------------------------------------------------------------------------------------------------------"""
        self.donor = []
        self.donor_n = 0
        self.acceptor = []
        self.acceptor_n = 0
        self.pre = pre
        self.post = post

        # initialize donor and acceptor
        self.site_init()

    def site_init(self):
        """-------------------------------------------------------------------------------------------------------------
        initialize the count matrices (donor and acceptor), and the total counts (donor_n, acceptor_n)

        :return: True
        -------------------------------------------------------------------------------------------------------------"""
        sitesize = self.pre + self.post
        self.donor = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for _ in range(sitesize)]
        self.acceptor = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for _ in range(sitesize)]
        self.donor_n = self.acceptor_n = 0

        return True

    def add_junction(self, donorpos, acceptorpos, sequence):
        """-------------------------------------------------------------------------------------------------------------
        add a donor/acceptor pair from sequence based con coordinates (from GFF)

        :param donorpos: int        end of donor exon
        :param acceptorpos: int     beginning of acceptor exon
        :param sequence: str        DNA sequence
        :return: int, int           number of donor/acceptor sequences
        -------------------------------------------------------------------------------------------------------------"""
        self.donor_n += 1
        self.acceptor_n += 1
        sitepos = 0
        for pos in range(donorpos - self.pre, donorpos + self.post):
            self.donor[sitepos][sequence[pos]] += 1
            sitepos += 1

        sitepos = 0
        for pos in range(acceptorpos - self.post, acceptorpos + self.pre):
            self.acceptor[sitepos][sequence[pos]] += 1
            sitepos += 1

        return self.donor_n, self.acceptor_n


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

    feature_n = gff.feature_get(['exon'])
    junction = defaultdict(list)
    feature_count = 0
    current = ''
    exon_set = []
    exon_set_n = 0
    for feature in gff.features:
        # find exons from multi exon genes and store in junction
        feature_count += 1
        print(f'{feature.start}\t{feature.end}\t{feature.strand}\t{feature.attribute['Parent']}')
        if current == feature.attribute['Parent']:
            exon_set.append(feature)
        else:
            if len(exon_set) > 1:
                # multiple exons, get acceptor and donor sites
                print(f'\tmultiple: {exon_set}')
                sequence = exon_set[0].seqid
                junction[sequence].append(exon_set)
            exon_set = [feature]
            exon_set_n += 1

        current = feature.attribute['Parent']

        if feature_count > 50: break

    print(f'\nexons: {feature_count}')
    print(f'multiple exons sets: {exon_set_n}')

    # read sequence in genome
    for sequence in genome:
        print(f'{sequence['id']}')
        id = sequence['id']
        for exon_set in junction[sequence['id']]:
            parent = exon_set[0].attribute['Parent']
            strand = exon_set[0].strand

            if strand == '-':
                continue

            print(f'{parent} strand: {strand}')
            for i in range(len(exon_set) - 1):
                donor = exon_set[i].end
                acceptor = exon_set[i + 1].start
                print(f'\tdonor {donor} {sequence['seq'][donor - 5:donor + 10]}\t'
                      f'acceptor {acceptor} {sequence['seq'][acceptor - 10:acceptor + 4]}')

    exit(0)
