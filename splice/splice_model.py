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
    base_complement = str.maketrans('ACGT', 'TGCA')

    def __init__(self, pre=5, post=10, fieldwidth=4, precision=0):
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
        self.fieldwidth = fieldwidth
        self.precision = precision

        # initialize donor and acceptor
        self.site_init()

    def site_init(self):
        """-------------------------------------------------------------------------------------------------------------
        initialize the count matrices (donor and acceptor), and the total counts (donor_n, acceptor_n)

        :return: True
        -------------------------------------------------------------------------------------------------------------"""
        sitesize = self.pre + self.post
        self.donor = [{'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0} for _ in range(sitesize)]
        self.acceptor = [{'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0} for _ in range(sitesize)]
        self.donor_n = self.acceptor_n = 0

        return True

    def add_junction(self, strand, donorpos, acceptorpos, sequence):
        """-------------------------------------------------------------------------------------------------------------
        add a donor/acceptor pair from sequence based on coordinates (from GFF). Input sequences are forced to
        uppercase. Currently, no checking for non A,C,G,T bases.

        :param donorpos: int        end of donor exon
        :param acceptorpos: int     beginning of acceptor exon
        :param strand: strand       strand, '+' or '-'
        :param sequence: str        DNA sequence
        :return: int, int           number of donor/acceptor sequences
        -------------------------------------------------------------------------------------------------------------"""
        sequence = sequence.upper()
        self.donor_n += 1
        self.acceptor_n += 1
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

        for jtype in ('donor', 'acceptor'):
            sitepos = 0
            pssm = getattr(self, jtype)

            site = sequence[left[jtype]:right[jtype]]
            if strand != '+':
                # reverse complement - strand
                site = site.translate(SpliceSite.base_complement)[::-1]

            for base in site:
                pssm[sitepos][base] += 1
                sitepos += 1

        return self.donor_n, self.acceptor_n

    def frequency(self):
        """------------------------------------------------------------------------------------------------------------
        Return a new SpliceSite object with counts converted to frequencies.

        :return: SpliceSite     positional probability of base | splice_site
        -------------------------------------------------------------------------------------------------------------"""
        frequency = SpliceSite()
        for jtype in ('donor', 'acceptor'):
            freqpos = 0
            pssm = getattr(self, jtype)
            freq = getattr(frequency, jtype)
            freq_n = 0.0

            # the frequency of bases at each position is calculated and stored in donor_n and acceptor_n. Currently,
            # this is overwritten with 1, but it could be used to get the average number of observations per column
            for pos in range(len(pssm)):
                total = sum(pssm[pos].values())
                print(total)

                for base in pssm[pos]:
                    freq[freqpos][base] = pssm[pos][base] / total
                    freq_n += freq[freqpos][base]

                freqpos += 1

            # After conversion, the frequency of all columns is 1
            setattr(frequency, jtype + '_n', 1)

        return frequency

    @staticmethod
    def complement(sequence):
        """-------------------------------------------------------------------------------------------------------------
        complements but does not reverse the sequence string A->T, C->G, G->C, T->A. No checking for ambiguous or
        incorrect bases

        :param sequence: str    DNA sequence, 5' to 3'
        :return: str            complementary strand in 3' to 5' order
        -------------------------------------------------------------------------------------------------------------"""

    def __str__(self):
        """-------------------------------------------------------------------------------------------------------------
        formatted version of site. since __str__() does not accept parameters, the format is included in the object
        itself (self.fmt)

        :return: str    formatted string with donor and accptor sites
        -------------------------------------------------------------------------------------------------------------"""
        fmt = f'{self.fieldwidth}.{self.precision}f'
        divider = f'{'|':>{self.fieldwidth}s}'
        dpos = {'donor': self.pre - 1, 'acceptor': self.post - 1}
        outstr = ''
        for site in ('donor', 'acceptor'):
            outstr += f'{site}:\n'
            sitedata = getattr(self, site)
            for base in sitedata[0]:
                outstr += f'\t{base}  '
                for pos, column in enumerate(sitedata):
                    outstr += f'{column[base]:{fmt}}'
                    if pos == dpos[site]:
                        outstr += divider
                outstr += '\n'

        return outstr


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

        # TODO remove debug
        if feature_count > 500: break

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

        # TODO remove debug
        break

        splice.fieldwidth = 7
        print(f'\n{splice}')

    print(f'{jtotal} splice junctions analyzed')

    frequency = splice.frequency()
    frequency.fieldwidth = 6
    frequency.precision = 3
    print(f'\n{frequency}')

    exit(0)
