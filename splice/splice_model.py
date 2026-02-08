"""=====================================================================================================================
splice_model.py

frequency model of splice donor and acceptor based on genome and annotation (GFF3)
frequencies will be P(base(pos)|acceptor) and P(base(pos)|donor)

Michael Gribskov 2/6/2026
====================================================================================================================="""
from include.gff.gff2 import GxfSet, GxfRecord
from include.sequence.fasta import Fasta
from collections import defaultdict

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
