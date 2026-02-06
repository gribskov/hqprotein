"""=====================================================================================================================
splice_model.py

frequency model of splice donor and acceptor based on genome and annotation (GFF3)
frequencies will be P(base(pos)|acceptor) and P(base(pos)|donor)

Michael Gribskov 2/6/2026
====================================================================================================================="""
from include.gff.gff2 import GxfSet, GxfRecord

# ======================================================================================================================
# main
# ======================================================================================================================
if __name__ == '__main__':
    gff_file = 'data/z.tritici.IP0323.reannot.gff3'
    gff = GxfSet(file=gff_file, fmt='gff')

    feature_n = gff.feature_get(['exon'])
    junction = []
    feature_count = 0
    current = ''
    exon_set = []
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
                junction.append(exon_set)
            exon_set = [feature]

        current = feature.attribute['Parent']

        # if feature_count > 50: break

    print(f'\nexons: {feature_count}')
    print(f'multiple exons sets: {len(junction)}')

    exit(0)
