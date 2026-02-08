"""=====================================================================================================================
splice_model.py

frequency model of splice donor and acceptor based on genome and annotation (GFF3)
frequencies will be P(base(pos)|acceptor) and P(base(pos)|donor)

Michael Gribskov 2/6/2026
====================================================================================================================="""
from include.gff2 import GxfSet, GxfRecord

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

        if feature_count > 50: break

    print(f'\nexons: {feature_count}')
    print(f'multiple exons sets: {len(junction)}')

    for j in junction:
        # TODO this needs to consider strand
        print(f'{j[0].attribute['Parent']} strand: {j[0].strand}')
        donor = j[0]
        for exon in range(1,len(j)):
            acceptor = j[exon]
            print(f'\tdonor {donor.end}\tacceptor {acceptor.start}\t')
            donor = acceptor

    exit(0)
