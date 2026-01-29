"""=====================================================================================================================
gene_regions_from_gtf.py

Extract regions from a fasta file with extra sequence on each end. Intended for getting gene sequences with edges but
should work for other features

2026-01-15 gribskov
====================================================================================================================="""
from collections import defaultdict
from fasta import Fasta
from include.gff2 import *

####################################################################################################
# Main
####################################################################################################
if __name__ == '__main__':
    padding = 500
    gffin = 'data/z.tritici.IP0323.reannot.gff3'
    gff = GxfSet(file=gffin, fmt='gff')
    feature_n = gff.feature_get(['gene'])
    print(f'{feature_n} features read from {gffin}')

    genes = defaultdict(list)
    for f in gff.features:
        # padded coordinate may be beyond the range of the chromosome
        begin = f.start - padding  # doesn't occur in test data so don't check for negative
        end = f.end + padding
        if begin <= 0:
            print(f'{f.seqid}\t{f.start}\t{f.end}\t{f.strand}\t{f.attribute['ID']}\t{begin}\t{end}')
        genes[f.seqid].append(f)
    print()

    fasta = Fasta('data/z.tritici.IP0323.fasta')
    while fasta.next():
        print(f'processing {fasta.id}')
        g = genes[fasta.id]

        for thisgene in g:
            fastaout = Fasta()
            fastaout.id = thisgene.attribute['ID']
            fastaout.doc = f'padding={padding} padded_loc=[{thisgene.start},{thisgene.end}] '
            fastaout.doc += f'strand={thisgene.strand}'
            seq = fasta.seq[thisgene.start-padding-1:thisgene.end+padding]
            if thisgene.strand == '-':
                seq = Fasta.reverseComplement(seq)
            fastaout.seq = seq
            print(fastaout.format(100))


    exit(0)
