"""=====================================================================================================================
# codon_model.py
#
# Michael Gribskov 1/28/2026
====================================================================================================================="""
from include.fasta import Fasta
from codon import Codon

dnafile = 'data/z.tritici.IP0323.reannot.cds.fasta'
cds = Fasta(filename=dnafile)

codon = Codon()

codon_total = 0
while cds.next():
    print(cds.id)
    codon_total += codon.add_from_dna(cds.seq)
    print(f'\tnumber of codons => {codon_total}')

    freq_sum = codon.update_frequencies()
    print(f'\tfrequency sum => {freq_sum}')

exit(0)
