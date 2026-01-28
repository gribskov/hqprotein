# hqprotein
## improved protein gene models

Anecdotally it seems that may gene models in complete genomes contain errors 
such as missing exons, included introns, artificial fusion to adjacent genes,
and artificial splitting into multiple predicted genes. Thee is a great deal of 
information that is not included in standard gene models that might be able to detect and correct these errors

## Possible information to include
For simplicity, we can discuss coding reading frames, **cORFs**, noncoding reading
frames that overlap a cORF, **oORFs**, noncoding reading frames such as introns
or UTRs, **nORFs**, and intergenic reading frames, **iORFs**.
* Amino acid composition - cORFs have a typical composition which is very 
different from the composition of oORFs, nORFs, or iORFs. 
or noncoding regions(introns, intergenic)
* codon preference - cORFs have a preference for specific
codons within a synonymous family, again this is different from 
oORFs, nORFs, or iORFs.
* codon usage - the frequency of each codon, combines amino acid 
composition and codon preference.
* DNA kmer composition - usually kmer composition is the basis for intron/exon 
prediction based on a hexamer (fifth order) hidden markov model. As this is
included in the base prediction, it may be superfluous. Would not distinguish 
cORFs, oORFs, nORFs, and iORFs.
* Protein kmer - Dipeptide composition should be different for cORFS and
oORFs, nORFs, or iORFs. Hexamer based models should already include dipeptide
frequencies and codon preference.
* 3-base periodicity. Fickett showed that coding regions have a distinct
three base periodicity due to the stronger  constraints on the first and
second bases of the codon. This allows coding regions (cORF and oORF) 
to be distinguished from noncoding regions (nORF and iORF).
  > Fickett JW. Recognition of protein coding regions in DNA sequences. Nucleic Acids Res. 1982 Sep 11;10(17):5303-18. doi: 10.1093/nar/10.17.5303.
* Splice donor/acceptor sites - use to identify possibly spliced regions
with ORFs. The canonical GT|AG found at the beginning|end of introns should
work, but what abut noncanonical splice sites?  GC|AG and AT|AC make 
up 0.89% and 0.10% of splices sites in human, respectively (species specific). 
  > Parada GE, Munita R, Cerda CA, Gysling K. A comprehensive survey of non-canonical splice sites in the human transcriptome. Nucleic Acids Res. 2014;42(16):10564-78. doi: 10.1093/nar/gku744.

## Approach
My hand-wavy idea is to calculate a posterior probability of coding on a 
base-by-base or codon-by-codon basis to identify ORFs or parts of ORFs that are
most likely to be coding.
Another possibility would be to train an ML model.