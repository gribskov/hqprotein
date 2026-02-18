# Conservation evidence

## Approach
Use DNA segments that match to real proteins from the BlastX search, and calculate the probability of that open reading
frame of being a true ORF. 

## Steps

1. Use input example data, for example the unannotated nucleotide sequences (scaffolds, contigs) of *Z. tritici* and BlastX against protein database.
Use a comprehensive database like Uniprot.
2. Process the BlastX results to get the associations of correct ORF 

## Ideas
A protein match that is split across multiple frames, likely indicated a framing error in the assembly.

BLASTX results should align with predicted exon boundaries, specifically look for splice sites.
