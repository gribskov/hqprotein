Maybe a little off target. This code uses the output from orthofinder to decide which of multiple transcripts is the 
most representative. Briefly the idea is to choose the transcipt that is most central to the orthogroup based on the 
MCL distance between the members of the OG.

1. use WorkingDirectory/SequenceIDs to map sequence IDs from individual genomes onto the numeric ID and species gene ID