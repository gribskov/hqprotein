# Conservation evidence

## Approach
Use DNA segments that match to real proteins from the BlastX search, and calculate the probability of that open reading
frame being a true ORF. 

## Steps

1. Use input example data, for example the unannotated nucleotide sequences (scaffolds, contigs) of *Z. tritici* and BlastX against the uniref50 protein database.
Use a comprehensive database like Uniprot50.

**Logic:** BlastX Gives Reading Frame Evidence.Each BlastX hit tells which reading frame that region of DNA is coding in. I can then associate those hits with candidate ORFs in the DNA.

3. Process the BlastX results to get the associations of correct ORF

**Logic:** ORFs Have Probabilities. Not every ORF is real. For each candidate ORF I can ask:

- How long is it?
- Does BlastX support it?
- What's the probability this ORF encodes a real protein by chance vs. being genuinely coding?

4. Combine the BlastX evidence with an ORF probability score to rank which ORFs are likely real.

5. Combine evidence across exons

Instead of judging each exon independently:

Script 1: Extract subsequences between splice junctions (each exon) and Generate all possible combinations of reading frames and exon pairings

Script 2: Treat the splice junctions as a system: what combination of frames and exon boundaries makes the overall multi-exon protein match most probable?

Script 3: Assign a joint probability to each combination

This turns the problem from "does this exon match a protein?" into "what gene model, considered as a whole, best explains all the BlastX evidence across all exons?"

Note: Use the latest high-quality *Z. tritici* annotation to confirm the predicted ORF.

## Ideas
A protein match that is split across multiple frames, likely indicated a framing error in the assembly.

BLASTX results should align with predicted exon boundaries, specifically look for splice sites.
