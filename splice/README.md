## Splice
model splice sequences in genome based on reference annotation

PSSM.py<br>
General object for storing position specific scoring matrices and scanning for matches in a 
sequence. 
Status: Tested and working 2026 March 22
* calculates frequencies as counts
* converts counts to probabilities
* logtransforms counts
* calculates positional Shannon entropy
* sharpens frequencies by raising to power
* writes and reads PSSMs
  * optional format control for values

splice_model.py<br>
give a gff and genome sequences, extracts transcripts and tabulates base frequencies from position
-5 to +10.
* writes donor and acceptor pssms for use in evaluating sequence regions

splice_scan.py<br>
uses donor and acceptor PSSM models to evaluate the probability that each position in a sequence
is a splice site
