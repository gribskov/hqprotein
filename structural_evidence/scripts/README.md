## Pipeline order

00\_wrapper.py               — Orchestrates all steps (TBD)
01\_parse\_gff.py             — Extract gene features and CDS coordinates from GFF3
02\_extract\_proteins.py      — Extract CDS from genome FASTA + translate to protein
03\_submit\_alphafold.py      — Sample proteins, generate AF3 input JSONs, submit jobs
04\_extract\_af3\_metrics.py   — Extract raw metrics from AF3 output (pLDDT, pTM, disorder, PAE)
05\_parse\_metrics.py         — Aggregate and summarize extracted metrics per protein
06\_validate\_quality.py      — Compare metric distributions across gene categories
07\_segment\_analysis.py      — Map confidence scores to individual exons
08\_calibrate\_bayesian.py    — Convert structural metrics into P(good model) for the Bayesian integrator





