## Pipeline order

00\_wrapper.py            — Orchestrates all steps (TBD)
01\_parse\_gff.py          — Extract gene features and CDS coordinates from GFF3
02\_extract\_proteins.py   — Extract CDS from genome FASTA + translate to protein
03\_submit\_alphafold.py   — Generate AlphaFold input, submit jobs, collect results
04\_parse\_cif.py          — Extract per-residue pLDDT from AlphaFold CIF output
05\_aggregate\_plddt.py    — Compute per-protein pLDDT statistics
06\_validate\_quality.py   — Compare pLDDT distributions between good and bad models
07\_segment\_analysis.py   — Map pLDDT scores to individual exons
08\_calibrate\_bayesian.py — Convert pLDDT into P(good model | pLDDT) for the Bayesian integrator





