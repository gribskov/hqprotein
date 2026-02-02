## structural_evidence
Structural evidence pillar for the Bayesian gene model integrator.

### Idea
Take predicted protein sequences (good, mediocre, and bad gene models),
run them through AlphaFold, and see how pLDDT scores relate to model quality.
Bad models (frameshifts, stop codons, wrong boundaries) should fold poorly.
This relationship becomes our structural evidence signal for the Bayesian integrator.

### Steps
1. Collect protein models of varying quality — use internal stop codon
   count as proxy (no stops = good, many stops = bad, the more stop
   codons the harder for AlphaFold to fold)
2. Run through AlphaFold, extract pLDDT scores
3. Calibrate: pLDDT vs model quality
4. Define P(structure evidence | real gene)
5. Wrong reading frame confidence comes from the intrinsic signals
   pillar — both pillars meet in the Bayesian integrator

### To Do
- Figure out data source / organism
- Parse AlphaFold output (pLDDT per residue)
- Summarize into gene-level score
- Calibration experiment (good vs corrupted models)
