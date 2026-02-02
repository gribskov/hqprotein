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

### Deliverables
- Cross-phylum benchmarking: does the pLDDT signal hold across diverse fungi or break in some lineages?
- Controlled corruption experiment: take gold-standard models, introduce specific errors (frameshift, missing exon, included intron, gene fusion, gene split), measure what pLDDT catches and what it misses
- Calibrated scoring: proper P(good model | pLDDT), ROC curves, not just binary better/worse
- Usable tool release: pipeline that takes GFF + genome and returns per-gene structural quality scores
- Biological insight: do certain gene categories (secreted, membrane, orphans) behave differently?

See Davison et al. 2025 (bioRxiv 10.1101/2025.10.21.683479) — they showed pLDDT works
for 2 fungi + 1 protist but only as proof-of-concept, no multi-evidence integration,
no controlled corruption, no cross-phylum analysis.

### Candidates 

| Organism | Phylum | Data quality |
|---|---|---|
| *Saccharomyces cerevisiae* | Ascomycota (Saccharomycotina) | SGD curated, AlphaFold proteome ready |
| *Schizosaccharomyces pombe* | Ascomycota (Taphrinomycotina) | PomBase, 100% annotated in UniProt |
| *Aspergillus fumigatus* | Ascomycota (Pezizomycotina) | >1000 curation events in FungiDB, used by Davison et al. |
| *Neurospora crassa* | Ascomycota (Pezizomycotina) | Broad + community curated, JGI MycoCosm |
| *Cryptococcus neoformans* | Basidiomycota (Agaricomycotina) | FungiDB + proteogenomic validation |
| *Candida albicans* | Ascomycota (Saccharomycotina) | FungiDB curated, AlphaFold proteome ready |

### To Do
- Figure out data source / organism
- Parse AlphaFold output (pLDDT per residue)
- Summarize into gene-level score
- Calibration experiment (good vs corrupted models)
