## Overall interpretation

**Microarray (GSE132903): sound.** The labels are correct, the contrast runs AD − control (14/14 known markers in the expected direction), and age and sex are balanced. Adjusting for them changes nothing (r = 1.00). Shuffled labels give almost no DEGs (0 in 199 of 200 permutations), so the FDR is well controlled. Split halves reproduce each other (median r = 0.94), and there are no outlier samples. A low pi0 (0.28; most genes shifted) matches PC1 tracking diagnosis (49% of variance, p = 3e-11), which points to a broad tissue-level change such as neuronal loss. This dataset is the anchor of the comparison.

**RNA-seq (GSE104704): sound, with caveats.** The raw counts were recovered exactly and the labels are correct. All 10 neuronal markers go down, but only 2 of 5 glial/inflammatory markers go up, so the glial signal is weak in this cohort. Age is balanced, and adjusting for it changes little (r = 0.98). Two WARNs:
- **Library size** differs between groups (p = 0.02). DESeq2 normalisation handles this.
- **The permutation test** reached the lowest p-value 20 permutations allow (0.048), and one shuffle gave 583 DEGs (51% of observed). With 22 samples some random splits partly reproduce the real structure.

Leaving out any single sample keeps fold changes at r ≥ 0.96 with the full result, so no single sample drives this result.

**SEA-AD (snRNA-seq pseudobulk): direction trustworthy, DEG list not.** The inputs are correct: raw integer UMIs, correct labels, and balanced age and sex. Three findings limit what can be claimed:
1. **Permutation test (FAIL).** 2 of 20 random label shuffles produced as many or more DEGs than the true labels (empirical p = 0.14). The *number and identity* of SEA-AD DEGs cannot be distinguished from those produced by arbitrary donor splits.
2. **Cell composition (WARN).** Genome-wide fold changes track each gene's correlation with neuronal fraction (rho = −0.83). Adjusting for neuronal fraction cuts the DEG count from 1,520 to 82. Much of the whole-tissue signal therefore reflects which cells were sampled.
3. **Influential donors (WARN).** Removing either of the two glia-heavy AD donors (H20.33.028 or H20.33.038, 54% and 26% neuronal) removes 70–82% of DEGs.

Five of the ten "controls" also carry intermediate or high AD neuropathology, which dilutes the contrast.

**Cross-platform comparison: concordance is real and not an artefact.**
- **Not chance:** overlaps of significant genes are 2.3–3.4× above chance (p ≤ 1e-11), and direction agreement is 95–100% where a coin flip gives 50%.
- **Not a gene-matching artefact:** genome-wide correlations sit far outside a gene-shuffle null (null r ≈ ±0.015 vs observed 0.25–0.40), with tight bootstrap CIs. Rank (Spearman) correlations are higher, so shrinkage does not create the agreement.
- **Survives dropping influential samples:** agreement stays 94–99% when the most influential sample in each small dataset is removed. Microarray–RNA-seq r drops 0.40 → 0.35 (flagged WARN), but direction agreement stays at 98%.

So SEA-AD's individual DEG calls are unreliable, but its *direction* of change agrees with the two independent bulk datasets far beyond chance.

**33-gene panel: partly circular.** The panel was hand-picked from earlier results, so its consistency is not independent evidence. Random genes are direction-consistent in only 14/31 cases, but random *microarray-significant* genes already reach 24.6/31. The panel (29/31) is only modestly better (p = 0.024). Present it as an illustration of the reproducible signature, not as a discovery.

### What this means for the write-up
1. The headline claim is supported: AD expression changes replicate in direction across platforms and regions (C2–C5 PASS). Cite the microarray–RNA-seq and microarray–SEA-AD agreement.
2. Do not report SEA-AD DEG counts or top SEA-AD genes as findings; use SEA-AD only for direction and cell-type context.
3. State that whole-tissue SEA-AD changes are largely compositional, and point to cell-type-resolved pseudobulk as the correct next analysis.
4. Describe the 33 genes as an illustrative, direction-consistent panel, not as an unbiased gene set.

### Follow-up: cell-type-resolved SEA-AD pseudobulk (`dataset3_seaad/scripts/03c–03d`)
UMIs were summed per donor *within* each of 10 cell types, and each cell type was tested separately with the same permutation and drop-donor checks.
- **Partly fixes the whole-tissue FAIL.** In 7 of 10 cell types (excitatory IT, Pvalb, Sst, Vip/Lamp5/Sncg/Pax6, microglia, oligodendrocytes, vascular), no shuffle produced as many DEGs as the true labels (p = 0.048, the floor for 20 shuffles). Excitatory deep neurons, astrocytes and OPCs still do not beat chance. DEG counts are modest (3–379 per cell type), and most are up-regulated.
- **DEG identity is fragile.** Dropping the two glia-heavy AD donors keeps the fold-change pattern (r = 0.87–0.98) but only 4–60% of the DEGs. Cell-type DEG lists from 10 vs 10 donors should be treated as exploratory.
- **The whole-tissue changes in the panel genes are compositional.** Within astrocytes, GFAP (+0.03), VIM (−0.04) and GJA1 (+0.43) barely change, against +1.3, +2.0 and +2.7 in whole tissue. TREM2 within microglia is +0.14, against +1.2 in whole tissue. SST within Sst neurons is −0.35, against −1.46 in whole tissue. The bulk "glial up / Sst down" signal is mainly more glia and fewer Sst neurons, not changed expression per cell.
- **Some neuronal genes do fall within neurons.** CARTPT (−1.14), NPTX2 (−0.74), CHGB (−0.57) and TAC1 (−0.54) decrease inside their main neuron type by about as much as in whole tissue (not significant at n = 10/10). This is a candidate cell-intrinsic component.
- **Bulk tissue matches whole-tissue pseudobulk best.** Agreement with the bulk datasets is highest for whole-tissue SEA-AD (rho 0.33/0.38) and much lower for any single cell type (≤ 0.22, best excitatory IT neurons). The bulk AD signature is largely a cell-composition signature. This is a result in its own right: it explains why AD transcriptomic changes replicate across platforms, and why they partly reflect neuron loss and gliosis.
