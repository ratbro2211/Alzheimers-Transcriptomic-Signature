# Validation report

Generated 2026-11-25 by `validation/scripts/05a-05d`. Verdicts: PASS = no problem; WARN = a caveat to report; FAIL = result not trustworthy; INFO = descriptive.

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

## Summary

| Dataset | PASS | WARN | FAIL | INFO |
|---|---|---|---|---|
| Microarray | 11 | 0 | 0 | 1 |
| RNA-seq | 8 | 3 | 0 | 0 |
| SEA-AD | 9 | 4 | 1 | 0 |
| SEA-AD cell-type | 4 | 3 | 3 | 1 |
| Comparison | 16 | 1 | 0 | 4 |

## Microarray

| ID | Check | Result | Verdict | Note |
|---|---|---|---|---|
| V1 | Title-derived labels match GEO 'diagnosis' field | 100% (AD 97 / ND 98) | **PASS** |  |
| V2 | Canonical AD markers in expected direction | 14/14 (binomial p = 1.2e-04) | **PASS** | tests that log2FC = AD - control, not reversed |
| V3a | Age balanced between AD and ND | median AD 86 vs ND 85, p = 0.89 | **PASS** | ages 90+ coded as 90 |
| V3b | Sex balanced between AD and ND | female AD 49% vs ND 49%, p = 1.00 | **PASS** |  |
| V4 | Age + sex adjustment changes little | r(log2FC) = 1.000; DEGs 724 -> 724; Jaccard 1.00 | **PASS** |  |
| V5 | p-value histogram: uniform tail + low-p spike | pi0 = 0.28; tail ratio (p>0.9 vs 0.5-0.6) = 0.99 | **PASS** | pi0 = estimated fraction of unchanged genes |
| V6 | Label-permutation null: DEGs under shuffled labels | observed 724; null median 0, max 5; 0/200 perms >= observed; empirical p = 0.005 | **PASS** | FAIL if empirical p > 0.05; WARN if any permutation reaches 50% of observed |
| V7a | PCA: condition captured on leading PCs | PC1 49% var; condition p on PC1-5: 3e-11, 9e-01, 2e-04, 1e-04, 4e-03 | **INFO** |  |
| V7b | PCA: sex/age not dominating PC1-2 | sex p PC1/PC2 = 0.33/0.085; age p = 0.16/0.086 | **PASS** |  |
| V7c | Outlier samples (Mahalanobis on PC1-5, p<0.001) | 0:  | **PASS** |  |
| V8 | Split-half reproducibility of log2FC (50 splits) | median r = 0.94 (range 0.84-0.97) | **PASS** | internal consistency of the AD signal |
| V9 | Probe choice: highest-mean vs most-significant probe | r = 0.960; DEGs 724 vs 761 (+5% with best-p probe) | **PASS** | most-significant-probe selection inflates DEG counts; highest-mean is unbiased |

## RNA-seq

| ID | Check | Result | Verdict | Note |
|---|---|---|---|---|
| V1 | Sample labels (AD vs age-matched Old) | AD 12 / control 10; Young excluded | **PASS** |  |
| V1b | Raw counts recovered exactly from GEO table | max deviation from integer = 2.1e-09 | **PASS** |  |
| V2 | Canonical AD markers in expected direction | 12/15 (neuronal down 10/10; glial up 2/5) | **PASS** |  |
| V3a | Age balanced between AD and control | median AD 66 vs control 68, p = 1.00 | **PASS** | sex not provided in GEO |
| V3b | Library size not confounded with condition | size factors 0.76-1.25; library size p = 0.02 | **WARN** |  |
| V4 | Age-adjusted model changes little | r(log2FC) = 0.979; DEGs 1152 -> 1195 | **PASS** |  |
| V5 | p-value histogram: uniform tail + low-p spike | pi0 = 0.57; tail ratio = 0.88 | **PASS** |  |
| V5b | Cook's distance outlier genes | 119 genes (0.7%); most often driven by 11-10T.RNA, 24-3T.RNA, 13-11T.RNA | **PASS** | DESeq2 replaces/flags these automatically |
| V6 | Label-permutation null: DEGs under shuffled labels | observed 1152; null median 0, max 583; 0/20 perms >= observed; empirical p = 0.048 | **WARN** | FAIL if empirical p > 0.05; WARN if any permutation reaches 50% of observed |
| V7 | Outlier samples (Mahalanobis PC1-3, p<0.01) | 1: 11-10T.RNA; condition p on PC1 = 0.25 | **WARN** |  |
| V8 | Leave-one-out: no single sample drives the result | min r = 0.961; DEG change -38% to +36%; most influential 16-14T.RNA | **PASS** |  |

## SEA-AD

| ID | Check | Result | Verdict | Note |
|---|---|---|---|---|
| V1 | Donor labels (cognitive status) | AD 10 / control 10 | **PASS** |  |
| V1b | Input is raw integer UMI counts | integer = TRUE; max per-donor gene count = 23665617 | **PASS** |  |
| V1c | Controls free of AD neuropathology | 5/10 controls have Intermediate/High ADNC (H20.33.014, H20.33.019, H20.33.025, H20.33.034, H21.33.025) | **WARN** | labels are clinical (dementia), not pathological; dilutes contrast |
| V2 | Canonical AD markers in expected direction | 14/16 (unshrunk log2FC); significant: 4 | **PASS** | apeglm shrinks most of these to ~0 |
| V3a | Age balanced | median AD 90 vs control 84, p = 0.50 | **PASS** |  |
| V3b | Sex balanced | female AD 5/10 vs control 7/10, p = 0.65 | **PASS** |  |
| V4 | Age + sex adjustment changes little | r(log2FC) = 0.971; DEGs 1520 -> 2547 | **PASS** |  |
| V5 | p-value histogram: uniform tail + low-p spike | pi0 = 0.55; tail ratio = 0.84 | **PASS** |  |
| V6a | Neuronal fraction balanced between groups | median AD 79% vs control 82%, p = 0.052 | **PASS** |  |
| V6b | DE log2FC explained by neuronal fraction | rho(log2FC, gene~neuronal-fraction corr) = -0.83 all genes, -0.44 DEGs | **WARN** | /rho/ > 0.5 = WARN; genome-wide fold changes strongly track neuronal fraction |
| V6c | Adjusting for neuronal fraction | r(log2FC) = 0.85; DEGs 1520 -> 82 | **WARN** | large change = composition drives much of the bulk signal |
| V7 | Label-permutation null: DEGs under shuffled labels | observed 1520; null median 2, max 5506; 2/20 perms >= observed; empirical p = 0.143 | **FAIL** | FAIL if empirical p > 0.05; WARN if any permutation reaches 50% of observed |
| V8 | Outlier donors (Mahalanobis PC1-3, p<0.01) | 0: ; PC1 vs neuronal fraction rho = 0.58 | **PASS** |  |
| V9 | Leave-one-out: no single donor drives the result | min r = 0.959; DEG change -82% to +95%; most influential H20.33.004 | **WARN** |  |

## SEA-AD cell-type

| ID | Check | Result | Verdict | Note |
|---|---|---|---|---|
| T1 | Permutation null + drop-2 donors: Excitatory IT | 379 DEGs (MLE 409); null median 0, max 38; perm p = 0.048; drop H20.33.028/038: r = 0.96, retained 16% | **PASS** | 10 AD / 10 control donors |
| T2 | Permutation null + drop-2 donors: Excitatory deep | 119 DEGs (MLE 133); null median 0, max 1588; perm p = 0.238; drop H20.33.028/038: r = 0.98, retained 42% | **FAIL** | 9 AD / 10 control donors |
| T3 | Permutation null + drop-2 donors: Pvalb/Chandelier | 41 DEGs (MLE 43); null median 0, max 12; perm p = 0.048; drop H20.33.028/038: r = 0.97, retained 21% | **PASS** | 10 AD / 10 control donors |
| T4 | Permutation null + drop-2 donors: Sst/Sst Chodl | 57 DEGs (MLE 59); null median 0, max 33; perm p = 0.048; drop H20.33.028/038: r = 0.96, retained 7% | **WARN** | 10 AD / 10 control donors |
| T5 | Permutation null + drop-2 donors: Vip/Lamp5/Sncg/Pax6 | 126 DEGs (MLE 142); null median 0, max 19; perm p = 0.048; drop H20.33.028/038: r = 0.96, retained 4% | **PASS** | 10 AD / 10 control donors |
| T6 | Permutation null + drop-2 donors: Astrocyte | 52 DEGs (MLE 53); null median 0, max 163; perm p = 0.095; drop H20.33.028/038: r = 0.94, retained 26% | **FAIL** | 10 AD / 10 control donors |
| T7 | Permutation null + drop-2 donors: Microglia/immune | 85 DEGs (MLE 90); null median 0, max 46; perm p = 0.048; drop H20.33.028/038: r = 0.93, retained 48% | **WARN** | 10 AD / 10 control donors |
| T8 | Permutation null + drop-2 donors: Oligodendrocyte | 109 DEGs (MLE 122); null median 0, max 18; perm p = 0.048; drop H20.33.028/038: r = 0.96, retained 60% | **PASS** | 10 AD / 10 control donors |
| T9 | Permutation null + drop-2 donors: OPC | 3 DEGs (MLE 3); null median 0, max 28; perm p = 0.190; drop H20.33.028/038: r = 0.87, retained 0% | **FAIL** | 10 AD / 10 control donors |
| T10 | Permutation null + drop-2 donors: Vascular | 12 DEGs (MLE 12); null median 1, max 7; perm p = 0.048; drop H20.33.028/038: r = 0.95, retained 42% | **WARN** | 10 AD / 9 control donors |
| T_conc | Agreement with microarray, by cell type (Pearson r) | Whole_tissue 0.27; Exc_IT 0.13; Exc_deep 0.11; Inh_Pvalb 0.06; Inh_Sst 0.03; Inh_CGE 0.06; Astrocyte 0.07; Microglia 0.06; Oligodendrocyte 0.06; OPC 0.01; Vascular 0.06 | **INFO** | which cell type's changes best match bulk tissue |

## Comparison

| ID | Check | Result | Verdict | Note |
|---|---|---|---|---|
| C1 | Gene-symbol matching across platforms | 13690 genes in all three; matched Microarray 69%, RNA-seq 77%, SEA-AD 45% | **PASS** | unmatched genes are mostly platform-specific (lncRNAs, unannotated probes) |
| C2 | Overlap above chance: Microarray vs RNA-seq | 102 shared vs 30.1 expected (3.4x), hypergeometric p = 1.4e-28 | **PASS** |  |
| C3 | Direction agreement vs coin-flip: Microarray vs RNA-seq | 101/102 (99%), binomial p = 4.1e-29 | **PASS** |  |
| C4 | Correlation robust (bootstrap CI, Spearman): Microarray vs RNA-seq | r = 0.40 [95% CI 0.38-0.42]; Spearman rho = 0.45 | **PASS** | Spearman is insensitive to apeglm shrinkage of magnitudes |
| C5 | Gene-shuffle null: Microarray vs RNA-seq | null r 95% range -0.017-0.017 (obs 0.40); null agreement 51% (obs 99%) | **PASS** |  |
| C2 | Overlap above chance: Microarray vs SEA-AD | 80 shared vs 35.3 expected (2.3x), hypergeometric p = 4.1e-12 | **PASS** |  |
| C3 | Direction agreement vs coin-flip: Microarray vs SEA-AD | 80/80 (100%), binomial p = 1.7e-24 | **PASS** |  |
| C4 | Correlation robust (bootstrap CI, Spearman): Microarray vs SEA-AD | r = 0.27 [95% CI 0.25-0.28]; Spearman rho = 0.33 | **PASS** | Spearman is insensitive to apeglm shrinkage of magnitudes |
| C5 | Gene-shuffle null: Microarray vs SEA-AD | null r 95% range -0.014-0.015 (obs 0.27); null agreement 49% (obs 100%) | **PASS** |  |
| C2 | Overlap above chance: RNA-seq vs SEA-AD | 76 shared vs 33.6 expected (2.3x), hypergeometric p = 1.7e-11 | **PASS** |  |
| C3 | Direction agreement vs coin-flip: RNA-seq vs SEA-AD | 72/76 (95%), binomial p = 3.6e-17 | **PASS** |  |
| C4 | Correlation robust (bootstrap CI, Spearman): RNA-seq vs SEA-AD | r = 0.25 [95% CI 0.24-0.27]; Spearman rho = 0.38 | **PASS** | Spearman is insensitive to apeglm shrinkage of magnitudes |
| C5 | Gene-shuffle null: RNA-seq vs SEA-AD | null r 95% range -0.015-0.015 (obs 0.25); null agreement 50% (obs 95%) | **PASS** |  |
| C6a | 33-gene panel vs random genes | panel 29/31 consistent; random genes mean 14.4, P(>= obs) = 0.0000 | **INFO** | random genes are consistent ~25% of the time by chance (2 of 8 sign patterns) |
| C6b | Panel vs random microarray-significant genes | random sig genes mean 24.6/31 consistent, P(>= 29) = 0.024 | **PASS** | random microarray-significant genes are already ~79% consistent: most of the panel's consistency reflects replication of significant genes in general |
| C7 | Drop most-influential sample: Microarray vs RNA-seq | r 0.40 -> 0.35; agreement 101/102 -> 57/58 | **WARN** | dropped: 16-14T.RNA, H20.33.004 |
| C7 | Drop most-influential sample: Microarray vs SEA-AD | r 0.27 -> 0.29; agreement 80/80 -> 115/116 | **PASS** | dropped: 16-14T.RNA, H20.33.004 |
| C7 | Drop most-influential sample: RNA-seq vs SEA-AD | r 0.25 -> 0.24; agreement 72/76 -> 92/99 | **PASS** | dropped: 16-14T.RNA, H20.33.004 |
| C8 | Effect-size ratio (/y///x/, sig in both): Microarray vs RNA-seq | median 0.97 | **INFO** | <1 = second dataset shows smaller fold changes (shrinkage, power, composition) |
| C8 | Effect-size ratio (/y///x/, sig in both): Microarray vs SEA-AD | median 1.87 | **INFO** | <1 = second dataset shows smaller fold changes (shrinkage, power, composition) |
| C8 | Effect-size ratio (/y///x/, sig in both): RNA-seq vs SEA-AD | median 1.37 | **INFO** | <1 = second dataset shows smaller fold changes (shrinkage, power, composition) |

## Figures

- `figures/pvalue_histograms.png`, `figures/permutation_null.png`, `figures/leave_one_out.png`,
  `figures/panel_circularity.png`, `figures/pca_microarray.png`, `figures/pca_rnaseq.png`, `figures/pca_seaad.png`
