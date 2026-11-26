# AD classifier on the cross-platform gene panel — results and interpretation

Script: `scripts/06_ml_classification.R` (seed 42; full log in `results/06_run_log.txt`).

## Setup
- **Features:** 31 of the 33 panel genes. CAMK2A and PENK have no annotated probe on the microarray, so they cannot be used by a model trained there.
- **Expression:** log2 scale in every dataset.
  - Microarray: log2 intensity.
  - RNA-seq: log2(DESeq2-normalised count + 1), from raw counts recovered exactly from the GEO table.
  - SEA-AD: log2(DESeq2-normalised pseudobulk UMI + 1).
- **Scaling:** each gene is z-scored *within each dataset*, which uses no labels. Microarray intensities and sequencing counts are on different scales, so a model cannot transfer without this. In LOO-CV the scaling is fitted on the training fold only.
- **Models:** L2 logistic regression (glmnet, lambda chosen by inner CV) and random forest (500 trees).
- **Threshold:** chosen once on microarray LOO-CV predictions so that sensitivity ≥ 0.85, then frozen. The RNA-seq and SEA-AD sets were never used for fitting, tuning or threshold choice.
- **Sample sizes:**
  - Microarray: 97 AD / 98 control.
  - RNA-seq: 12 AD / 10 age-matched controls.
  - SEA-AD: 10 dementia / 10 no dementia.

## Results

| Dataset | Model | AUC [95% CI] | Sensitivity | Specificity | Balanced acc. |
|---|---|---|---|---|---|
| Microarray, LOO-CV (n = 195) | Logistic L2 | **0.87** [0.82–0.92] | 0.86 | 0.72 | 0.79 |
| | Random forest | **0.86** [0.81–0.91] | 0.86 | 0.74 | 0.80 |
| RNA-seq, independent (n = 22) | Logistic L2 | **0.90** [0.77–1.00] | 0.92 (11/12) | 0.70 (7/10) | 0.81 |
| | Random forest | **0.95** [0.86–1.00] | 0.92 (11/12) | 0.60 (6/10) | 0.76 |
| SEA-AD, independent (n = 20) | Logistic L2 | 0.63 [0.36–0.90] | 0.90 (9/10) | 0.50 (5/10) | 0.70 |
| | Random forest | 0.63 [0.37–0.89] | 0.50 (5/10) | 0.50 (5/10) | 0.50 |

Sensitivity, specificity and balanced accuracy are at the frozen threshold; results at 0.5 are in `results/crossdataset_results.csv`.

- **Microarray and RNA-seq:** AUCs are well above chance (Wilcoxon p ≈ 1e-18 and ≤ 0.001).
- **SEA-AD:** not distinguishable from chance (AUC CI includes 0.5; p = 0.34–0.35).
- **Learning curve:** LOO-CV AUC rises from ~0.80 (n = 39) to ~0.86 (n = 117) and then levels off (0.85–0.87 at 156–195). More microarray samples would add little; the limit is the information in 31 genes, not sample size.
- **Feature importance:** VGF ranks 1st (RF) and 3rd (LR), and GFAP 4th and 2nd, so both expected markers rank top. SNAP25 ranks 2nd in RF but 21st in LR. SST ranks 3rd in RF but last in LR. NRGN, NPY and NPTX2 are also consistently high.
  - Why the two models disagree: the neuronal genes are highly correlated, so L2 regression spreads weight across them. Some coefficients flip sign (VSNL1, TAC1, SYT4 fall in AD but get positive coefficients). Read LR coefficients as *conditional* weights, not as per-gene directions. The column `marginal_direction_in_AD` gives the true direction.

## What the cross-dataset performance means
The microarray model transferred well to RNA-seq, a different platform, brain region (lateral temporal lobe) and cohort, catching 11 of 12 AD cases. The signal the panel captures (lower neuronal/neuropeptide genes, higher GFAP) is therefore not a microarray artefact.

It did **not** transfer to SEA-AD whole-tissue pseudobulk. This matches the validation report. SEA-AD's whole-tissue AD signal is weak, depends on which cells each donor contributed, and hinges on a few glia-heavy donors (validation V6–V9). SEA-AD groups are also defined by cognition: 5 of the 10 "controls" carry intermediate or high AD pathology, so many controls are biologically AD-like. A failure here says more about the SEA-AD contrast than about the model.

**Clinical framing (important):** every dataset is post-mortem brain tissue. The classifier separates AD from control *brains*; it is not a diagnostic test for living patients, and its sensitivity and specificity are not clinical performance figures. The sensitivity ≥ 0.85 threshold shows how the trade-off behaves on transfer. It reached 0.92 sensitivity on RNA-seq at a cost of 60–70% specificity.

## Limitations
1. **Panel selection is circular.** The 33 genes were chosen from earlier analyses because their direction agreed across *all three* datasets. Within-microarray LOO-CV is optimistic, and so is the RNA-seq "independent" test: these samples helped choose the features, though they never touched model fitting or the threshold. A clean test needs genes chosen on the microarray alone, or a fourth cohort (e.g. GSE5281, Mayo/ROSMAP).
2. **Small test sets.** With n = 22 and n = 20, AUC 95% CIs span about ±0.13–0.27, and one sample moves sensitivity or specificity by 8–10 points.
3. **Per-dataset z-scoring needs a batch of samples** and assumes a broadly similar case mix (here near-balanced). It cannot score a single new sample, and an extreme case/control ratio would shift the scaling.
4. **Composition, not only regulation.** The top features (VGF, NRGN, SNAP25, GFAP) mark neuron loss and gliosis. The cell-type analysis showed that much of their tissue-level change is compositional. The classifier mostly detects tissue that has lost neurons, which is real AD biology but not a cell-intrinsic molecular signature.
5. **Region and platform differences** (MTG vs LTL; array vs sequencing) were handled only by within-dataset scaling.

## Paragraph for a graduate-school application
> To test whether a reproducible transcriptomic signature of Alzheimer's disease generalises across measurement technologies, I trained regularised logistic regression and random-forest classifiers on a 31-gene panel using 195 post-mortem middle temporal gyrus microarrays (leave-one-out AUC 0.87). Without any retraining, the frozen model classified an independent RNA-seq cohort from a different temporal-lobe region with AUC 0.90–0.95, identifying 11 of 12 AD brains. It failed, however, on single-nucleus pseudobulk data (AUC 0.63). Tracing this failure showed that the panel's tissue-level signal is driven largely by changes in cell-type composition (neuron loss and gliosis) rather than by per-cell expression changes, and that the single-nucleus cohort's cognitively defined controls often carried AD pathology. A learning curve showed that performance plateaued at about 120 training samples, pointing to feature information, not sample size, as the limit. The project taught me to treat generalisation failures as data, and to separate what a tissue is made of from how its cells behave, a distinction I want to model quantitatively in complex biological systems.

## Files
- `results/feature_matrices/*.csv`: samples × 31 genes (log2, before z-scoring) + `condition` (AD = 1).
- `results/loocv_results.csv`, `results/loocv_predictions.csv`
- `results/crossdataset_results.csv`, `results/crossdataset_predictions.csv`
- `results/learning_curve.csv` (+ `_all_repeats.csv`)
- `results/feature_importance.csv`
- `figures/fig1_roc_curves.pdf`, `fig2_learning_curve.pdf`, `fig3_feature_importance.pdf`, `fig4_confusion_matrices.pdf`
