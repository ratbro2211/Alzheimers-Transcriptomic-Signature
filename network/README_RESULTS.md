# Network analysis of the 33-gene AD panel

Network-level analysis of the panel genes, built on the differential-expression outputs in `../analysis/`.
Run from `network/`:
`Rscript scripts/p2_00_prepare_expression.R` → `python scripts/p2_01_string_network.py` → `p2_02`, `p2_03`, `Rscript scripts/p2_04_wgcna.R`, `p2_05`, `p2_06` → `p2_07`.

**Terminology:** the STRING graph is a *functional-association / protein-interaction network*, not a gene regulatory network. Regulatory (TF → target) edges come only from TRRUST (section 3).

### Inputs
- **Network files** (`data/`): STRING v12.5 edges and ID map for the 33 genes, and TRRUST v2 (human).
- **`p2_00` rebuilds the expression matrices,** because `analysis/*/results/` hold DE tables, not per-sample matrices. It is written in R because probe mapping needs `illuminaHumanv4.db` and normalisation needs DESeq2.
  - Microarray: log2 intensity, highest-mean probe per gene (19,743 genes).
  - RNA-seq: log2(DESeq2-normalised + 1) from raw counts. AD 12 vs age-matched Old 10.
  - SEA-AD: log2(DESeq2-normalised UMI + 1) from the donor pseudobulk. 10 vs 10, grouped by cognitive status.
- **WGCNA** uses the 5,000 most variable microarray genes plus the panel. A single network of all ~20k genes does not fit in 8.6 GB of RAM.
- **TF enrichment** uses Fisher's exact test against two backgrounds: all microarray genes, and the stricter TRRUST targets on the array.
- **Cell-type resolution** reuses the cell-type DESeq2 results in `analysis/dataset3_seaad/results/celltype/` (10 SEA-AD cell classes).

## 0. STRING network (`p2_01`)
- **Network:** 131 edges (score ≥ 0.4). 39 have experimental or curated-database evidence. Density 0.25, and NECAB1 is isolated.
- **Top-5 hubs: SNAP25, GFAP, NEFL, SST, CAMK2A.** SNAP25, NEFL, GFAP and CAMK2A remain hubs with experimental/database edges only. The script stops if VGF ever enters the top 5, which would indicate a wrong network input.
- **Modularity:** the 3-module partition gives Q = 0.23 vs a degree-preserving null of −0.02 ± 0.03 (z = 8.6, p < 0.001). Clustering is 0.62 vs 0.39.
- **Robustness:** removing the highest-betweenness genes fragments the network after 7 removals; random removal still leaves 75% connected at that point.

## 1. Differential co-expression (`p2_02`)

| | Microarray (97/98) raw | Microarray, cell-type-adjusted | RNA-seq (12/10) | SEA-AD (10/10) |
|---|---|---|---|---|
| Gene pairs changing (FDR < 0.05) | **108 / 465** | **3 / 465** | 0 / 528 | 32 / 528 |
| Mean \|r\| ctrl → AD (perm. p) | 0.50 → 0.51 (0.55) | 0.22 → 0.21 (0.37) | 0.49 → 0.43 (0.48) | 0.39 → 0.54 (0.32) |

- **Cell-type proxies** come from marker genes outside the panel. Neuron proxy: −0.76 SD in AD (p < 10⁻⁴). SST-interneuron proxy: −0.46 SD (p < 10⁻⁴).
- **Raw microarray, per gene:**
  - SNAP25 gains connectivity (Δk +0.10, FDR 0.02), as do NECAB1 and VSNL1.
  - VIM decouples (Δk −0.26, FDR 0.02).
  - The VIM–GFAP edge weakens (r 0.56 → 0.16).
- **After adjustment** no gene's Δk survives FDR. The ranking does not replicate in RNA-seq (ρ = 0.36) or SEA-AD (ρ = −0.18).
- **Conclusion:** the apparent network rewiring is a cell-composition effect (neuron and SST-interneuron loss), not cell-intrinsic rewiring.

## 2. Module coherence and WGCNA (`p2_05`, `p2_04`)
- **Coherence** (vs 2,000 expression-matched random gene sets):
  - Neuropeptide and Synaptic modules are coherent in all three datasets (p ≤ 0.001; microarray mean r 0.74 and 0.52 vs null ≈ 0.01).
  - The Glial module is coherent in SEA-AD (r 0.73, p = 0.001), weakly in the microarray (0.36, p = 0.03), and not in RNA-seq (p = 0.13).
- **Separability:** Synaptic and Neuropeptide are not separable (between-module r = 0.60 in the microarray; 16 of 19 synaptic genes sit closer to the neuropeptide genes). Glial and neuronal genes are anti-correlated (r ≈ −0.35).
- **WGCNA** (signed; power 16, scale-free R² = 0.82; 13 modules):
  - Neuronal panel genes fall into three modules that go **down** in AD: turquoise (r = −0.41), yellow (−0.40) and red (−0.38), all FDR < 10⁻⁷.
  - Glial genes fall into brown (vascular/astrocyte, **up**, r = +0.51, FDR 6 × 10⁻¹³), green (astro-vascular, not associated) and magenta (microglia, up, FDR 0.04).
- **Answer:** the data support a neuronal axis vs a glial axis, not three independent modules.

## 3. TF enrichment, TRRUST v2 (`p2_03`)

| Gene set | Background | Top TF (overlap) | FDR | Other TFs at FDR < 0.05 |
|---|---|---|---|---|
| 33-gene panel | TRRUST targets on array | **REST** (VGF, CHGB, TAC1, CARTPT) | 3.1 × 10⁻⁴ | FOS, JUN, ETV4 |
| 33-gene panel | all microarray genes | **REST** | 4 × 10⁻⁶ | about 20 TFs (background inflation) |
| Concordant DOWN (79 genes) | TRRUST targets on array | **REST** (VGF, CHGB, BDNF, PCSK1, TUBB3) | 2 × 10⁻⁵ | NFATC1 |
| Concordant UP (105 genes) | either | none | — | — |

**REST log2FC (AD − control):**
- Microarray: +0.149 (padj 2 × 10⁻⁶).
- RNA-seq: +0.155 (padj 0.32).
- SEA-AD: +0.96 unshrunk (+1.66 apeglm; padj 0.068).

REST is higher in AD in all three datasets but significant only in the microarray. Treat it as a hypothesis: Lu et al. (*Nature* 2014) report reduced *nuclear* REST protein in AD.

## 4. Cell-type resolution (`p2_06`)
- **Specificity:** panel genes are cell-type specific (τ up to 0.99). SST, NPY and PENK peak in Sst neurons, TAC1 in Pvalb neurons, GFAP and GJA1 in astrocytes, TREM2 in microglia and VIM in vascular cells.
- **Composition** (dementia vs no dementia; relative proportions, because SEA-AD sorts nuclei at fixed ratios):
  - Sst neurons ↓ 2× (p = 0.002), Pvalb ↓ (p = 0.03).
  - Astrocytes ↑ 2.8× (p = 0.03), vascular cells ↑ (p = 0.009).
- **Whole tissue vs within the expressing cell type** (unshrunk log2FC):
  - **Glial genes:** median **+1.68** in tissue vs **+0.09** within astrocytes or microglia. GFAP +1.33 → +0.03, VIM +2.04 → −0.04, GJA1 +2.73 → +0.43, TREM2 +1.22 → +0.14.
  - **SST, NPY, TAC1:** median −1.09 in tissue vs −0.35 within their cell types. SST itself is −1.46 in tissue vs −0.35 inside Sst neurons.
  - **CARTPT, NPTX2, CHGB:** fall by about as much inside excitatory neurons as in tissue (−1.14, −0.74, −0.57), but not significantly at n = 10 vs 10. These are the candidate cell-intrinsic changes.

## 5. Figures (`figures/`, PDF)

| File | Content |
|---|---|
| `fig1_network.pdf` | STRING network grouped by module; size = degree; colour = consensus log2FC; solid edges = experimental/database evidence |
| `fig2_differential_network.pdf` | Δr (AD − control) on STRING edges, raw vs composition-adjusted |
| `fig3_modules.pdf` | Coherence vs null, within/between-module r, WGCNA mapping |
| `fig4_celltype.pdf` | Cell-type specificity, composition shift, whole tissue vs within cell type |
| `fig5_tf_network.pdf` | TRRUST regulators of the panel (TRRUST-target background) |
| `fig6_robustness.pdf` | Attack tolerance; hub degree under stricter evidence |
| `wgcna_dendrogram.pdf` | WGCNA dendrogram |

## 6. Limitations
- **Curated panel:** the 33 genes were chosen for direction consistency across the same three datasets, so their coherence and consistency are partly circular.
- **Literature bias in STRING:** STRING includes text-mining evidence that favours well-studied genes. The experimental/database-only check partly addresses this.
- **Small cohorts:** only the microarray (n = 195) has the power to detect differences in co-expression.
- **Approximate composition adjustment:** bulk data cannot fully separate composition from regulation, and the marker-proxy adjustment is approximate.
- **SEA-AD labels:** SEA-AD groups are defined by cognitive status, and 5 of 10 controls have intermediate or high AD pathology.
