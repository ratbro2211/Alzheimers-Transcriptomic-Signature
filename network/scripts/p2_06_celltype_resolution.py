"""
p2_06_celltype_resolution.py
Cell-type resolution of the 33-gene panel from the analysis/ outputs.
No DESeq2 is rerun: within-cell-type results come from
analysis/dataset3_seaad/results/celltype/panel33_by_celltype.csv.

  1. Specificity: mean log2 CPM per cell type (10 SEA-AD classes) from the donor x
     cell-type raw-UMI pseudobulk (script 03c); tau specificity index; top cell type.
  2. Composition: per-donor proportions of each class, dementia vs no dementia.
  3. Bulk vs within-cell-type log2FC for every panel gene (unshrunk MLE for both,
     so the comparison is like with like), from panel33_by_celltype.csv.
Run from network/.
"""
import subprocess, numpy as np, pandas as pd
from scipy import stats

F = "../analysis/dataset3_seaad/results"
P33 = pd.read_csv("../analysis/comparison/results/panel_33gene_summary.csv").gene.tolist()
MOD = {g: "Glial" for g in ["GFAP", "VIM", "GJA1", "TREM2"]}
MOD.update({g: "Neuropeptide" for g in ["VGF", "SST", "NPY", "CARTPT", "TAC1", "PENK", "CHGB", "SCG2", "VSNL1"]})
MOD.update({g: "Synaptic" for g in P33 if g not in MOD})

# 1. specificity: export per-cell-type CPM of the panel genes from the 03c rds via R (rds is R-only)
subprocess.run(["Rscript", "-e",
    f'pb <- readRDS("{F}/celltype/seaad_celltype_pseudobulk.rds"); g <- c({",".join(chr(34)+x+chr(34) for x in P33)});'
    'lib <- colSums(pb$counts); cpm <- sweep(pb$counts[g, ], 2, lib, "/") * 1e6;'
    'out <- data.frame(sample = colnames(cpm), celltype = pb$coldata$celltype, n_nuclei = pb$coldata$n_nuclei, t(cpm), check.names = FALSE);'
    'write.csv(out, "results/celltype_cpm_33genes_by_sample.csv", row.names = FALSE)'], check=True)
cp = pd.read_csv("results/celltype_cpm_33genes_by_sample.csv")
cp = cp[cp.n_nuclei >= 20]
M = cp.groupby("celltype")[P33].mean().T                    # genes x cell types, mean CPM
M.round(3).to_csv("results/celltype_mean_cpm_33genes.csv")
tau = (1 - M.div(M.max(1), axis=0)).sum(1) / (M.shape[1] - 1)
spec = pd.DataFrame({"module": pd.Series(MOD), "top_celltype": M.idxmax(1),
                     "second_celltype": M.apply(lambda r: r.nlargest(2).index[1], axis=1),
                     "top_mean_CPM": M.max(1).round(1), "tau_specificity": tau.round(3)}).loc[P33]
spec.to_csv("results/celltype_specificity_33genes.csv")
print("Top cell type per panel gene (mean CPM, donor x cell-type pseudobulk):")
print(spec.sort_values(["module", "top_celltype"]).to_string())

# 2. composition (10 classes)
G = {"L2/3 IT": "Exc_IT", "L4 IT": "Exc_IT", "L5 IT": "Exc_IT", "L6 IT": "Exc_IT", "L6 IT Car3": "Exc_IT", "L5 ET": "Exc_deep",
     "L5/6 NP": "Exc_deep", "L6 CT": "Exc_deep", "L6b": "Exc_deep", "Pvalb": "Inh_Pvalb", "Chandelier": "Inh_Pvalb", "Sst": "Inh_Sst",
     "Sst Chodl": "Inh_Sst", "Vip": "Inh_CGE", "Lamp5": "Inh_CGE", "Lamp5 Lhx6": "Inh_CGE", "Sncg": "Inh_CGE", "Pax6": "Inh_CGE",
     "Astrocyte": "Astrocyte", "Immune": "Microglia", "Oligodendrocyte": "Oligodendrocyte", "OPC": "OPC",
     "Endothelial": "Vascular", "VLMC & Perivascular": "Vascular"}
C = pd.read_csv(f"{F}/seaad_celltype_counts_per_donor.csv", index_col=0).T.groupby(G).sum().T
meta = pd.read_csv(f"{F}/seaad_pseudobulk_coldata.csv").set_index("donor_id").loc[C.index]
ad = (meta.condition == "AD").values
Pp = C.div(C.sum(1), axis=0)
comp = pd.DataFrame({"mean_prop_ctrl": Pp[~ad].mean(), "mean_prop_dementia": Pp[ad].mean()})
comp["log2_ratio"] = np.log2(comp.mean_prop_dementia / comp.mean_prop_ctrl)
comp["p_mannwhitney"] = [stats.mannwhitneyu(Pp.loc[ad, c], Pp.loc[~ad, c]).pvalue for c in Pp]
comp.sort_values("p_mannwhitney").round(5).to_csv("results/celltype_composition_dementia.csv")
print("\nComposition, dementia vs no dementia (SEA-AD sorting ratios fixed by design -> relative proportions):")
print(comp.sort_values("p_mannwhitney").round(4).to_string())

# 3. bulk vs within-cell-type (from the cell-type DESeq2 in analysis/; not rerun)
W = pd.read_csv(f"{F}/celltype/panel33_by_celltype.csv")
W["module"] = W.gene.map(MOD)
W["fraction_of_bulk_change_within_celltype"] = (W.home_log2FC_mle / W.whole_tissue_log2FC_mle).round(2)
W.to_csv("results/celltype_within_celltype_DE.csv", index=False)
print("\nWhole-tissue vs within-home-cell-type log2FC (unshrunk), panel genes:")
print(W[["gene", "module", "home_celltype", "whole_tissue_log2FC_mle", "home_log2FC_mle", "home_padj"]].round(3).to_string(index=False))
gl = W[W.module == "Glial"]; npg = W[W.gene.isin(["SST", "NPY", "TAC1"])]
print(f"\nGlial genes: median whole-tissue {gl.whole_tissue_log2FC_mle.median():+.2f} vs within astrocyte/microglia {gl.home_log2FC_mle.median():+.2f}")
print(f"Interneuron neuropeptides (SST, NPY, TAC1): whole tissue {npg.whole_tissue_log2FC_mle.median():+.2f} vs within cell type {npg.home_log2FC_mle.median():+.2f}")
