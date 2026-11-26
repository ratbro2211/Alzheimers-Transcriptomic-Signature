"""
p2_05_module_coherence.py
Are the three hand-defined modules (Glial / Neuropeptide / Synaptic) truly co-expressed?

For each dataset and each module:
  - mean pairwise Spearman r within the module (all samples, and AD / control separately)
  - empirical p-value vs 2000 random gene sets of the same size, matched to the module
    genes' mean-expression distribution (expression-matched null: highly expressed genes
    correlate more, so an unmatched null would be too lenient)
  - silhouette of the 3-module partition on the 1 - r distance (does each gene sit closer
    to its own module than to the others?)
Also reports the between-module correlation (Glial vs neuronal anticorrelation is the
expected signature of cell-composition variation).
Run from network/
"""
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import silhouette_samples

rng = np.random.default_rng(3)
NNULL = 2000
G33 = ("VGF SST SCG2 CARTPT CHGB VSNL1 TAC1 SNAP25 NPTX2 NEFL GABRA1 NRGN GABRG2 SYT4 CCK NEFM SHANK2 "
       "GRIN2A NECAB1 LRRTM1 CBLN2 CAMK2A SYN2 HOMER1 NPY SCN1A PENK NLGN1 CADM1 VIM TREM2 GJA1 GFAP").split()
MOD = {g: "Glial" for g in ["GFAP", "VIM", "GJA1", "TREM2"]}
MOD.update({g: "Neuropeptide" for g in ["VGF", "SST", "NPY", "CARTPT", "TAC1", "PENK", "CHGB", "SCG2", "VSNL1"]})
MOD.update({g: "Synaptic" for g in G33 if g not in MOD})

def spearman(M):  # genes x samples
    R = np.apply_along_axis(stats.rankdata, 1, M); return np.corrcoef(R)

def mean_offdiag(C):
    iu = np.triu_indices_from(C, 1); return np.nanmean(C[iu])

rows, sil_rows, between = [], [], []
for name in ["ARR", "RNA", "SEA"]:
    X = pd.read_csv(f"data/expr_{name}.csv.gz", index_col=0)
    L = pd.read_csv(f"data/labels_{name}.csv", index_col=0).group.reindex(X.columns)
    X = X[X.var(1) > 0]
    mu = X.mean(1); bins = pd.qcut(mu, 20, labels=False, duplicates="drop")
    pool = {b: X.index[bins == b].to_numpy() for b in np.unique(bins)}
    for m in ["Glial", "Neuropeptide", "Synaptic"]:
        genes = [g for g in G33 if MOD[g] == m and g in X.index]
        gb = bins[genes].values
        for subset, cols in [("all", X.columns), ("ctrl", X.columns[L == "ND"]), ("AD", X.columns[L == "AD"])]:
            obs = mean_offdiag(spearman(X.loc[genes, cols].values))
            Xs = X[cols]
            null = []
            for _ in range(NNULL):
                rs = [rng.choice(pool[b]) for b in gb]
                null.append(mean_offdiag(spearman(Xs.loc[rs].values)))
            null = np.array(null)
            rows.append(dict(dataset=name, module=m, subset=subset, n_genes=len(genes), n_samples=len(cols),
                             mean_r=round(obs, 3), null_mean=round(null.mean(), 3), null_95pct=round(np.quantile(null, .95), 3),
                             z=round((obs - null.mean()) / null.std(), 2),
                             p_empirical=(np.sum(null >= obs) + 1) / (NNULL + 1)))
    # silhouette + between-module correlation (all samples)
    genes = [g for g in G33 if g in X.index]
    C = spearman(X.loc[genes].values)
    lab = np.array([MOD[g] for g in genes])
    s = silhouette_samples(1 - C, lab, metric="precomputed")
    for g, v, l in zip(genes, s, lab): sil_rows.append(dict(dataset=name, gene=g, module=l, silhouette=round(v, 3)))
    idx = {m: [i for i, l in enumerate(lab) if l == m] for m in ["Glial", "Neuropeptide", "Synaptic"]}
    for a in idx:
        for b in idx:
            between.append(dict(dataset=name, module_A=a, module_B=b, mean_r=round(np.nanmean(C[np.ix_(idx[a], idx[b])][~np.eye(len(idx[a]), dtype=bool)] if a == b else C[np.ix_(idx[a], idx[b])]), 3)))
    print(f"{name} done")

R = pd.DataFrame(rows); R.to_csv("results/module_coherence.csv", index=False)
S = pd.DataFrame(sil_rows); S.to_csv("results/module_silhouette.csv", index=False)
B = pd.DataFrame(between).pivot_table(index=["dataset", "module_A"], columns="module_B", values="mean_r")
B.to_csv("results/module_between_correlation.csv")
pd.set_option("display.width", 200)
print("\nWithin-module coherence vs expression-matched random gene sets:")
print(R.to_string(index=False))
print("\nMean silhouette per module (>0: genes closer to own module):")
print(S.groupby(["dataset", "module"]).silhouette.mean().round(3).unstack())
print("\nGenes with negative silhouette (sit closer to another module), microarray:")
print(S[(S.dataset == "ARR") & (S.silhouette < 0)].to_string(index=False))
print("\nMean correlation between/within modules:"); print(B)
