"""
p2_02_differential_network.py
"Dynamic" (condition-specific / differential) network analysis of the 33 genes.

A static STRING network cannot change between AD and control. What CAN change is the
co-expression (coordination) of the genes across individuals. For each dataset we
build a control network and an AD network from Spearman correlations and ask:

  1. Global coherence: does the mean |r| among the 33 genes change in AD?
  2. Per-gene differential connectivity  dk_i = k_i(AD) - k_i(ctrl),
     k_i = mean |r_ij| over the other genes; p-value from 2000 label permutations.
  3. Differential edges: Fisher-z test of r_AD vs r_ctrl for every pair (BH-FDR);
     "rewiring degree" = number of differential edges per gene.
  4. STRING-constrained networks: STRING edges weighted by co-expression in each
     condition; change in weighted strength of each hub.
  5. Cell-composition confound: bulk co-expression of neuronal genes is partly
     driven by how many neurons are in each sample (neuron loss in AD). We repeat
     1-2 in the microarray after regressing out marker-based proxies for neuron,
     astrocyte, microglia, oligodendrocyte and SST-interneuron abundance
     (proxies built from marker genes NOT in the 33-gene set).

Primary dataset: microarray (97 AD / 98 ND). RNA-seq (12/10) and SEA-AD (10/10)
are underpowered for correlation differences and are reported as replication only.
Run from network/
"""
import numpy as np, pandas as pd, networkx as nx
from scipy import stats
from statsmodels.stats.multitest import multipletests

rng = np.random.default_rng(7)
NPERM = 2000
G33 = ("VGF SST SCG2 CARTPT CHGB VSNL1 TAC1 SNAP25 NPTX2 NEFL GABRA1 NRGN GABRG2 SYT4 CCK NEFM SHANK2 "
       "GRIN2A NECAB1 LRRTM1 CBLN2 CAMK2A SYN2 HOMER1 NPY SCN1A PENK NLGN1 CADM1 VIM TREM2 GJA1 GFAP").split()
MARKERS = {  # cell-type abundance proxies (none of these is in the 33-gene set)
    "neuron": ["RBFOX3", "STMN2", "SYT1", "GAP43", "SNAP91", "ELAVL4", "INA", "RTN1"],
    "astrocyte": ["AQP4", "ALDH1L1", "SLC1A3", "SLC1A2", "GLUL", "SOX9", "ALDOC"],
    "microglia": ["CX3CR1", "P2RY12", "CSF1R", "C1QB", "C1QA", "TYROBP", "AIF1"],
    "oligodendrocyte": ["MBP", "PLP1", "MOG", "MOBP", "MAG", "CNP", "ST18"],
    "sst_interneuron": ["CORT", "NOS1", "CHODL"],
}

def load(name):
    X = pd.read_csv(f"data/expr_{name}.csv.gz", index_col=0)
    L = pd.read_csv(f"data/labels_{name}.csv", index_col=0).group.reindex(X.columns)
    return X, L

def rank_corr(M):                      # M: samples x genes -> Spearman matrix
    R = np.apply_along_axis(stats.rankdata, 0, M)
    return np.corrcoef(R, rowvar=False)

def conn(C):                           # mean |r| to other genes
    A = np.abs(C.copy()); np.fill_diagonal(A, np.nan); return np.nanmean(A, 1)

def analyse(X, L, genes, tag):
    genes = [g for g in genes if g in X.index]
    M = X.loc[genes].T.values; y = (L == "AD").values
    Cad, Cnd = rank_corr(M[y]), rank_corr(M[~y])
    dk = conn(Cad) - conn(Cnd)
    glob = lambda Ca, Cn: np.nanmean(np.abs(Ca[np.triu_indices_from(Ca, 1)])) - np.nanmean(np.abs(Cn[np.triu_indices_from(Cn, 1)]))
    g_obs = glob(Cad, Cnd)
    null_dk, null_g = [], []
    for _ in range(NPERM):
        yp = rng.permutation(y); a, n = rank_corr(M[yp]), rank_corr(M[~yp])
        null_dk.append(conn(a) - conn(n)); null_g.append(glob(a, n))
    null_dk = np.array(null_dk)
    p_dk = (np.sum(np.abs(null_dk) >= np.abs(dk), 0) + 1) / (NPERM + 1)
    p_g = (np.sum(np.abs(null_g) >= abs(g_obs)) + 1) / (NPERM + 1)
    # differential edges (Fisher z)
    na, nn = y.sum(), (~y).sum()
    iu = np.triu_indices(len(genes), 1)
    za, zn = np.arctanh(np.clip(Cad[iu], -.999, .999)), np.arctanh(np.clip(Cnd[iu], -.999, .999))
    zstat = (za - zn) / np.sqrt(1 / (na - 3) + 1 / (nn - 3))
    pe = 2 * stats.norm.sf(np.abs(zstat)); fe = multipletests(pe, method="fdr_bh")[1]
    E = pd.DataFrame({"A": np.array(genes)[iu[0]], "B": np.array(genes)[iu[1]], "r_ctrl": Cnd[iu], "r_AD": Cad[iu],
                      "delta_r": Cad[iu] - Cnd[iu], "z": zstat, "p": pe, "fdr": fe})
    rew = pd.concat([E[E.fdr < .05].A, E[E.fdr < .05].B]).value_counts().reindex(genes).fillna(0).astype(int)
    G = pd.DataFrame({f"k_ctrl_{tag}": conn(Cnd), f"k_AD_{tag}": conn(Cad), f"dk_{tag}": dk,
                      f"p_dk_{tag}": p_dk, f"rewired_edges_{tag}": rew.values}, index=genes)
    G[f"fdr_dk_{tag}"] = multipletests(p_dk, method="fdr_bh")[1]
    mean_abs = lambda C: np.nanmean(np.abs(C[iu]))
    glob_row = dict(dataset=tag, n_AD=int(na), n_ctrl=int(nn), n_genes=len(genes),
                    mean_abs_r_ctrl=round(mean_abs(Cnd), 3), mean_abs_r_AD=round(mean_abs(Cad), 3),
                    delta=round(g_obs, 3), p_perm=p_g, n_diff_edges_fdr05=int((E.fdr < .05).sum()), n_pairs=len(E))
    return G, E, glob_row, (Cnd, Cad, genes)

def residualise(X, L, genes, proxies):
    """Regress each gene on cell-type proxies + group (pooled), keep residual + group effect removed
    only for proxies, i.e. adjusted expression = y - B_proxy * proxy."""
    Z = proxies.values; Zc = Z - Z.mean(0)
    D = np.column_stack([np.ones(len(L)), (L == "AD").astype(float).values, Zc])
    out = {}
    for g in genes:
        if g not in X.index: continue
        b, *_ = np.linalg.lstsq(D, X.loc[g].values, rcond=None)
        out[g] = X.loc[g].values - Zc @ b[2:]
    return pd.DataFrame(out, index=X.columns).T

if __name__ == "__main__":
    allG, allE, glob = [], {}, []
    mats = {}
    for name in ["ARR", "RNA", "SEA"]:
        X, L = load(name)
        G, E, gr, m = analyse(X, L, G33, name); allG.append(G); allE[name] = E; glob.append(gr); mats[name] = m
        print(f"{name}: done ({gr['n_genes']} genes)")

    # composition-adjusted microarray analysis
    X, L = load("ARR")
    prox = {}
    for ct, mk in MARKERS.items():
        mk = [g for g in mk if g in X.index]
        if len(mk) < 2: continue
        Zs = X.loc[mk].sub(X.loc[mk].mean(1), axis=0).div(X.loc[mk].std(1), axis=0)
        prox[ct] = Zs.mean(0)
        print(f"  proxy {ct}: {mk}")
    prox = pd.DataFrame(prox)
    pr = pd.DataFrame({c: stats.mannwhitneyu(prox[c][L == "AD"], prox[c][L == "ND"]).pvalue for c in prox}, index=["p_AD_vs_ND"]).T
    pr["mean_AD_minus_ND"] = [prox[c][L == "AD"].mean() - prox[c][L == "ND"].mean() for c in prox]
    print("Cell-type proxy shifts in AD (microarray):\n", pr.round(4))
    pr.to_csv("results/celltype_proxy_shift_microarray.csv")
    Xadj = residualise(X, L, G33, prox)
    G, E, gr, m = analyse(Xadj, L, G33, "ARRadj"); allG.append(G); allE["ARRadj"] = E; glob.append(gr); mats["ARRadj"] = m

    GL = pd.DataFrame(glob); GL.to_csv("results/diffnet_global_coherence.csv", index=False)
    print("\nGlobal co-expression coherence of the 33 genes (mean |Spearman r|):"); print(GL.to_string(index=False))

    T = pd.concat(allG, axis=1)
    # STRING-constrained condition-specific strength (microarray, raw and adjusted)
    se = pd.read_csv("results/string_edges.csv")
    for tag in ["ARR", "ARRadj"]:
        Cn, Ca, genes = mats[tag]; ix = {g: i for i, g in enumerate(genes)}
        for cond, C in [("ctrl", Cn), ("AD", Ca)]:
            g = nx.Graph(); g.add_nodes_from(genes)
            for r in se.itertuples():
                if r.A in ix and r.B in ix: g.add_edge(r.A, r.B, w=abs(C[ix[r.A], ix[r.B]]))
            T[f"string_strength_{cond}_{tag}"] = pd.Series(dict(g.degree(weight="w")))
        T[f"string_dstrength_{tag}"] = T[f"string_strength_AD_{tag}"] - T[f"string_strength_ctrl_{tag}"]
    H = pd.read_csv("results/hub_genes_string.csv", index_col=0)[["module", "degree_all", "betweenness_all", "hub_rank_all"]]
    T = H.join(T)
    T["replicates_RNA"] = np.sign(T.dk_RNA) == np.sign(T.dk_ARR)
    T["replicates_SEA"] = np.sign(T.dk_SEA) == np.sign(T.dk_ARR)
    T = T.sort_values("dk_ARR")
    T.round(4).to_csv("results/differential_connectivity_33genes.csv")
    for k, E in allE.items(): E.sort_values("p").round(5).to_csv(f"results/differential_edges_{k}.csv", index=False)

    pd.set_option("display.width", 250)
    show = ["module", "degree_all", "k_ctrl_ARR", "k_AD_ARR", "dk_ARR", "fdr_dk_ARR", "rewired_edges_ARR",
            "dk_ARRadj", "fdr_dk_ARRadj", "dk_RNA", "dk_SEA", "string_dstrength_ARR"]
    print("\nPer-gene differential connectivity (sorted by dk in microarray):")
    print(T[show].round(3).to_string())
    print("\nSign agreement of dk with microarray: RNA-seq", T.replicates_RNA.sum(), "/", T.dk_RNA.notna().sum(),
          "| SEA-AD", T.replicates_SEA.sum(), "/", T.dk_SEA.notna().sum())
    print("Spearman(dk_ARR, dk_RNA) =", round(stats.spearmanr(T.dk_ARR, T.dk_RNA, nan_policy="omit")[0], 2),
          "| Spearman(dk_ARR, dk_SEA) =", round(stats.spearmanr(T.dk_ARR, T.dk_SEA, nan_policy="omit")[0], 2))
    print("\nTop differential STRING-supported edges (microarray, raw):")
    E = allE["ARR"]; s = set(map(frozenset, zip(se.A, se.B)))
    E["in_STRING"] = [frozenset((a, b)) in s for a, b in zip(E.A, E.B)]
    print(E[E.in_STRING].sort_values("p").head(10).round(4).to_string(index=False))
    np.savez("results/coexpr_matrices.npz", **{f"{k}_{c}": v for k, (Cn, Ca, g) in mats.items() for c, v in [("ctrl", Cn), ("AD", Ca)]},
             **{f"{k}_genes": np.array(g) for k, (Cn, Ca, g) in mats.items()})
