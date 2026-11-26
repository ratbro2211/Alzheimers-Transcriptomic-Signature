"""
p2_01_string_network.py
Topology of the STRING v12.5 functional-association network of the 33 genes
(downloaded via the STRING API into data/, combined score >= 0.4).

NOTE ON TERMINOLOGY: STRING edges are undirected functional associations (co-expression,
experiments, curated databases, text mining). This is a protein-protein / functional
association network, NOT a gene regulatory network. Regulatory (TF -> target) edges are
added separately in p2_03 (TRRUST).

Outputs: results/string_edges.csv, results/hub_genes_string.csv,
         results/network_null_model.csv, results/attack_tolerance.csv
Run from network/
"""
import numpy as np, pandas as pd, networkx as nx
from networkx.algorithms.community import louvain_communities, modularity

rng = np.random.default_rng(1)
G33 = ("VGF SST SCG2 CARTPT CHGB VSNL1 TAC1 SNAP25 NPTX2 NEFL GABRA1 NRGN GABRG2 SYT4 CCK NEFM SHANK2 "
       "GRIN2A NECAB1 LRRTM1 CBLN2 CAMK2A SYN2 HOMER1 NPY SCN1A PENK NLGN1 CADM1 VIM TREM2 GJA1 GFAP").split()
MOD = {g: "Glial" for g in ["GFAP", "VIM", "GJA1", "TREM2"]}
MOD.update({g: "Neuropeptide" for g in ["VGF", "SST", "NPY", "CARTPT", "TAC1", "PENK", "CHGB", "SCG2", "VSNL1"]})
MOD.update({g: "Synaptic" for g in G33 if g not in MOD})

e = pd.read_csv("data/string_network_33genes_score400.tsv", sep="\t")
e = e.rename(columns={"preferredName_A": "A", "preferredName_B": "B"})
e = e.drop_duplicates(subset=["A", "B"])
e["exp_or_db"] = (e.escore > 0) | (e.dscore > 0)
e["textmining_only"] = (e[["nscore", "fscore", "pscore", "ascore", "escore", "dscore"]].sum(1) == 0)
e[["A", "B", "score", "escore", "dscore", "tscore", "ascore", "exp_or_db", "textmining_only"]].to_csv("results/string_edges.csv", index=False)
print(f"STRING edges (score>=0.4): {len(e)};  with experimental/database evidence: {e.exp_or_db.sum()};  text-mining only: {e.textmining_only.sum()}")

def build(edges):
    g = nx.Graph(); g.add_nodes_from(G33)
    for r in edges.itertuples(): g.add_edge(r.A, r.B, weight=r.score)
    return g

def participation(g, part):
    """Guimera-Amaral participation coefficient: 0 = all links inside own module, ->1 = spread across modules."""
    pc = {}
    for n in g:
        k = g.degree(n)
        if k == 0: pc[n] = 0.0; continue
        cnt = pd.Series([part[m] for m in g.neighbors(n)]).value_counts()
        pc[n] = 1 - ((cnt / k) ** 2).sum()
    return pc

def centralities(g, tag):
    bc = nx.betweenness_centrality(g, normalized=True)
    wbc = nx.betweenness_centrality(g, weight=lambda u, v, d: 1 - d["weight"] + 1e-3, normalized=True)
    cc = nx.closeness_centrality(g)
    giant = g.subgraph(max(nx.connected_components(g), key=len))
    ev = nx.eigenvector_centrality_numpy(giant, weight="weight"); ev = {n: ev.get(n, 0) for n in g}
    pc = participation(g, MOD)
    df = pd.DataFrame({"module": pd.Series(MOD), "degree": dict(g.degree()),
                       "strength": dict(g.degree(weight="weight")), "betweenness": bc,
                       "betweenness_weighted": wbc, "closeness": cc, "eigenvector": ev,
                       "participation_coef": pc}).loc[G33]
    # composite hub score = mean rank across degree, betweenness, eigenvector
    df["hub_rank"] = df[["degree", "betweenness", "eigenvector"]].rank(ascending=False).mean(1)
    return df.add_suffix(f"_{tag}").rename(columns={f"module_{tag}": "module"})

g_all = build(e)
g_exp = build(e[e.exp_or_db])
g_hc = build(e[e.score >= 0.7])
H = centralities(g_all, "all").join(centralities(g_exp, "expdb").drop(columns="module")) \
                              .join(centralities(g_hc, "hc700").drop(columns="module"))
H = H.sort_values("hub_rank_all")
top5 = list(H.index[:5])
assert "VGF" not in top5, f"STOP: VGF in top-5 hubs {top5} - network input is wrong (hand-typed edges?)"
print("Top-5 hubs:", top5, "(sanity check passed: VGF not in top 5)")
H.round(4).to_csv("results/hub_genes_string.csv")
cols = ["module", "degree_all", "betweenness_all", "eigenvector_all", "participation_coef_all",
        "degree_expdb", "betweenness_expdb", "degree_hc700"]
print("\nTop 12 hubs (full STRING >=0.4), with robustness columns:")
print(H[cols].head(12).round(3).to_string())
iso = [n for n in G33 if g_all.degree(n) == 0]
print(f"\nIsolated genes (no STRING edge >=0.4): {iso}")
print(f"Density {nx.density(g_all):.3f}; giant component {len(max(nx.connected_components(g_all), key=len))}/33; "
      f"avg clustering {nx.average_clustering(g_all):.3f}")

# ---- Modularity of the predefined 3-module partition + Louvain, vs degree-preserving nulls ----
gc = g_all.subgraph(max(nx.connected_components(g_all), key=len)).copy()
pre = [{n for n in gc if MOD[n] == m} for m in ["Glial", "Neuropeptide", "Synaptic"]]
pre = [s for s in pre if s]
Q_pre = modularity(gc, pre)
louv = louvain_communities(gc, seed=1, resolution=1.0)
Q_louv = modularity(gc, louv)
null_pre, null_louv, null_clust = [], [], []
for i in range(1000):
    r = gc.copy()
    nx.double_edge_swap(r, nswap=4 * r.number_of_edges(), max_tries=10 ** 5, seed=int(rng.integers(1e9)))
    null_pre.append(modularity(r, pre))
    null_louv.append(modularity(r, louvain_communities(r, seed=1)))
    null_clust.append(nx.average_clustering(r))
def z(x, nl): return (x - np.mean(nl)) / np.std(nl), (np.sum(np.array(nl) >= x) + 1) / (len(nl) + 1)
rows = []
for name, val, nl in [("Q predefined 3 modules", Q_pre, null_pre), ("Q Louvain", Q_louv, null_louv),
                      ("avg clustering", nx.average_clustering(gc), null_clust)]:
    zz, p = z(val, nl); rows.append(dict(statistic=name, observed=round(val, 3), null_mean=round(np.mean(nl), 3),
                                         null_sd=round(np.std(nl), 3), z=round(zz, 2), p_empirical=p))
NM = pd.DataFrame(rows); NM.to_csv("results/network_null_model.csv", index=False)
print("\nDegree-preserving null model (1000 rewirings, giant component):"); print(NM.to_string(index=False))
louv_ordered = sorted(louv, key=lambda c: (-len(c), min(c)))   # deterministic numbering (size, then first gene)
lv = pd.Series({n: i for i, c in enumerate(louv_ordered) for n in c}, name="louvain_community")
lv = lv.loc[sorted(lv.index, key=lambda g: (lv[g], g))]            # deterministic row order
lv.to_csv("results/louvain_communities.csv")
print("\nLouvain communities:")
for i, c in enumerate(louv_ordered):
    print(f"  C{i}: {sorted(c)}")
from sklearn.metrics import adjusted_rand_score
common = list(lv.index)
print(f"Adjusted Rand index Louvain vs predefined modules: {adjusted_rand_score([MOD[n] for n in common], lv[common]):.2f}")

# ---- Robustness: random failure vs targeted hub attack (complex-systems classic) ----
def attack(g, order):
    g = g.copy(); out = [len(max(nx.connected_components(g), key=len)) / 33]
    for n in order[:-1]:
        g.remove_node(n); out.append(len(max(nx.connected_components(g), key=len)) / 33 if len(g) else 0)
    return out
deg_order = list(H.index)                                   # by composite hub rank
bt = []                                                      # adaptive betweenness attack
g = g_all.copy()
while len(g):
    b = nx.betweenness_centrality(g); n = max(b, key=b.get); bt.append(n); g.remove_node(n)
rand = np.mean([attack(g_all, list(rng.permutation(G33))) for _ in range(500)], axis=0)
A = pd.DataFrame({"n_removed": range(33), "random_failure": rand, "hub_attack": attack(g_all, deg_order),
                  "adaptive_betweenness_attack": attack(g_all, bt)})
A.to_csv("results/attack_tolerance.csv", index=False)
k = int(np.argmax(A.adaptive_betweenness_attack.values <= 0.5))
print(f"\nGiant component falls below 50% after removing {k} genes by adaptive betweenness attack "
      f"({bt[:k]}), vs random failure: {A.random_failure.iloc[k]:.2f} of nodes still connected")
nx.write_graphml(g_all, "results/string_network.graphml")
