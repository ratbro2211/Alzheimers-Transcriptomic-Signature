"""
p2_07_figures.py - network figures (PDF)

fig1_network.png          STRING network: module-grouped layout, node size = degree,
                          node colour = consensus AD log2FC, edge style = evidence type
fig2_differential_net.png AD-vs-control co-expression change on STRING edges, raw vs
                          cell-composition-adjusted (microarray, n = 195)
fig3_modules.png          module coherence vs matched null, between/within-module r, WGCNA mapping
fig4_celltype.png         SEA-AD cell-type specificity, composition shift, bulk vs within-cell-type LFC
fig5_tf_network.png       TRRUST regulators of the 33-gene set
fig6_robustness.png       attack tolerance + hub robustness across evidence filters
Run from network/
"""
import numpy as np, pandas as pd, networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull

# ---- palette (validated reference instance; categorical slots 1-3 pass all-pairs CVD) ----
SURF, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
CAT = ["#2a78d6", "#eb6834", "#1baf7a"]
DIV = LinearSegmentedColormap.from_list("div", ["#1f5fb0", "#7fa9e0", "#f0efec", "#ec8f86", "#c42b2b"])
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
                     "font.family": "DejaVu Sans", "font.size": 9, "text.color": INK, "axes.labelcolor": INK2,
                     "axes.edgecolor": MUTED, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.titleweight": "bold",
                     "axes.titlesize": 10, "axes.titlelocation": "left"})
FR = "../analysis"
G33 = ("VGF SST SCG2 CARTPT CHGB VSNL1 TAC1 SNAP25 NPTX2 NEFL GABRA1 NRGN GABRG2 SYT4 CCK NEFM SHANK2 "
       "GRIN2A NECAB1 LRRTM1 CBLN2 CAMK2A SYN2 HOMER1 NPY SCN1A PENK NLGN1 CADM1 VIM TREM2 GJA1 GFAP").split()
MOD = {g: "Glial" for g in ["GFAP", "VIM", "GJA1", "TREM2"]}
MOD.update({g: "Neuropeptide" for g in ["VGF", "SST", "NPY", "CARTPT", "TAC1", "PENK", "CHGB", "SCG2", "VSNL1"]})
MOD.update({g: "Synaptic" for g in G33 if g not in MOD})
MODS = ["Synaptic", "Neuropeptide", "Glial"]
MCOL = dict(zip(MODS, CAT))

E = pd.read_csv("results/string_edges.csv")
H = pd.read_csv("results/hub_genes_string.csv", index_col=0)
panel = pd.read_csv(f"{FR}/comparison/results/panel_33gene_summary.csv", index_col="gene")
lfc = panel.loc[G33, ["log2FC_micro", "log2FC_rnaseq", "log2FC_seaad"]].mean(1)   # SEA-AD apeglm-shrunk

G = nx.Graph(); G.add_nodes_from(G33)
for r in E.itertuples(): G.add_edge(r.A, r.B, w=r.score, exp=r.exp_or_db)

# module-grouped layout: each module on its own circle (nodes ordered by degree), labels placed radially
anchor = {"Synaptic": np.array([1.45, 0.0]), "Neuropeptide": np.array([-1.35, 0.55]), "Glial": np.array([-1.0, -1.55])}
radius = {"Synaptic": 1.05, "Neuropeptide": 0.62, "Glial": 0.38}
deg0 = dict(G.degree())
pos, ang = {}, {}
for m in MODS:
    mem = sorted([g for g in G33 if MOD[g] == m], key=lambda g: -deg0[g])
    # interleave high/low degree so hubs are spread round the circle
    mem = [x for pair in zip(mem[::2], mem[1::2][::-1] + [None]) for x in pair if x][:len(mem)]
    for i, g in enumerate(mem):
        a = np.pi / 2 + 2 * np.pi * i / len(mem); ang[g] = a
        pos[g] = anchor[m] + radius[m] * np.array([np.cos(a), np.sin(a)])

def label(ax, g, size=7.5, weight="normal", color=None, off=0.13):
    a = ang[g]; x, y = pos[g] + off * np.array([np.cos(a), np.sin(a)])
    ha = "left" if np.cos(a) > 0.3 else "right" if np.cos(a) < -0.3 else "center"
    va = "bottom" if np.sin(a) > 0.3 else "top" if np.sin(a) < -0.3 else "center"
    ax.text(x, y, g, ha=ha, va=va, fontsize=size, fontweight=weight, color=color or INK, zorder=4)

def hulls(ax, label=True):
    for m in MODS:
        pts = np.array([pos[g] for g in G33 if MOD[g] == m])
        c = pts.mean(0); pts = c + (pts - c) * 1.0
        ring = np.vstack([pts + d for d in 0.42 * np.array([[1, 0], [-1, 0], [0, 1], [0, -1], [.7, .7], [-.7, .7], [.7, -.7], [-.7, -.7]])])
        hv = ring[ConvexHull(ring).vertices]
        ax.fill(hv[:, 0], hv[:, 1], color=MCOL[m], alpha=0.07, lw=0, zorder=0)
        ax.plot(*np.vstack([hv, hv[:1]]).T, color=MCOL[m], lw=1.2, alpha=.8, zorder=0)
        if label:
            top = hv[np.argmax(hv[:, 1])]
            ax.text(top[0], top[1] + 0.06, f"{m} module", color=INK2, fontsize=9, ha="center", fontweight="bold")

def draw_nodes(ax, color, size, cmap=None, norm=None, edgecolor=INK):
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=color, node_size=size, cmap=cmap,
                           vmin=None if norm else None, edgecolors=edgecolor, linewidths=0.8)

# ============================== FIG 1 =========================================
fig, ax = plt.subplots(figsize=(10.5, 8.2))
hulls(ax)
for u, v, d in G.edges(data=True):
    x, y = zip(pos[u], pos[v])
    ax.plot(x, y, color=INK2 if d["exp"] else MUTED, lw=0.4 + 2.6 * (d["w"] - .4) / .6,
            ls="-" if d["exp"] else (0, (3, 2)), alpha=.75 if d["exp"] else .45, zorder=1)
norm = TwoSlopeNorm(vmin=-1.2, vcenter=0, vmax=1.2)
deg = H.loc[G33, "degree_all"]
sc = ax.scatter([pos[g][0] for g in G33], [pos[g][1] for g in G33], s=90 + 38 * deg.values,
                c=lfc.values, cmap=DIV, norm=norm, edgecolors=INK, linewidths=.8, zorder=3)
top5 = list(H.sort_values("hub_rank_all").index[:5])
for g in G33:
    label(ax, g, size=9 if g in top5 else 7.5, weight="bold" if g in top5 else "normal", off=0.1 + 0.004 * deg[g])
cb = fig.colorbar(sc, ax=ax, shrink=.45, pad=.01); cb.set_label("Consensus AD log2FC\n(mean of 3 datasets)", color=INK2)
cb.outline.set_visible(False)
leg = [Line2D([], [], color=INK2, lw=2, label="experimental / curated-database evidence"),
       Line2D([], [], color=MUTED, lw=1.5, ls=(0, (3, 2)), label="co-expression / text-mining only"),
       Line2D([], [], color=INK2, lw=0.6, label="thin = STRING score 0.4"),
       Line2D([], [], color=INK2, lw=3, label="thick = STRING score 1.0")]
for k in [2, 8, 16]:
    leg.append(Line2D([], [], marker="o", ls="", markerfacecolor="#f0efec", markeredgecolor=INK,
                      markersize=np.sqrt(90 + 38 * k) / 1.0, label=f"degree {k}"))
ax.legend(handles=leg, loc="lower right", frameon=False, fontsize=7.5, labelspacing=1.1, borderpad=0, ncol=2, columnspacing=1.5)
ax.set_title("STRING v12.5 functional-association network of the 33 AD genes  (131 edges, score >= 0.4)", pad=22)
ax.text(0, 1.01, f"Bold = top-5 hubs by composite rank ({', '.join(top5)}).  Node colour: blue = lower, red = higher in AD.",
        transform=ax.transAxes, fontsize=8, color=INK2, va="bottom")
ax.axis("off"); ax.set_aspect("equal")
fig.tight_layout(); fig.savefig("figures/fig1_network.pdf", dpi=220); plt.close()

# ============================== FIG 2 =========================================
Z = np.load("results/coexpr_matrices.npz", allow_pickle=True)
fig, axs = plt.subplots(1, 2, figsize=(15, 6.0))
dnorm = TwoSlopeNorm(vmin=-.5, vcenter=0, vmax=.5)
for ax, tag, title in [(axs[0], "ARR", "Raw bulk co-expression"),
                       (axs[1], "ARRadj", "After removing cell-type-abundance proxies")]:
    genes = list(Z[f"{tag}_genes"]); ix = {g: i for i, g in enumerate(genes)}
    dR = Z[f"{tag}_AD"] - Z[f"{tag}_ctrl"]
    de = pd.read_csv(f"results/differential_edges_{tag}.csv")
    sig = {frozenset((a, b)) for a, b, f in zip(de.A, de.B, de.fdr) if f < .05}
    hulls(ax, label=(tag == "ARR"))
    nsig = 0
    for u, v in G.edges():
        if u not in ix or v not in ix: continue
        d = dR[ix[u], ix[v]]; s = frozenset((u, v)) in sig; nsig += s
        x, y = zip(pos[u], pos[v])
        ax.plot(x, y, color=DIV(dnorm(d)), lw=3.2 if s else 1.0, alpha=1 if s else .6, zorder=2 if s else 1)
    ax.scatter([pos[g][0] for g in G33], [pos[g][1] for g in G33], s=130, c="#f0efec",
               edgecolors=[INK if g in ix else MUTED for g in G33], linewidths=.8, zorder=3)
    for g in G33: label(ax, g, color=INK if g in ix else MUTED, off=0.1)
    ntot = (de.fdr < .05).sum()
    ax.set_title(f"{title}\n{nsig} STRING edges and {ntot}/465 gene pairs change (FDR<0.05)")
    ax.axis("off"); ax.set_aspect("equal")
sm = plt.cm.ScalarMappable(cmap=DIV, norm=dnorm)
cb = fig.colorbar(sm, ax=axs, shrink=.5, pad=.01); cb.set_label("Change in Spearman r (AD - control)", color=INK2); cb.outline.set_visible(False)
fig.suptitle("Differential co-expression on the STRING network, GSE132903 microarray (97 AD / 98 control)\n"
             "Thick edges: FDR<0.05 (Fisher z).  Red = genes more tightly coupled in AD.  Grey labels: no array probe.",
             x=0.02, y=0.99, ha="left", fontsize=10, fontweight="bold", color=INK)
fig.savefig("figures/fig2_differential_network.pdf", dpi=220, bbox_inches="tight"); plt.close()

# ============================== FIG 3 =========================================
C = pd.read_csv("results/module_coherence.csv")
B = pd.read_csv("results/module_between_correlation.csv")
W = pd.read_csv("results/wgcna_33gene_membership.csv")
fig = plt.figure(figsize=(14, 4.6)); gs = fig.add_gridspec(1, 3, width_ratios=[1.25, .9, 1.2], wspace=.45)
ax = fig.add_subplot(gs[0])
dsn = {"ARR": "Microarray (n=195)", "RNA": "RNA-seq (n=22)", "SEA": "SEA-AD (n=20)"}
c = C[C.subset == "all"]
for i, m in enumerate(MODS):
    for j, d in enumerate(dsn):
        r = c[(c.module == m) & (c.dataset == d)].iloc[0]
        y = i * 4 + j
        ax.plot([r.null_mean, r.null_95pct], [y, y], color=GRID, lw=5, solid_capstyle="round", zorder=1)
        ax.scatter(r.mean_r, y, s=60, color=MCOL[m], edgecolor=SURF, lw=1.5, zorder=3, marker="os^"[j])
        ax.text(max(r.mean_r, r.null_95pct) + .03, y, f"p={r.p_empirical:.3f}" if r.p_empirical >= .001 else "p<0.001",
                va="center", fontsize=7.5, color=INK2)
ax.set_yticks([i * 4 + j for i in range(3) for j in range(3)], [f"{m} - {dsn[d]}" for m in MODS for d in dsn], fontsize=7.5)
ax.invert_yaxis(); ax.set_xlim(-.1, 1.0); ax.axvline(0, color=GRID, lw=.8)
ax.set_xlabel("Mean within-module Spearman r")
ax.set_title("A  Within-module coherence\n(grey bar: matched-random null, mean to 95th pct)")
ax = fig.add_subplot(gs[1])
b = B[B.dataset == "ARR"].set_index("module_A")[MODS].loc[MODS]
im = ax.imshow(b.values, cmap=DIV, vmin=-.8, vmax=.8)
for i in range(3):
    for j in range(3): ax.text(j, i, f"{b.values[i, j]:.2f}", ha="center", va="center", fontsize=9, color=INK)
ax.set_xticks(range(3), MODS, rotation=30, ha="right"); ax.set_yticks(range(3), MODS)
ax.spines[:].set_visible(False); ax.tick_params(length=0)
ax.set_title("B  Mean r within / between\nmodules (microarray)")
ax = fig.add_subplot(gs[2])
order = ["yellow", "turquoise", "red", "blue", "brown", "green", "magenta"]
lab = {"yellow": "yellow\n(neuronal, AD-)", "turquoise": "turquoise\n(neuronal, AD-)", "red": "red\n(neuronal, AD-)",
       "blue": "blue\n(mixed, AD+)", "brown": "brown\n(vascular/astro, AD+)", "green": "green\n(vascular/astro)",
       "magenta": "magenta\n(microglia, AD+)"}
ct = pd.crosstab(W.predefined_module, W.wgcna_module).reindex(index=MODS, columns=order).fillna(0)
left = np.zeros(len(order))
for m in MODS:
    ax.barh(range(len(order)), ct.loc[m].values, left=left, color=MCOL[m], height=.65, edgecolor=SURF, lw=2, label=m)
    left += ct.loc[m].values
ax.set_yticks(range(len(order)), [lab[o] for o in order], fontsize=7.5); ax.invert_yaxis()
ax.set_xlabel("Number of the 33 genes")
ax.legend(frameon=False, fontsize=8, loc="lower right", title="Predefined module", title_fontsize=8)
ax.set_title("C  Unsupervised WGCNA modules (5,000 genes)\ncontaining the 33 genes")
fig.savefig("figures/fig3_modules.pdf", dpi=220, bbox_inches="tight"); plt.close()

# ============================== FIG 4 =========================================
M = pd.read_csv("results/celltype_mean_cpm_33genes.csv", index_col=0)
subs = ["Exc_IT", "Exc_deep", "Inh_Pvalb", "Inh_Sst", "Inh_CGE", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC", "Vascular"]
CTL = {"Exc_IT": "Excitatory IT", "Exc_deep": "Excitatory deep", "Inh_Pvalb": "Pvalb/Chandelier", "Inh_Sst": "Sst/Sst Chodl", "Inh_CGE": "Vip/Lamp5/Sncg/Pax6", "Astrocyte": "Astrocyte", "Microglia": "Microglia/immune", "Oligodendrocyte": "Oligodendrocyte", "OPC": "OPC", "Vascular": "Vascular"}
ordr = [g for m in ["Glial", "Neuropeptide", "Synaptic"] for g in G33 if MOD[g] == m]
lin = M.loc[ordr, subs]; rel = lin.div(lin.max(1), axis=0)
comp = pd.read_csv("results/celltype_composition_dementia.csv", index_col=0)
Wd = pd.read_csv("results/celltype_within_celltype_DE.csv")
fig = plt.figure(figsize=(15, 8.2)); gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], height_ratios=[1, 1], wspace=.3, hspace=.45)
ax = fig.add_subplot(gs[:, 0])
seq = LinearSegmentedColormap.from_list("seq", ["#f4f3f0", "#a9c6ee", "#2a78d6", "#0d3b78"])
ax.imshow(rel.values, cmap=seq, aspect="auto", vmin=0, vmax=1)
ax.set_xticks(range(len(subs)), [CTL[c] for c in subs], rotation=55, ha="right", fontsize=7.5)
ax.set_yticks(range(len(ordr)), ordr, fontsize=7.5)
for i, g in enumerate(ordr): ax.get_yticklabels()[i].set_color(MCOL[MOD[g]])
ax.spines[:].set_visible(False); ax.tick_params(length=0)
for y in [3.5, 12.5]: ax.axhline(y, color=SURF, lw=3)
ax.set_title("A  Where the 33 genes are expressed (SEA-AD MTG, donor x cell-type pseudobulk)\nexpression relative to each gene's maximum; label colour = predefined module")
ax = fig.add_subplot(gs[0, 1])
cs = ["Inh_Sst", "Inh_Pvalb", "Inh_CGE", "Exc_IT", "Exc_deep", "Astrocyte", "Microglia", "Oligodendrocyte", "Vascular"]
cc = comp.loc[cs]
ax.barh(range(len(cs)), cc.log2_ratio, color=[DIV(.15) if v < 0 else DIV(.85) for v in cc.log2_ratio], height=.6)
for i, (v, p) in enumerate(zip(cc.log2_ratio, cc.p_mannwhitney)):
    ax.text(v + (.06 if v >= 0 else -.06), i, f"p={p:.3f}", va="center", ha="left" if v >= 0 else "right", fontsize=7.5, color=INK2)
ax.set_yticks(range(len(cs)), [CTL[c] for c in cs]); ax.invert_yaxis(); ax.axvline(0, color=INK2, lw=.8); ax.set_xlim(-1.9, 2.2)
ax.set_xlabel("log2 (proportion in dementia / no dementia)")
ax.set_title("B  Cell-type composition shift (10 vs 10 donors)")
ax = fig.add_subplot(gs[1, 1])
keyg = ["GFAP", "GJA1", "VIM", "TREM2", "SST", "NPY", "TAC1", "CARTPT", "NPTX2", "CHGB"]
rows = [Wd[Wd.gene == g].iloc[0] for g in keyg]
y = np.arange(len(rows))
ax.barh(y - .19, [r.whole_tissue_log2FC_mle for r in rows], .36, color=MUTED, label="bulk (all nuclei)")
ax.barh(y + .19, [r.home_log2FC_mle for r in rows], .36, color=CAT[0], label="within the expressing cell type")
ax.set_yticks(y, [f"{r.gene}  ({CTL[r.home_celltype]})" for r in rows], fontsize=7.5); ax.invert_yaxis()
ax.axvline(0, color=INK2, lw=.8); ax.set_xlabel("log2FC dementia vs no dementia")
ax.legend(frameon=False, fontsize=7.5, loc="lower right")
ax.set_title("C  Bulk change vs change inside the cell type\n(glial 'up-regulation' is mostly more glial cells)")
fig.savefig("figures/fig4_celltype.pdf", dpi=220, bbox_inches="tight"); plt.close()

# ============================== FIG 5 =========================================
T = pd.read_csv("results/tf_enrichment_trrust.csv"); T = T[T.background == "TRRUST targets on microarray"]   # stricter background; the all-genes background inflates (see README)
tfs = T[(T.gene_set == "33-gene panel") & (T.fdr < .05)]
rest = T[(T.gene_set == "concordant DOWN in AD") & (T.TF == "REST")].iloc[0]
g = nx.DiGraph(); modecol = {"A": CAT[2], "R": CAT[1], "U": MUTED}
edges = []
for r in tfs.itertuples():
    for t in r.targets.split(";"): edges.append((r.TF, t.split("(")[0], t[-2]))
for t in rest.targets.split(";"):
    e = ("REST", t.split("(")[0], t[-2])
    if e[:2] not in [x[:2] for x in edges]: edges.append(e)
for a, b, m in edges: g.add_edge(a, b, mode=m)
tfl = list(dict.fromkeys(a for a, _, _ in edges)); tg = sorted({b for _, b, _ in edges}, key=lambda x: (x not in G33, x))
P = {t: (0, -i * 1.6 - 1) for i, t in enumerate(tfl)}
P.update({t: (4, -i * 0.62) for i, t in enumerate(tg)})
fig, ax = plt.subplots(figsize=(8.5, 7.5))
for a, b, d in g.edges(data=True):
    ax.annotate("", xy=P[b], xytext=P[a], arrowprops=dict(arrowstyle="-|>" if d["mode"] != "R" else "-[",
                color=modecol[d["mode"]], lw=1.4, shrinkA=14, shrinkB=8, alpha=.9))
lf = panel["log2FC_micro"]
for t in tfl:
    ax.scatter(*P[t], s=900, marker="s", color="#f0efec", edgecolor=INK, zorder=3)
    ax.text(*P[t], t, ha="center", va="center", fontweight="bold", fontsize=9, zorder=4)
for t in tg:
    in33 = t in G33
    ax.scatter(*P[t], s=70, color=MCOL[MOD[t]] if in33 else "#f0efec", edgecolor=INK, zorder=3)
    ax.text(P[t][0] + .2, P[t][1], t + ("" if in33 else "  (signature)"), va="center", fontsize=8, color=INK if in33 else INK2)
info = {r.TF: (f"FDR={r.fdr:.1e}" if r.fdr < 1e-3 else f"FDR={r.fdr:.3f}") for r in tfs.itertuples()}
for t in tfl:
    ax.text(P[t][0] - .55, P[t][1], info.get(t, f"FDR={rest.fdr:.1e}\n(signature)"), ha="right", va="center", fontsize=7.5, color=INK2)
ax.legend(handles=[Line2D([], [], color=CAT[2], lw=1.6, label="activation"), Line2D([], [], color=CAT[1], lw=1.6, label="repression"),
                   Line2D([], [], color=MUTED, lw=1.6, label="mode unknown"),
                   Line2D([], [], marker="o", ls="", color=MCOL["Synaptic"], label="33-gene set: Synaptic"),
                   Line2D([], [], marker="o", ls="", color=MCOL["Neuropeptide"], label="33-gene set: Neuropeptide"),
                   Line2D([], [], marker="o", ls="", color=MCOL["Glial"], label="33-gene set: Glial"),
                   Line2D([], [], marker="o", ls="", markerfacecolor="#f0efec", markeredgecolor=INK, label="other concordant DOWN gene")],
          frameon=False, fontsize=7.5, loc="lower left", bbox_to_anchor=(0, -0.02))
ax.set_xlim(-2.2, 6.5); ax.axis("off")
ax.set_title("TRRUST v2 regulators enriched among the 33 genes (FDR<0.05; background = TRRUST targets on array)\nREST also top TF for the concordant DOWN-in-AD signature (" + str(int((pd.read_csv("results/concordant_signature.csv").mean_log2FC < 0).sum())) + " genes, " + str(int(rest.set_size)) + " in TRRUST)")
fig.savefig("figures/fig5_tf_network.pdf", dpi=220, bbox_inches="tight"); plt.close()

# ============================== FIG 6 =========================================
A = pd.read_csv("results/attack_tolerance.csv")
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6), gridspec_kw={"width_ratios": [1, 1.1], "wspace": .35})
ax = axs[0]
for col, c, lab in [("random_failure", CAT[2], "random removal (mean of 500)"), ("hub_attack", CAT[1], "remove hubs (static rank)"),
                    ("adaptive_betweenness_attack", CAT[0], "remove highest-betweenness (adaptive)")]:
    ax.plot(A.n_removed, A[col], color=c, lw=2, label=lab)
ax.axhline(.5, color=GRID, lw=1, ls="--")
ax.set_xlabel("Genes removed"); ax.set_ylabel("Largest connected component (fraction of 33)")
ax.legend(frameon=False, fontsize=8); ax.set_xlim(0, 20); ax.set_ylim(0, 1.02)
ax.set_title("A  Robustness: random failure vs targeted attack\n(bridge genes, not high-degree genes, hold the network together)")
ax = axs[1]
top = H.sort_values("hub_rank_all").head(12)
y = np.arange(len(top))
for col, c, mk, lab in [("degree_all", CAT[0], "o", "all STRING evidence (>=0.4)"),
                        ("degree_hc700", CAT[2], "^", "high confidence only (>=0.7)"),
                        ("degree_expdb", CAT[1], "s", "experimental/database only")]:
    ax.scatter(top[col], y, color=c, marker=mk, s=45, label=lab, zorder=3, edgecolor=SURF, lw=1)
for i, g in enumerate(top.index): ax.plot([top.loc[g, ["degree_all", "degree_hc700", "degree_expdb"]].min(), top.loc[g, "degree_all"]], [i, i], color=GRID, lw=2, zorder=1)
ax.set_yticks(y, top.index); ax.invert_yaxis(); ax.set_xlabel("Degree")
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.set_title("B  Hub degree under stricter evidence filters\n(SNAP25, NEFL, GFAP, CAMK2A stay hubs)")
fig.savefig("figures/fig6_robustness.pdf", dpi=220, bbox_inches="tight"); plt.close()
print("Saved figures/fig1-6")
