"""
p2_03_tf_enrichment.py
TF enrichment with TRRUST v2 (reused local file), one-sided Fisher's exact test.

  - Background 1 = all genes measured on the microarray (19,743). Background 2
    (stricter) = TRRUST targets on the array; a whole-genome background inflates
    enrichment for TFs whose curated targets are mostly neuronal.
  - Gene sets:
      * 33-gene panel
      * concordant signature: same sign of log2FC in all three datasets AND
        significant (padj < 0.05, |log2FC| > 0.5) in >= 2 datasets, split DOWN / UP
  - TF log2FC from the DEG tables (microarray limma; RNA-seq and SEA-AD apeglm-shrunk,
    so SEA-AD values are compressed towards 0); an unshrunk SEA-AD estimate for REST
    is added from a DESeq2 MLE refit.
Run from network/.
"""
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

F = "../analysis"
P33 = pd.read_csv(f"{F}/comparison/results/panel_33gene_summary.csv").gene.tolist()
D = {k: pd.read_csv(f"{F}/{p}").drop_duplicates("gene").set_index("gene") for k, p in
     {"ARR": "dataset1_microarray/results/microarray_DEGs.csv", "RNA": "dataset2_rnaseq/results/rnaseq_DEGs.csv",
      "SEA": "dataset3_seaad/results/seaad_DEGs.csv"}.items()}
tr = pd.read_csv("data/trrust_rawdata.human.tsv", sep="\t", header=None, names=["TF", "target", "mode", "pmid"])
arr = pd.read_csv("data/expr_ARR.csv.gz", index_col=0); lab = pd.read_csv("data/labels_ARR.csv", index_col=0).group

# concordant signature
common = D["ARR"].index.intersection(D["RNA"].index).intersection(D["SEA"].index)
L = pd.DataFrame({k: D[k].log2FC[common] for k in D}); S = pd.DataFrame({k: D[k].direction[common] != "NS" for k in D})
cons = (np.sign(L).nunique(axis=1) == 1) & (S.sum(1) >= 2)
sig = L[cons].mean(1)
down, up = set(sig.index[sig < 0]), set(sig.index[sig > 0])
pd.DataFrame({"mean_log2FC": sig, "n_sig": S.sum(1)[cons]}).sort_values("mean_log2FC").to_csv("results/concordant_signature.csv")
print(f"Concordant signature: {cons.sum()} genes (down {len(down)}, up {len(up)}); panel genes in it: {sorted(set(P33) & set(sig.index))}")

sets = {"33-gene panel": set(P33), "concordant DOWN in AD": down, "concordant UP in AD": up}
def run(bg_name, bg):
    t = tr[tr.target.isin(bg)]
    reg = t.groupby("TF").target.apply(set); reg = reg[reg.apply(len) >= 5]
    rows = []
    for sname, Sg in sets.items():
        Sg = Sg & bg; N, n = len(bg), len(Sg)
        for tf, T in reg.items():
            k = len(Sg & T)
            if k == 0: continue
            table = [[k, n - k], [len(T) - k, N - n - len(T) + k]]
            p = stats.fisher_exact(table, alternative="greater")[1]
            modes = t[(t.TF == tf) & t.target.isin(Sg & T)].groupby("target")["mode"].first()
            rows.append(dict(background=bg_name, gene_set=sname, set_size=n, TF=tf, regulon_size=len(T), overlap=k,
                             fold_enrichment=round(k / (n * len(T) / N), 2), p=p,
                             targets=";".join(f"{g}({modes[g][0]})" for g in sorted(Sg & T))))
    R = pd.DataFrame(rows)
    R["fdr"] = R.groupby("gene_set").p.transform(lambda p: multipletests(p, method="fdr_bh")[1])
    return R.sort_values(["gene_set", "p"])

bg_all = set(arr.index)
R = pd.concat([run("all microarray genes", bg_all), run("TRRUST targets on microarray", bg_all & set(tr.target))])
for k in D:
    R[f"TF_log2FC_{k}"] = R.TF.map(D[k].log2FC); R[f"TF_padj_{k}"] = R.TF.map(D[k].padj)
def tf_r(tf, targets, grp):
    if tf not in arr.index: return np.nan
    X = arr.loc[:, lab.reindex(arr.columns) == grp]
    rs = [stats.spearmanr(X.loc[tf], X.loc[g])[0] for g in targets if g in X.index]
    return np.mean(rs) if rs else np.nan
tg = [[x.split("(")[0] for x in s.split(";")] for s in R.targets]; ok = (R.p < 0.05).values
R["mean_r_TF_targets_ctrl"] = [tf_r(t, g, "ND") if o else np.nan for t, g, o in zip(R.TF, tg, ok)]
R["mean_r_TF_targets_AD"] = [tf_r(t, g, "AD") if o else np.nan for t, g, o in zip(R.TF, tg, ok)]
R.to_csv("results/tf_enrichment_trrust.csv", index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
show = ["TF", "overlap", "regulon_size", "fold_enrichment", "p", "fdr", "TF_log2FC_ARR", "TF_log2FC_RNA", "TF_log2FC_SEA", "targets"]
for b in R.background.unique():
    for s in sets:
        print(f"\n=== {s}  |  background: {b} ===")
        print(R[(R.background == b) & (R.gene_set == s)].head(6)[show].round(4).to_string(index=False))

# REST in all three datasets (+ unshrunk SEA-AD estimate, which apeglm does not compress)
rest = pd.DataFrame({k: {"log2FC": D[k].log2FC.get("REST"), "padj": D[k].padj.get("REST"), "direction": D[k].direction.get("REST")} for k in D}).T
try:
    import subprocess, io
    out = subprocess.run(["Rscript", "-e",
        'suppressMessages(library(DESeq2)); pb <- readRDS("../analysis/dataset3_seaad/results/seaad_pseudobulk.rds");'
        'd <- DESeqDataSetFromMatrix(pb$counts, pb$coldata, ~condition); d$condition <- relevel(d$condition,"control");'
        'd <- d[rowSums(counts(d)>=10)>=5,]; d <- DESeq(d, quiet=TRUE); r <- results(d, name="condition_AD_vs_control");'
        'cat(r["REST","log2FoldChange"], r["REST","padj"])'], capture_output=True, text=True, timeout=900)
    v = out.stdout.split(); rest.loc["SEA (unshrunk MLE)"] = [float(v[0]), float(v[1]), ""]
except Exception as ex:
    print("unshrunk SEA-AD REST estimate unavailable:", ex)
rest.to_csv("results/REST_log2FC_three_datasets.csv")
print("\nREST log2FC (AD - control):"); print(rest)
cov = tr[tr.target.isin(P33)].groupby("target").TF.apply(lambda x: ",".join(sorted(set(x))))
cov.to_csv("results/trrust_regulators_of_33genes.csv")
print(f"\nTRRUST covers {len(cov)}/33 panel genes")
