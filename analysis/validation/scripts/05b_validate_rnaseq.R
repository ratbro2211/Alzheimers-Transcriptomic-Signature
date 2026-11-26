# =============================================================================
# 05b_validate_rnaseq.R — validation of Dataset 2 (GSE104704, DESeq2 + apeglm)
# Rebuilds raw counts exactly as script 02 and tests correctness and bias.
# Writes results/checks_rnaseq.csv, leave-one-out tables, figures.
# =============================================================================
# Repository root: taken from AD_PROJECT_ROOT, otherwise found by walking up from the
# working directory to the folder that contains the .ad_project_root marker file.
ad_root <- function() {
  r <- Sys.getenv("AD_PROJECT_ROOT"); if (nzchar(r)) return(normalizePath(r))
  d <- normalizePath(getwd())
  repeat {
    if (file.exists(file.path(d, ".ad_project_root"))) return(d)
    p <- dirname(d); if (p == d) stop("Run inside the repository or set AD_PROJECT_ROOT"); d <- p
  }
}

set.seed(42); options(width = 200)
proj <- ad_root()
vdir <- file.path(proj, "analysis", "validation")
suppressPackageStartupMessages({ library(DESeq2); library(apeglm); library(ggplot2) })
sink(file.path(vdir, "results", "05b_log.txt"), split = TRUE)
checks <- data.frame()
add <- function(id, check, value, verdict, note = "") {
  checks <<- rbind(checks, data.frame(dataset = "RNA-seq", id = id, check = check, value = value, verdict = verdict, note = note))
  cat(sprintf("[%s] %-4s %-48s %s  %s\n", verdict, id, check, value, note))
}
PADJ <- 0.05; LFC <- 0.5
MARK_DN <- c("VGF","SST","NPTX2","BDNF","SCG2","CHGB","PCSK1","CORT","CRH","RGS4")
MARK_UP <- c("GFAP","CD44","SERPINA3","CHI3L1","C4B","AQP4")

# ---- rebuild data as in script 02 -------------------------------------------
gdir <- file.path(proj, "data", "GSE104704")
tab  <- read.delim(file.path(gdir, "GSE104704_RNA-Seq_Table.txt.gz"), row.names = 1, check.names = FALSE)
norm <- as.matrix(tab[, setdiff(colnames(tab), c("Young_vs_Old", "Old_vs_AD", "Young_vs_AD"))])
raw  <- sweep(norm, 2, apply(norm, 2, function(v) min(v[v > 0])), "/")
dev  <- max(abs(raw - round(raw))); counts <- round(raw); storage.mode(counts) <- "integer"
con <- gzfile(file.path(gdir, "GSE104704_series_matrix.txt.gz"), "rt"); hdr <- readLines(con); close(con)
fld <- function(tag) gsub('"', '', strsplit(hdr[grep(tag, hdr)][1], "\t")[[1]][-1])
titles <- fld("^!Sample_title"); ages <- as.numeric(sub(".*: ", "", fld("age \\(years\\)")))
meta <- data.frame(sample = sub("\\.[^.]+$", "", titles), grp = sub(".*\\.", "", titles), age = ages)
meta <- meta[meta$grp %in% c("Old", "AD"), ]; rownames(meta) <- meta$sample
meta$condition <- factor(ifelse(meta$grp == "AD", "AD", "control"), levels = c("control", "AD"))
counts <- counts[, meta$sample]
deg <- read.csv(file.path(proj, "analysis/dataset2_rnaseq/results/rnaseq_DEGs.csv"))

run <- function(cn, cd, design = ~ condition, shrink = FALSE) {
  dds <- DESeqDataSetFromMatrix(cn, cd, design); dds$condition <- relevel(dds$condition, "control")
  dds <- dds[rowSums(counts(dds) >= 10) >= min(table(dds$condition)), ]
  dds <- DESeq(dds, quiet = TRUE)
  r <- results(dds, name = "condition_AD_vs_control", alpha = 0.05)
  out <- data.frame(gene = rownames(r), lfc_mle = r$log2FoldChange, pvalue = r$pvalue, padj = r$padj)
  if (shrink) out$log2FC <- lfcShrink(dds, coef = "condition_AD_vs_control", type = "apeglm", quiet = TRUE)$log2FoldChange
  attr(out, "dds") <- dds; out
}
nsig <- function(d, col = "lfc_mle") sum(!is.na(d$padj) & d$padj < PADJ & abs(d[[col]]) > LFC)

# ---- V1 labels + counts --------------------------------------------------------
add("V1", "Sample labels (AD vs age-matched Old)", sprintf("AD %d / control %d; Young excluded", sum(meta$condition == "AD"), sum(meta$condition == "control")),
    ifelse(all(table(meta$condition) == c(10, 12)), "PASS", "FAIL"))
add("V1b", "Raw counts recovered exactly from GEO table", sprintf("max deviation from integer = %.1e", dev), ifelse(dev < 1e-6, "PASS", "FAIL"))

full <- run(counts, meta)
dds <- attr(full, "dds")
stopifnot(nrow(full) == nrow(deg))

# ---- V2 marker direction -------------------------------------------------------
m <- deg[deg$gene %in% c(MARK_DN, MARK_UP), ]; m$expected <- ifelse(m$gene %in% MARK_UP, 1, -1)
ok2 <- sum(sign(m$log2FC) == m$expected)
okd <- sum(sign(m$log2FC[m$expected < 0]) < 0); oku <- sum(sign(m$log2FC[m$expected > 0]) > 0)
add("V2", "Canonical AD markers in expected direction", sprintf("%d/%d (neuronal down %d/%d; glial up %d/%d)", ok2, nrow(m),
    okd, sum(m$expected < 0), oku, sum(m$expected > 0)), ifelse(ok2 / nrow(m) >= 0.8, "PASS", ifelse(ok2 / nrow(m) >= 0.6, "WARN", "FAIL")))

# ---- V3 covariates / library size ---------------------------------------------
pa <- wilcox.test(age ~ condition, meta)$p.value
add("V3a", "Age balanced between AD and control", sprintf("median AD %.0f vs control %.0f, p = %.2f",
    median(meta$age[meta$condition == "AD"]), median(meta$age[meta$condition == "control"]), pa), ifelse(pa < 0.05, "WARN", "PASS"), "sex not provided in GEO")
sf <- sizeFactors(dds); lib <- colSums(counts)
pl <- wilcox.test(lib ~ meta$condition)$p.value
add("V3b", "Library size not confounded with condition", sprintf("size factors %.2f-%.2f; library size p = %.2f", min(sf), max(sf), pl), ifelse(pl < 0.05, "WARN", "PASS"))
adj <- run(counts, meta, ~ age + condition)
r_age <- cor(full$lfc_mle, adj$lfc_mle[match(full$gene, adj$gene)], use = "complete.obs")
add("V4", "Age-adjusted model changes little", sprintf("r(log2FC) = %.3f; DEGs %d -> %d", r_age, nsig(full), nsig(adj)),
    ifelse(r_age > 0.95, "PASS", ifelse(r_age > 0.85, "WARN", "FAIL")))

# ---- V5 p-values / Cook's outliers --------------------------------------------
p <- full$pvalue[!is.na(full$pvalue)]; pi0 <- min(1, mean(p > 0.5) / 0.5); tr <- mean(p > 0.9) / mean(p > 0.5 & p <= 0.6)
add("V5", "p-value histogram: uniform tail + low-p spike", sprintf("pi0 = %.2f; tail ratio = %.2f", pi0, tr), ifelse(tr > 1.3 | tr < 0.7, "WARN", "PASS"))
write.csv(data.frame(dataset = "RNA-seq", p = p), file.path(vdir, "results", "pvals_rnaseq.csv"), row.names = FALSE)
cooks <- assays(dds)[["cooks"]]; cut <- qf(0.99, 2, ncol(dds) - 2)
nck <- sum(apply(cooks, 1, max) > cut, na.rm = TRUE)
worst <- names(sort(table(colnames(cooks)[apply(cooks, 1, which.max)][apply(cooks, 1, max) > cut]), decreasing = TRUE))[1:3]
add("V5b", "Cook's distance outlier genes", sprintf("%d genes (%.1f%%); most often driven by %s", nck, 100 * nck / nrow(cooks), paste(worst, collapse = ", ")),
    ifelse(nck / nrow(cooks) < 0.05, "PASS", "WARN"), "DESeq2 replaces/flags these automatically")

# ---- V6 permutation null ------------------------------------------------------
NPERM <- 20
nullsig <- sapply(seq_len(NPERM), function(i) { mp <- meta; mp$condition <- sample(mp$condition); nsig(run(counts, mp)) })
pe <- (sum(nullsig >= nsig(full)) + 1) / (length(nullsig) + 1)
add("V6", "Label-permutation null: DEGs under shuffled labels",
    sprintf("observed %d; null median %.0f, max %d; %d/%d perms >= observed; empirical p = %.3f",
            nsig(full), median(nullsig), max(nullsig), sum(nullsig >= nsig(full)), length(nullsig), pe),
    if (pe > 0.05) "FAIL" else if (any(nullsig >= 0.5 * nsig(full))) "WARN" else "PASS",
    "FAIL if empirical p > 0.05; WARN if any permutation reaches 50% of observed")
write.csv(data.frame(dataset = "RNA-seq", null_DEGs = nullsig, observed = nsig(full)), file.path(vdir, "results", "permnull_rnaseq.csv"), row.names = FALSE)

# ---- V7 PCA / outliers --------------------------------------------------------
v <- assay(vst(dds, blind = TRUE)); pc <- prcomp(t(v[order(-apply(v, 1, var))[1:2000], ]))
ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
d <- sqrt(rowSums(scale(pc$x[, 1:3])^2)); out <- names(d)[d > sqrt(qchisq(0.99, 3))]
add("V7", "Outlier samples (Mahalanobis PC1-3, p<0.01)", sprintf("%d: %s; condition p on PC1 = %.2g", length(out), paste(out, collapse = ", "),
    wilcox.test(pc$x[, 1] ~ meta$condition)$p.value), ifelse(length(out) == 0, "PASS", "WARN"))
pcd <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], condition = meta$condition, s = rownames(pc$x))
g <- ggplot(pcd, aes(PC1, PC2, colour = condition, label = s)) + geom_point(size = 2.5) +
  ggrepel::geom_text_repel(size = 2.6, show.legend = FALSE, seed = 42) +
  scale_colour_manual(values = c(control = "#457B9D", AD = "#E63946")) +
  labs(title = "RNA-seq PCA (VST, top 2,000 genes)", x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2])) +
  theme_bw() + theme(plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(vdir, "figures", "pca_rnaseq.png"), g, width = 6.5, height = 5, dpi = 300, bg = "white")

# ---- V8 leave-one-sample-out --------------------------------------------------
loo <- do.call(rbind, lapply(meta$sample, function(s) {
  k <- setdiff(meta$sample, s); r <- run(counts[, k], meta[k, ])
  mm <- match(full$gene, r$gene)
  data.frame(left_out = s, condition = meta[s, "condition"], n_DEGs = nsig(r),
             r_lfc_vs_full = cor(full$lfc_mle, r$lfc_mle[mm], use = "complete.obs"))
}))
loo$pct_change_DEGs <- round(100 * (loo$n_DEGs / nsig(full) - 1), 1)
write.csv(loo, file.path(vdir, "results", "loo_rnaseq.csv"), row.names = FALSE)
infl <- loo$left_out[which.max(abs(loo$pct_change_DEGs))]
add("V8", "Leave-one-out: no single sample drives the result", sprintf("min r = %.3f; DEG change %+.0f%% to %+.0f%%; most influential %s",
    min(loo$r_lfc_vs_full), min(loo$pct_change_DEGs), max(loo$pct_change_DEGs), infl),
    ifelse(min(loo$r_lfc_vs_full) > 0.9 & max(abs(loo$pct_change_DEGs)) < 50, "PASS", "WARN"))
# shrunk results without the most influential sample, for the comparison check
k <- setdiff(meta$sample, infl); r <- run(counts[, k], meta[k, ], shrink = TRUE)
r$direction <- ifelse(!is.na(r$padj) & r$padj < PADJ & r$log2FC > LFC, "Up", ifelse(!is.na(r$padj) & r$padj < PADJ & r$log2FC < -LFC, "Down", "NS"))
write.csv(r[, c("gene", "log2FC", "pvalue", "padj", "direction")], file.path(vdir, "results", sprintf("rnaseq_DEGs_drop_%s.csv", infl)), row.names = FALSE)
writeLines(infl, file.path(vdir, "results", "rnaseq_most_influential.txt"))

write.csv(checks, file.path(vdir, "results", "checks_rnaseq.csv"), row.names = FALSE)
sink()
