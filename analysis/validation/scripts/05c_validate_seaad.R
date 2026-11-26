# =============================================================================
# 05c_validate_seaad.R — validation of Dataset 3 (SEA-AD pseudobulk, DESeq2)
# Uses the pseudobulk saved by script 03 (raw UMI sums) and tests correctness,
# covariate balance, outlier influence and cell-composition bias.
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
vdir <- file.path(proj, "analysis", "validation"); d3 <- file.path(proj, "analysis", "dataset3_seaad")
suppressPackageStartupMessages({ library(DESeq2); library(apeglm); library(ggplot2) })
sink(file.path(vdir, "results", "05c_log.txt"), split = TRUE)
checks <- data.frame()
add <- function(id, check, value, verdict, note = "") {
  checks <<- rbind(checks, data.frame(dataset = "SEA-AD", id = id, check = check, value = value, verdict = verdict, note = note))
  cat(sprintf("[%s] %-4s %-48s %s  %s\n", verdict, id, check, value, note))
}
PADJ <- 0.05; LFC <- 0.5
MARK_DN <- c("VGF","SST","NPTX2","BDNF","SCG2","CHGB","PCSK1","CORT","CRH","RGS4")
MARK_UP <- c("GFAP","CD44","SERPINA3","CHI3L1","C4B","AQP4")

pb <- readRDS(file.path(d3, "results", "seaad_pseudobulk.rds")); counts <- pb$counts; meta <- pb$coldata
meta$age <- as.numeric(sub("\\+", "", meta$age))
cc <- read.csv(file.path(d3, "results", "seaad_celltype_counts_per_donor.csv"), row.names = 1, check.names = FALSE)[rownames(meta), ]
neur <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L5 ET","L5/6 NP","L6 CT","L6b","Pvalb","Chandelier","Sst","Sst Chodl","Vip","Lamp5","Lamp5 Lhx6","Sncg","Pax6")
meta$neuronal_frac <- rowSums(cc[, intersect(neur, colnames(cc))]) / rowSums(cc)
meta$nf_z <- as.numeric(scale(meta$neuronal_frac))
deg <- read.csv(file.path(d3, "results", "seaad_DEGs.csv"))

run <- function(cn, cd, design = ~ condition, shrink = FALSE) {
  dds <- DESeqDataSetFromMatrix(cn, cd, design); dds$condition <- relevel(dds$condition, "control")
  dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]; dds <- DESeq(dds, quiet = TRUE)
  r <- results(dds, name = "condition_AD_vs_control", alpha = 0.05)
  out <- data.frame(gene = rownames(r), lfc_mle = r$log2FoldChange, pvalue = r$pvalue, padj = r$padj)
  if (shrink) out$log2FC <- lfcShrink(dds, coef = "condition_AD_vs_control", type = "apeglm", quiet = TRUE)$log2FoldChange
  attr(out, "dds") <- dds; out
}
nsig <- function(d) sum(!is.na(d$padj) & d$padj < PADJ & abs(d$lfc_mle) > LFC)

# ---- V1 labels + raw counts --------------------------------------------------
add("V1", "Donor labels (cognitive status)", sprintf("AD %d / control %d", sum(meta$condition == "AD"), sum(meta$condition == "control")),
    ifelse(all(table(meta$condition) == 10), "PASS", "FAIL"))
add("V1b", "Input is raw integer UMI counts", sprintf("integer = %s; max per-donor gene count = %d", all(counts == round(counts)), max(counts)),
    ifelse(all(counts == round(counts)) && max(counts) > 1000, "PASS", "FAIL"))
hi_ctrl <- meta$donor_id[meta$condition == "control" & meta$ADNC %in% c("High", "Intermediate")]
add("V1c", "Controls free of AD neuropathology", sprintf("%d/10 controls have Intermediate/High ADNC (%s)", length(hi_ctrl), paste(hi_ctrl, collapse = ", ")),
    ifelse(length(hi_ctrl) == 0, "PASS", "WARN"), "labels are clinical (dementia), not pathological; dilutes contrast")

full <- run(counts, meta); dds <- attr(full, "dds")

# ---- V2 marker direction -----------------------------------------------------
m <- deg[deg$gene %in% c(MARK_DN, MARK_UP), ]; m$expected <- ifelse(m$gene %in% MARK_UP, 1, -1)
mm <- full[match(m$gene, full$gene), ]
ok2 <- sum(sign(mm$lfc_mle) == m$expected); nsm <- sum(!is.na(m$padj) & m$padj < PADJ)
add("V2", "Canonical AD markers in expected direction", sprintf("%d/%d (unshrunk log2FC); significant: %d", ok2, nrow(m), nsm),
    ifelse(ok2 / nrow(m) >= 0.8, "PASS", ifelse(ok2 / nrow(m) >= 0.6, "WARN", "FAIL")), "apeglm shrinks most of these to ~0")

# ---- V3 covariates ------------------------------------------------------------
pa <- wilcox.test(age ~ condition, meta)$p.value; ps <- fisher.test(table(meta$sex, meta$condition))$p.value
add("V3a", "Age balanced", sprintf("median AD %.0f vs control %.0f, p = %.2f", median(meta$age[meta$condition == "AD"]),
    median(meta$age[meta$condition == "control"]), pa), ifelse(pa < 0.05, "WARN", "PASS"))
add("V3b", "Sex balanced", sprintf("female AD %d/10 vs control %d/10, p = %.2f", sum(meta$sex[meta$condition == "AD"] == "Female"),
    sum(meta$sex[meta$condition == "control"] == "Female"), ps), ifelse(ps < 0.05, "WARN", "PASS"))
adj <- run(counts, meta, ~ age + sex + condition)
r_cov <- cor(full$lfc_mle, adj$lfc_mle[match(full$gene, adj$gene)], use = "complete.obs")
add("V4", "Age + sex adjustment changes little", sprintf("r(log2FC) = %.3f; DEGs %d -> %d", r_cov, nsig(full), nsig(adj)),
    ifelse(r_cov > 0.95, "PASS", ifelse(r_cov > 0.85, "WARN", "FAIL")))

# ---- V5 p-values ---------------------------------------------------------------
p <- full$pvalue[!is.na(full$pvalue)]; pi0 <- min(1, mean(p > 0.5) / 0.5); tr <- mean(p > 0.9) / mean(p > 0.5 & p <= 0.6)
add("V5", "p-value histogram: uniform tail + low-p spike", sprintf("pi0 = %.2f; tail ratio = %.2f", pi0, tr), ifelse(tr > 1.3 | tr < 0.7, "WARN", "PASS"))
write.csv(data.frame(dataset = "SEA-AD", p = p), file.path(vdir, "results", "pvals_seaad.csv"), row.names = FALSE)

# ---- V6 composition bias (the key risk for whole-tissue pseudobulk) -------------
pn <- wilcox.test(neuronal_frac ~ condition, meta)$p.value
add("V6a", "Neuronal fraction balanced between groups", sprintf("median AD %.0f%% vs control %.0f%%, p = %.3f",
    100 * median(meta$neuronal_frac[meta$condition == "AD"]), 100 * median(meta$neuronal_frac[meta$condition == "control"]), pn),
    ifelse(pn < 0.05, "WARN", "PASS"))
v <- assay(vst(dds, blind = TRUE))
gene_nf <- apply(v, 1, function(x) cor(x, meta$neuronal_frac, method = "spearman"))
rc <- cor(full$lfc_mle, gene_nf[full$gene], use = "complete.obs", method = "spearman")
sg <- !is.na(full$padj) & full$padj < PADJ & abs(full$lfc_mle) > LFC
rcs <- cor(full$lfc_mle[sg], gene_nf[full$gene[sg]], method = "spearman")
add("V6b", "DE log2FC explained by neuronal fraction", sprintf("rho(log2FC, gene~neuronal-fraction corr) = %.2f all genes, %.2f DEGs", rc, rcs),
    ifelse(abs(rc) > 0.5 | abs(rcs) > 0.5, "WARN", "PASS"), "|rho| > 0.5 = WARN; genome-wide fold changes strongly track neuronal fraction")
cmp <- run(counts, meta, ~ nf_z + condition)
r_nf <- cor(full$lfc_mle, cmp$lfc_mle[match(full$gene, cmp$gene)], use = "complete.obs")
add("V6c", "Adjusting for neuronal fraction", sprintf("r(log2FC) = %.2f; DEGs %d -> %d", r_nf, nsig(full), nsig(cmp)),
    ifelse(r_nf > 0.85, "PASS", "WARN"), "large change = composition drives much of the bulk signal")

# ---- V7 permutation null ------------------------------------------------------
NPERM <- 20
nullsig <- sapply(seq_len(NPERM), function(i) { mp <- meta; mp$condition <- sample(mp$condition); nsig(run(counts, mp)) })
pe <- (sum(nullsig >= nsig(full)) + 1) / (length(nullsig) + 1)
add("V7", "Label-permutation null: DEGs under shuffled labels",
    sprintf("observed %d; null median %.0f, max %d; %d/%d perms >= observed; empirical p = %.3f",
            nsig(full), median(nullsig), max(nullsig), sum(nullsig >= nsig(full)), length(nullsig), pe),
    if (pe > 0.05) "FAIL" else if (any(nullsig >= 0.5 * nsig(full))) "WARN" else "PASS",
    "FAIL if empirical p > 0.05; WARN if any permutation reaches 50% of observed")
write.csv(data.frame(dataset = "SEA-AD", null_DEGs = nullsig, observed = nsig(full)), file.path(vdir, "results", "permnull_seaad.csv"), row.names = FALSE)

# ---- V8 PCA / outliers --------------------------------------------------------
pc <- prcomp(t(v[order(-apply(v, 1, var))[1:2000], ])); ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
d <- sqrt(rowSums(scale(pc$x[, 1:3])^2)); out <- names(d)[d > sqrt(qchisq(0.99, 3))]
add("V8", "Outlier donors (Mahalanobis PC1-3, p<0.01)", sprintf("%d: %s; PC1 vs neuronal fraction rho = %.2f", length(out), paste(out, collapse = ", "),
    cor(pc$x[, 1], meta$neuronal_frac, method = "spearman")), ifelse(length(out) == 0, "PASS", "WARN"))
pcd <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], condition = meta$condition, nf = meta$neuronal_frac, s = rownames(pc$x))
g <- ggplot(pcd, aes(PC1, PC2, colour = condition, size = nf, label = s)) + geom_point(alpha = .85) +
  ggrepel::geom_text_repel(size = 2.6, show.legend = FALSE, seed = 42) + scale_size(range = c(1.5, 5), name = "neuronal\nfraction") +
  scale_colour_manual(values = c(control = "#457B9D", AD = "#E63946")) +
  labs(title = "SEA-AD pseudobulk PCA (VST, top 2,000 genes)", x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2])) +
  theme_bw() + theme(plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(vdir, "figures", "pca_seaad.png"), g, width = 7, height = 5, dpi = 300, bg = "white")

# ---- V9 leave-one-donor-out ----------------------------------------------------
loo <- do.call(rbind, lapply(rownames(meta), function(s) {
  k <- setdiff(rownames(meta), s); r <- run(counts[, k], meta[k, ])
  data.frame(left_out = s, condition = meta[s, "condition"], n_nuclei = meta[s, "n_cells"], neuronal_frac = round(meta[s, "neuronal_frac"], 3),
             n_DEGs = nsig(r), r_lfc_vs_full = cor(full$lfc_mle, r$lfc_mle[match(full$gene, r$gene)], use = "complete.obs"))
}))
loo$pct_change_DEGs <- round(100 * (loo$n_DEGs / nsig(full) - 1), 1)
write.csv(loo, file.path(vdir, "results", "loo_seaad.csv"), row.names = FALSE)
infl <- loo$left_out[which.max(abs(loo$pct_change_DEGs))]
add("V9", "Leave-one-out: no single donor drives the result", sprintf("min r = %.3f; DEG change %+.0f%% to %+.0f%%; most influential %s",
    min(loo$r_lfc_vs_full), min(loo$pct_change_DEGs), max(loo$pct_change_DEGs), infl),
    ifelse(min(loo$r_lfc_vs_full) > 0.9 & max(abs(loo$pct_change_DEGs)) < 50, "PASS", "WARN"))
k <- setdiff(rownames(meta), infl); r <- run(counts[, k], meta[k, ], shrink = TRUE)
r$direction <- ifelse(!is.na(r$padj) & r$padj < PADJ & r$log2FC > LFC, "Up", ifelse(!is.na(r$padj) & r$padj < PADJ & r$log2FC < -LFC, "Down", "NS"))
write.csv(r[, c("gene", "log2FC", "pvalue", "padj", "direction")], file.path(vdir, "results", sprintf("seaad_DEGs_drop_%s.csv", infl)), row.names = FALSE)
writeLines(infl, file.path(vdir, "results", "seaad_most_influential.txt"))

write.csv(checks, file.path(vdir, "results", "checks_seaad.csv"), row.names = FALSE)
sink()
