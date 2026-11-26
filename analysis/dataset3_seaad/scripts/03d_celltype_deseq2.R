# =============================================================================
# 03d_celltype_deseq2.R — cell-type-resolved DE in SEA-AD + validation
# -----------------------------------------------------------------------------
# For each cell type: DESeq2 on donor-level raw-UMI pseudobulk (dementia vs no
# dementia), apeglm-shrunk log2FC, then the same tests that failed for the
# whole-tissue analysis:
#   - label-permutation null (20 shuffles; empirical p)
#   - dropping the two glia-heavy AD donors (H20.33.028, H20.33.038)
# and concordance of each cell type's fold changes with the two bulk datasets.
# Inputs : results/celltype/seaad_celltype_pseudobulk.rds (script 03c),
#          dataset1/2 DEG tables, dataset3 whole-tissue DEG table.
# Outputs: results/celltype/*, figures/celltype/*, and
#          ../validation/results/checks_seaad_celltype.csv
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
proj <- file.path(ad_root(), "analysis")
d3 <- file.path(proj, "dataset3_seaad"); rdir <- file.path(d3, "results", "celltype"); fdir <- file.path(d3, "figures", "celltype")
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(DESeq2); library(apeglm); library(ggplot2); library(pheatmap); library(parallel) })
sink(file.path(rdir, "03d_run_log.txt"), split = TRUE)

PADJ <- 0.05; LFC <- 0.5; MIN_NUCLEI <- 20; MIN_DONORS <- 4; NPERM <- 20; CORES <- 4
DROP <- c("H20.33.028", "H20.33.038")                        # glia-heavy AD donors flagged in validation V9
CT_LABEL <- c(Exc_IT = "Excitatory IT", Exc_deep = "Excitatory deep", Inh_Pvalb = "Pvalb/Chandelier", Inh_Sst = "Sst/Sst Chodl",
              Inh_CGE = "Vip/Lamp5/Sncg/Pax6", Astrocyte = "Astrocyte", Microglia = "Microglia/immune",
              Oligodendrocyte = "Oligodendrocyte", OPC = "OPC", Vascular = "Vascular")
COL_UP <- "#E63946"; COL_DN <- "#457B9D"
PANEL33 <- c("VGF","SST","SCG2","CARTPT","CHGB","VSNL1","TAC1","SNAP25","NPTX2","NEFL","GABRA1","NRGN","GABRG2","SYT4","CCK","NEFM",
             "SHANK2","GRIN2A","NECAB1","LRRTM1","CBLN2","CAMK2A","SYN2","HOMER1","NPY","SCN1A","PENK","NLGN1","CADM1","VIM","TREM2","GJA1","GFAP")
th <- theme_bw(base_size = 11) + theme(plot.background = element_rect(fill = "white", colour = NA), plot.title = element_text(face = "bold"))

pb <- readRDS(file.path(rdir, "seaad_celltype_pseudobulk.rds")); C0 <- pb$counts; M0 <- pb$coldata

fit <- function(cn, cd, shrink = FALSE) {
  dds <- DESeqDataSetFromMatrix(cn, cd, ~ condition); dds$condition <- relevel(dds$condition, "control")
  dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
  dds <- DESeq(dds, quiet = TRUE)
  r <- results(dds, name = "condition_AD_vs_control", alpha = PADJ)
  o <- data.frame(gene = rownames(r), baseMean = r$baseMean, lfc_mle = r$log2FoldChange, pvalue = r$pvalue, padj = r$padj)
  if (shrink) o$log2FC <- lfcShrink(dds, coef = "condition_AD_vs_control", type = "apeglm", quiet = TRUE)$log2FoldChange
  o
}
nsig <- function(o) sum(!is.na(o$padj) & o$padj < PADJ & abs(o$lfc_mle) > LFC)

# ---- per cell type ------------------------------------------------------------
# Resume: if a previous run already saved the per-cell-type results (the slow part:
# 20 permutations + drop-2 refit per cell type), reuse them. Delete
# results/celltype/celltype_summary.csv to force a full recomputation.
summ <- list(); res <- list()
sumf <- file.path(rdir, "celltype_summary.csv")
if (file.exists(sumf)) {
  S <- read.csv(sumf)
  for (ct in S$celltype) res[[ct]] <- read.csv(file.path(rdir, sprintf("DEGs_%s.csv", ct)))
  cat("Reusing saved per-cell-type results:\n"); print(S[, c("celltype", "DEGs", "perm_p", "drop2_retained")], row.names = FALSE)
}
if (!file.exists(sumf)) for (ct in names(CT_LABEL)) {
  cd <- M0[M0$celltype == ct & M0$n_nuclei >= MIN_NUCLEI, ]
  nA <- sum(cd$condition == "AD"); nC <- sum(cd$condition == "control")
  if (nA < MIN_DONORS || nC < MIN_DONORS) {
    cat(sprintf("%-16s skipped: %d AD / %d control donors with >= %d nuclei\n", ct, nA, nC, MIN_NUCLEI)); next }
  cn <- C0[, cd$sample]
  full <- fit(cn, cd, shrink = TRUE)
  full$direction <- ifelse(!is.na(full$padj) & full$padj < PADJ & full$log2FC > LFC, "Up",
                    ifelse(!is.na(full$padj) & full$padj < PADJ & full$log2FC < -LFC, "Down", "NS"))
  write.csv(full[order(full$padj), c("gene", "log2FC", "lfc_mle", "baseMean", "pvalue", "padj", "direction")],
            file.path(rdir, sprintf("DEGs_%s.csv", ct)), row.names = FALSE)
  obs <- nsig(full)
  # permutation null: shuffle condition among this cell type's donors
  nul <- unlist(mclapply(seq_len(NPERM), function(i) { set.seed(1000 + i); cp <- cd; cp$condition <- sample(cp$condition); nsig(fit(cn, cp)) },
                         mc.cores = CORES))
  pe <- (sum(nul >= obs) + 1) / (NPERM + 1)
  # drop the two glia-heavy AD donors
  k <- !(cd$donor_id %in% DROP); dr <- fit(cn[, k], cd[k, ])
  keep_frac <- if (obs > 0) length(intersect(full$gene[!is.na(full$padj) & full$padj < PADJ & abs(full$lfc_mle) > LFC],
                                              dr$gene[!is.na(dr$padj) & dr$padj < PADJ & abs(dr$lfc_mle) > LFC])) / obs else NA
  r_drop <- cor(full$lfc_mle, dr$lfc_mle[match(full$gene, dr$gene)], use = "complete.obs")
  summ[[ct]] <- data.frame(celltype = ct, label = CT_LABEL[[ct]], donors_AD = nA, donors_control = nC,
                           median_nuclei = median(cd$n_nuclei), genes_tested = nrow(full),
                           DEGs = sum(full$direction != "NS"), up = sum(full$direction == "Up"), down = sum(full$direction == "Down"),
                           DEGs_mle = obs, perm_null_median = median(nul), perm_null_max = max(nul), perm_p = round(pe, 3),
                           drop2_r_lfc = round(r_drop, 3), drop2_DEGs = nsig(dr), drop2_retained = round(keep_frac, 2))
  res[[ct]] <- full
  cat(sprintf("%-16s %2d AD/%2d ctrl  genes %5d  DEGs %4d (up %4d / down %4d)  perm p = %.3f (null max %d)  drop-2 r = %.2f, retained %.0f%%\n",
              ct, nA, nC, nrow(full), sum(full$direction != "NS"), sum(full$direction == "Up"), sum(full$direction == "Down"),
              pe, max(nul), r_drop, 100 * keep_frac))
}
if (length(summ)) { S <- do.call(rbind, summ); write.csv(S, sumf, row.names = FALSE) }

# ---- concordance with the bulk datasets ----------------------------------------
rd <- function(f) { d <- read.csv(f); d[!duplicated(d$gene) & !is.na(d$padj), ] }
B <- list(Microarray = rd(file.path(proj, "dataset1_microarray/results/microarray_DEGs.csv")),
          `RNA-seq` = rd(file.path(proj, "dataset2_rnaseq/results/rnaseq_DEGs.csv")))
wt <- rd(file.path(d3, "results", "seaad_DEGs.csv"))
conc <- do.call(rbind, lapply(names(c(list(Whole_tissue = wt), res)), function(ct) {
  x <- if (ct == "Whole_tissue") wt else res[[ct]]
  x <- x[!is.na(x$padj), ]; x$sig <- x$direction != "NS"
  do.call(rbind, lapply(names(B), function(b) {
    y <- B[[b]]; m <- merge(x[, c("gene", "log2FC", "sig")], data.frame(gene = y$gene, lf = y$log2FC, s2 = y$direction != "NS"), by = "gene")
    both <- m$sig & m$s2
    data.frame(celltype = ct, bulk = b, n_genes = nrow(m), pearson = round(cor(m$log2FC, m$lf), 3),
               spearman = round(cor(m$log2FC, m$lf, method = "spearman"), 3), sig_both = sum(both),
               same_direction = sum(both & sign(m$log2FC) == sign(m$lf)),
               pct_same = if (sum(both)) round(100 * mean(sign(m$log2FC[both]) == sign(m$lf[both])), 1) else NA)
  }))
}))
write.csv(conc, file.path(rdir, "concordance_with_bulk.csv"), row.names = FALSE)

# ---- 33-gene panel by cell type -------------------------------------------------
lfm <- sapply(res, function(x) x$log2FC[match(PANEL33, x$gene)]); rownames(lfm) <- PANEL33
pjm <- sapply(res, function(x) x$padj[match(PANEL33, x$gene)]); rownames(pjm) <- PANEL33
mle <- sapply(res, function(x) x$lfc_mle[match(PANEL33, x$gene)]); rownames(mle) <- PANEL33
bm  <- sapply(res, function(x) x$baseMean[match(PANEL33, x$gene)]); rownames(bm) <- PANEL33
home <- apply(bm, 1, function(v) if (all(is.na(v))) NA else names(which.max(v)))   # cell type with highest expression
P <- data.frame(gene = PANEL33, home_celltype = home,
                home_log2FC = mapply(function(g, h) if (is.na(h)) NA else lfm[g, h], PANEL33, home),
                home_log2FC_mle = mapply(function(g, h) if (is.na(h)) NA else mle[g, h], PANEL33, home),
                home_padj = mapply(function(g, h) if (is.na(h)) NA else pjm[g, h], PANEL33, home),
                whole_tissue_log2FC_mle = NA)
wtm <- fit(readRDS(file.path(d3, "results", "seaad_pseudobulk.rds"))$counts, readRDS(file.path(d3, "results", "seaad_pseudobulk.rds"))$coldata)
P$whole_tissue_log2FC_mle <- wtm$lfc_mle[match(P$gene, wtm$gene)]
P$n_celltypes_sig <- rowSums(!is.na(pjm) & pjm < PADJ)
write.csv(P, file.path(rdir, "panel33_by_celltype.csv"), row.names = FALSE)
write.csv(cbind(gene = PANEL33, as.data.frame(round(lfm, 3))), file.path(rdir, "panel33_log2FC_matrix.csv"), row.names = FALSE)

# ---- validation checks (feed the validation report) ----------------------------
chk <- do.call(rbind, lapply(seq_len(nrow(S)), function(i) { s <- S[i, ]
  v <- if (s$DEGs_mle == 0) "INFO" else if (s$perm_p > 0.05) "FAIL" else if (s$perm_null_max >= 0.5 * s$DEGs_mle) "WARN" else "PASS"
  data.frame(dataset = "SEA-AD cell-type", id = paste0("T", i), check = paste0("Permutation null + drop-2 donors: ", s$label),
             value = sprintf("%d DEGs (MLE %d); null median %.0f, max %d; perm p = %.3f; drop H20.33.028/038: r = %.2f, retained %s",
                             s$DEGs, s$DEGs_mle, s$perm_null_median, s$perm_null_max, s$perm_p, s$drop2_r_lfc,
                             ifelse(is.na(s$drop2_retained), "n/a", sprintf("%.0f%%", 100 * s$drop2_retained))),
             verdict = v, note = sprintf("%d AD / %d control donors", s$donors_AD, s$donors_control)) }))
cm <- conc[conc$bulk == "Microarray", ]
chk <- rbind(chk, data.frame(dataset = "SEA-AD cell-type", id = "T_conc", check = "Agreement with microarray, by cell type (Pearson r)",
                             value = paste(sprintf("%s %.2f", cm$celltype, cm$pearson), collapse = "; "), verdict = "INFO",
                             note = "which cell type's changes best match bulk tissue"))
write.csv(chk, file.path(proj, "validation", "results", "checks_seaad_celltype.csv"), row.names = FALSE)

# =============================== FIGURES =========================================
lab <- function(x) CT_LABEL[x]
# F1 DEGs per cell type with permutation p
d1 <- rbind(data.frame(ct = S$celltype, dir = "Up", n = S$up), data.frame(ct = S$celltype, dir = "Down", n = -S$down))
d1$ct <- factor(lab(d1$ct), levels = rev(lab(S$celltype)))
pl <- data.frame(ct = factor(lab(S$celltype), levels = levels(d1$ct)), y = S$up,
                 txt = sprintf("perm p = %.3f%s", S$perm_p, ifelse(S$perm_p <= 0.05, " (best possible, 0/20 shuffles)", "  (not above chance)")))
g1 <- ggplot(d1, aes(n, ct, fill = dir)) + geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = "grey30") +
  geom_text(data = pl, aes(x = pmax(y, 0) + max(abs(d1$n)) * 0.03, y = ct, label = txt), inherit.aes = FALSE, hjust = 0, size = 3.2) +
  scale_fill_manual(values = c(Up = COL_UP, Down = COL_DN), labels = c(Up = "Up in AD", Down = "Down in AD")) +
  scale_x_continuous(labels = abs, expand = expansion(mult = c(0.05, 0.45))) +
  labs(title = "SEA-AD cell-type-resolved DE: dementia vs no dementia", subtitle = "Genes with padj < 0.05 and |log2FC| > 0.5 (apeglm); label = label-permutation p (20 shuffles)",
       x = "number of genes (down ← | → up)", y = NULL, fill = NULL) + th + theme(legend.position = "top")
ggsave(file.path(fdir, "fig1_DEGs_per_celltype.png"), g1, width = 9, height = 5.5, dpi = 300, bg = "white")

# F2 panel heatmap by cell type
lim <- 1.5; z <- pmin(pmax(lfm, -lim), lim); z[is.na(z)] <- 0
st <- matrix(ifelse(!is.na(pjm) & pjm < PADJ, "*", ""), nrow(pjm)); st[is.na(lfm)] <- "·"
colnames(z) <- lab(colnames(z))
hm <- pheatmap(z, cluster_rows = TRUE, cluster_cols = FALSE, display_numbers = st, number_color = "black", fontsize_number = 10,
               color = colorRampPalette(c(COL_DN, "white", COL_UP))(101), breaks = seq(-lim, lim, length.out = 102),
               border_color = "white", angle_col = 45, fontsize_row = 8, cellwidth = 34, silent = TRUE,
               main = "33-gene panel within each cell type\nlog2FC (apeglm, dementia − control); * padj < 0.05; · not tested")
ggsave(file.path(fdir, "fig2_panel33_by_celltype.png"), hm$gtable, width = 8, height = 10, dpi = 300, bg = "white")

# F3 concordance with bulk datasets
conc$ct <- factor(ifelse(conc$celltype == "Whole_tissue", "Whole tissue (script 03)", lab(conc$celltype)),
                  levels = rev(c("Whole tissue (script 03)", lab(S$celltype))))
g3 <- ggplot(conc, aes(spearman, ct, fill = bulk)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey30") + scale_fill_manual(values = c(Microarray = "#2196F3", `RNA-seq` = "#FF9800")) +
  labs(title = "Which SEA-AD cell type's changes match bulk tissue?", subtitle = "Spearman correlation of log2FC with each bulk dataset (all shared genes)",
       x = "Spearman rho", y = NULL, fill = "Bulk dataset") + th + theme(legend.position = "top")
ggsave(file.path(fdir, "fig3_concordance_with_bulk.png"), g3, width = 8, height = 5.5, dpi = 300, bg = "white")

# F4 whole tissue vs within home cell type (unshrunk MLE for both, like with like)
keyg <- c("GFAP", "GJA1", "VIM", "TREM2", "SST", "NPY", "TAC1", "VGF", "NPTX2", "SNAP25", "CHGB", "CARTPT")
d4 <- P[P$gene %in% keyg, ]
d4 <- rbind(data.frame(gene = d4$gene, what = "Whole tissue", lfc = d4$whole_tissue_log2FC_mle, home = d4$home_celltype),
            data.frame(gene = d4$gene, what = "Within its main cell type", lfc = d4$home_log2FC_mle, home = d4$home_celltype))
d4$glab <- factor(sprintf("%s (%s)", d4$gene, lab(d4$home)), levels = unique(sprintf("%s (%s)", d4$gene, lab(d4$home))[order(match(d4$gene, keyg))]))
g4 <- ggplot(d4, aes(lfc, glab, fill = what)) + geom_col(position = position_dodge(0.75), width = 0.7) + geom_vline(xintercept = 0, colour = "grey30") +
  scale_fill_manual(values = c("Whole tissue" = "#AAAAAA", "Within its main cell type" = "#2a78d6")) +
  scale_y_discrete(limits = rev) +
  labs(title = "Whole-tissue change vs change inside the expressing cell type", subtitle = "SEA-AD, unshrunk log2FC (dementia − control)",
       x = "log2 fold change", y = NULL, fill = NULL) + th + theme(legend.position = "top")
ggsave(file.path(fdir, "fig4_tissue_vs_celltype.png"), g4, width = 8.5, height = 6, dpi = 300, bg = "white")

cat("\n==== SUMMARY ====\n"); print(S[, c("label", "donors_AD", "donors_control", "genes_tested", "DEGs", "up", "down", "perm_p", "perm_null_max", "drop2_r_lfc", "drop2_retained")], row.names = FALSE)
cat("\n==== CONCORDANCE WITH BULK ====\n"); print(conc[, 1:8], row.names = FALSE)
cat("\n==== 33-GENE PANEL IN HOME CELL TYPE ====\n"); print(P, row.names = FALSE, digits = 3)
print(sessioninfo::session_info())
sink()
