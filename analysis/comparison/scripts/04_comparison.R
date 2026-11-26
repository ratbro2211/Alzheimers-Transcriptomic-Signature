# =============================================================================
# 04_comparison.R  —  Cross-platform comparison of the three AD datasets
# -----------------------------------------------------------------------------
# Inputs (all log2FC = AD − control, produced by scripts 01–03):
#   dataset1_microarray/results/microarray_DEGs.csv   GSE132903, limma (unshrunk log2FC)
#   dataset2_rnaseq/results/rnaseq_DEGs.csv           GSE104704, DESeq2 + apeglm
#   dataset3_seaad/results/seaad_DEGs.csv             SEA-AD pseudobulk, DESeq2 + apeglm
# "Significant" everywhere = padj < 0.05 AND |log2FC| > 0.5 (the `direction` column),
# except the heatmap asterisk, which marks padj < 0.05 as specified.
# Outputs go to comparison/results and comparison/figures only.
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

set.seed(42)
options(width = 200)   # print tables on one line in the log and summary
proj   <- file.path(ad_root(), "analysis")
outdir <- file.path(proj, "comparison"); resdir <- file.path(outdir, "results"); figdir <- file.path(outdir, "figures")
sink(file.path(resdir, "04_run_log.txt"), split = TRUE)

if (!requireNamespace("UpSetR", quietly = TRUE))
  install.packages("UpSetR", repos = "https://cloud.r-project.org", quiet = TRUE)
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(pheatmap); library(UpSetR) })

DS   <- c(micro = "Microarray", rnaseq = "RNA-seq", seaad = "SEA-AD")
DCOL <- c(Microarray = "#2196F3", "RNA-seq" = "#FF9800", "SEA-AD" = "#9C27B0")
COL_UP <- "#E63946"; COL_DN <- "#457B9D"; COL_NS <- "#BBBBBB"
PANEL33 <- c("VGF","SST","SCG2","CARTPT","CHGB","VSNL1","TAC1","SNAP25","NPTX2","NEFL","GABRA1","NRGN",
             "GABRG2","SYT4","CCK","NEFM","SHANK2","GRIN2A","NECAB1","LRRTM1","CBLN2","CAMK2A","SYN2",
             "HOMER1","NPY","SCN1A","PENK","NLGN1","CADM1","VIM","TREM2","GJA1","GFAP")
GLIAL <- c("GFAP", "VIM", "GJA1", "TREM2")
theme_fig <- theme_bw(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA), plot.title = element_text(face = "bold"))

# =============================== PART A =======================================
# ---- 1. load and merge (inner join on gene symbol) ---------------------------
rd <- function(f) { d <- read.csv(file.path(proj, f)); d[!duplicated(d$gene), c("gene", "log2FC", "padj", "direction")] }
L <- list(micro  = rd("dataset1_microarray/results/microarray_DEGs.csv"),
          rnaseq = rd("dataset2_rnaseq/results/rnaseq_DEGs.csv"),
          seaad  = rd("dataset3_seaad/results/seaad_DEGs.csv"))
for (k in names(L)) cat(sprintf("%-10s %6d genes tested, %4d significant\n", DS[k], nrow(L[[k]]), sum(L[[k]]$direction != "NS")))
for (k in names(L)) names(L[[k]])[-1] <- paste0(names(L[[k]])[-1], "_", k)
M <- Reduce(function(a, b) merge(a, b, by = "gene"), L)
M <- M[complete.cases(M[, grep("^log2FC_", names(M))]), ]
cat(sprintf("\nGenes present in all three datasets (inner join): %d\n", nrow(M)))
sig <- function(k) !is.na(M[[paste0("padj_", k)]]) & M[[paste0("direction_", k)]] != "NS"

# ---- 2. pairwise statistics ----------------------------------------------------
pairs <- list(c("micro", "rnaseq"), c("micro", "seaad"), c("rnaseq", "seaad"))
PW <- do.call(rbind, lapply(pairs, function(p) {
  a <- M[[paste0("log2FC_", p[1])]]; b <- M[[paste0("log2FC_", p[2])]]
  both <- sig(p[1]) & sig(p[2])
  same <- sum(sign(a[both]) == sign(b[both]))
  data.frame(pair = paste(DS[p[1]], "vs", DS[p[2]]),
             r_all_genes = round(cor(a, b), 3),
             r_sig_both = if (sum(both) > 2) round(cor(a[both], b[both]), 3) else NA,
             n_sig_both = sum(both), n_same_direction = same,
             pct_same_direction = if (sum(both)) round(100 * same / sum(both), 1) else NA)
}))
write.csv(PW, file.path(resdir, "pairwise_correlations.csv"), row.names = FALSE)

# =============================== PART B =======================================
# All 33 genes, looked up in each dataset separately (not restricted to the merge)
get <- function(k, col, g) L[[k]][[paste0(col, "_", k)]][match(g, L[[k]]$gene)]
P <- data.frame(gene = PANEL33,
                log2FC_micro = get("micro", "log2FC", PANEL33),   padj_micro = get("micro", "padj", PANEL33),
                log2FC_rnaseq = get("rnaseq", "log2FC", PANEL33), padj_rnaseq = get("rnaseq", "padj", PANEL33),
                log2FC_seaad = get("seaad", "log2FC", PANEL33),   padj_seaad = get("seaad", "padj", PANEL33))
dirs <- sapply(names(DS), function(k) get(k, "direction", PANEL33))
P$n_datasets_sig <- rowSums(!is.na(dirs) & dirs != "NS")
lf <- as.matrix(P[, c("log2FC_micro", "log2FC_rnaseq", "log2FC_seaad")])
P$n_datasets_measured <- rowSums(!is.na(lf))
# TRUE only if measured in all three and all signs agree; NA if a dataset lacks the gene
P$direction_consistent <- ifelse(P$n_datasets_measured < 3, NA, apply(sign(lf), 1, function(s) length(unique(s)) == 1))
P$consensus <- ifelse(rowMeans(sign(lf), na.rm = TRUE) < 0, "down", "up")
write.csv(P, file.path(resdir, "panel_33gene_summary.csv"), row.names = FALSE)

# =============================== PART C =======================================
save_fig <- function(p, name, w, h) ggsave(file.path(figdir, name), p, width = w, height = h, dpi = 300, bg = "white")

# ---- fig1–3: pairwise scatter plots --------------------------------------------
scatter <- function(p, fname) {
  x <- paste0("log2FC_", p[1]); y <- paste0("log2FC_", p[2])
  d <- data.frame(gene = M$gene, x = M[[x]], y = M[[y]])
  both <- sig(p[1]) & sig(p[2])
  d$cls <- ifelse(!both, "Not significant in both",
           ifelse(sign(d$x) == sign(d$y), ifelse(d$x > 0, "Up in both", "Down in both"), "Opposite direction"))
  d$cls <- factor(d$cls, levels = c("Not significant in both", "Up in both", "Down in both", "Opposite direction"))
  d <- d[order(d$cls != "Not significant in both"), ]
  lab <- d[d$gene %in% PANEL33, ]
  r <- PW[PW$pair == paste(DS[p[1]], "vs", DS[p[2]]), ]
  lim <- range(c(d$x, d$y)); lim <- c(max(lim[1], -4), min(lim[2], 4))
  g <- ggplot(d, aes(x, y)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) + geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
    geom_point(aes(colour = cls), size = ifelse(d$cls == "Not significant in both", 0.6, 1.6),
               alpha = ifelse(d$cls == "Not significant in both", 0.35, 0.9)) +
    geom_point(data = lab, shape = 21, size = 2.3, fill = NA, colour = "black", stroke = 0.6) +
    geom_text_repel(data = lab, aes(label = gene), size = 2.8, max.overlaps = Inf, seed = 42,
                    min.segment.length = 0, segment.colour = "grey40") +
    scale_colour_manual(values = c("Not significant in both" = COL_NS, "Up in both" = COL_UP,
                                   "Down in both" = COL_DN, "Opposite direction" = "black"), drop = FALSE) +
    annotate("label", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 3.4,
             label = sprintf("r (all genes) = %.2f, n = %d\nr (sig in both) = %.2f, n = %d\nsame direction: %d/%d (%.0f%%)",
                             r$r_all_genes, nrow(d), r$r_sig_both, r$n_sig_both, r$n_same_direction, r$n_sig_both, r$pct_same_direction)) +
    coord_cartesian(xlim = lim, ylim = lim) +
    labs(title = sprintf("%s vs %s — log2 fold change (AD − control)", DS[p[1]], DS[p[2]]),
         subtitle = "Circled & labelled: 33-gene panel. Dashed: y = x.",
         x = sprintf("%s log2FC", DS[p[1]]), y = sprintf("%s log2FC", DS[p[2]]), colour = NULL) +
    theme_fig + theme(legend.position = "bottom") +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
  save_fig(g, fname, 8, 8.4)
}
scatter(c("micro", "rnaseq"), "fig1_scatter_micro_vs_rnaseq.png")
scatter(c("micro", "seaad"),  "fig2_scatter_micro_vs_seaad.png")
scatter(c("rnaseq", "seaad"), "fig3_scatter_rnaseq_vs_seaad.png")

# ---- fig4: direction agreement (same = solid, opposite = hatched) -------------
# ggpattern is not installed, so the hatching is drawn as 45-degree line segments
# clipped to the "opposite" rectangle. Bars use one neutral fill because each bar
# is a PAIR of datasets; dataset colours appear in the x-axis labels instead.
B <- data.frame(pair = PW$pair, same = PW$n_same_direction, n = PW$n_sig_both, pct = PW$pct_same_direction)
B$opp <- B$n - B$same
w <- 0.32; xi <- seq_len(nrow(B)); ymax <- max(B$n) * 1.12
asp <- ymax / (nrow(B) + 0.4)                        # y-units per x-unit, keeps lines ~45 degrees on screen
clip_hatch <- function(i) {
  if (B$opp[i] == 0) return(NULL)
  xl <- i - w; xr <- i + w; yb <- B$same[i]; yt <- B$n[i]; sp <- max(B$n) / 45
  cs <- seq(yb - 2 * w * asp, yt, by = sp)          # line: y = c + asp * (x - xl)
  do.call(rbind, lapply(cs, function(c0) {
    x0 <- max(xl, xl + (yb - c0) / asp); x1 <- min(xr, xl + (yt - c0) / asp)
    if (x1 <= x0) return(NULL)
    data.frame(x = x0, xend = x1, y = c0 + asp * (x0 - xl), yend = c0 + asp * (x1 - xl))
  }))
}
hatch <- do.call(rbind, lapply(xi, clip_hatch))
axis_lab <- c(sprintf("<span style='color:%s'>Microarray</span><br>vs<br><span style='color:%s'>RNA-seq</span>", DCOL[1], DCOL[2]),
              sprintf("<span style='color:%s'>Microarray</span><br>vs<br><span style='color:%s'>SEA-AD</span>", DCOL[1], DCOL[3]),
              sprintf("<span style='color:%s'>RNA-seq</span><br>vs<br><span style='color:%s'>SEA-AD</span>", DCOL[2], DCOL[3]))
use_md <- requireNamespace("ggtext", quietly = TRUE)
g4 <- ggplot(B) +
  geom_rect(aes(xmin = xi - w, xmax = xi + w, ymin = 0, ymax = same), fill = "#4D4D4D") +
  geom_rect(aes(xmin = xi - w, xmax = xi + w, ymin = same, ymax = n), fill = "white", colour = "#4D4D4D", linewidth = 0.5) +
  { if (!is.null(hatch)) geom_segment(data = hatch, aes(x = x, xend = xend, y = y, yend = yend), colour = "#4D4D4D", linewidth = 0.45) } +
  geom_text(aes(x = xi, y = same / 2, label = sprintf("%d / %d\nsame direction\n(%.0f%%)", same, n, pct)),
            colour = "white", fontface = "bold", size = 3.8, lineheight = 0.95) +
  geom_text(aes(x = xi, y = n, label = sprintf("%d opposite (%.0f%%)", opp, 100 - pct)), vjust = -0.7, size = 3.4) +
  scale_x_continuous(breaks = xi, labels = if (use_md) axis_lab else gsub(" vs ", "\nvs\n", B$pair)) +
  scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
  labs(title = "Direction agreement of genes significant in both datasets",
       subtitle = "Solid = same direction; hatched = opposite direction (padj < 0.05 and |log2FC| > 0.5 in both)",
       x = NULL, y = "Genes significant in both datasets") + theme_fig +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
if (use_md) g4 <- g4 + theme(axis.text.x = ggtext::element_markdown(size = 11, lineheight = 1.1))
save_fig(g4, "fig4_direction_agreement_bars.png", 7.5, 6)

# ---- fig5: 33-gene heatmap (log2FC, datasets as columns) -----------------------
H <- as.matrix(P[, c("log2FC_micro", "log2FC_rnaseq", "log2FC_seaad")]); rownames(H) <- P$gene
colnames(H) <- c("Microarray\n(MTG, n=195)", "RNA-seq\n(LTL, n=22)", "SEA-AD\n(MTG, n=20)")
S <- as.matrix(P[, c("padj_micro", "padj_rnaseq", "padj_seaad")])
star <- matrix(ifelse(!is.na(S) & S < 0.05, "*", ""), nrow(S)); star[is.na(H)] <- "n/a"
# Cell-type annotation: known glial genes = Glial-up; genes lower in AD in every dataset
# where measured = Neuronal-down; anything else (sign disagrees across datasets) = Mixed.
ann <- ifelse(P$gene %in% GLIAL, "Glial-up",
       ifelse(apply(sign(H), 1, function(s) all(s[!is.na(s)] < 0)), "Neuronal-down", "Mixed"))
annrow <- data.frame(`cell-type group` = factor(ann, levels = c("Neuronal-down", "Glial-up", "Mixed")), row.names = P$gene, check.names = FALSE)
Hc <- H; Hc[is.na(Hc)] <- 0                                      # NA treated as 0 for clustering only
lim5 <- max(abs(H), na.rm = TRUE); lim5 <- min(ceiling(lim5 * 2) / 2, 2.5)
hm <- pheatmap(pmin(pmax(Hc, -lim5), lim5), cluster_rows = TRUE, cluster_cols = FALSE,
               color = colorRampPalette(c(COL_DN, "white", COL_UP))(101), breaks = seq(-lim5, lim5, length.out = 102),
               display_numbers = star, number_color = "black", fontsize_number = 11,
               annotation_row = annrow, annotation_names_row = FALSE,
               annotation_colors = list(`cell-type group` = c("Neuronal-down" = COL_DN, "Glial-up" = COL_UP, "Mixed" = "#AAAAAA")),
               angle_col = 0, fontsize_row = 9, fontsize_col = 9, cellwidth = 70, border_color = "white", silent = TRUE,
               main = "33-gene panel: log2FC (AD − control) in each dataset\n* padj < 0.05; n/a = not measured (shown as 0)")
ggsave(file.path(figdir, "fig5_33gene_heatmap.png"), hm$gtable, width = 7, height = 10, dpi = 300, bg = "white")

# ---- fig6: UpSet plot of significant genes ------------------------------------
sets <- lapply(L, function(d) d$gene[!is.na(d[[grep("^padj_", names(d), value = TRUE)]]) & d[[grep("^direction_", names(d), value = TRUE)]] != "NS"])
names(sets) <- DS[names(sets)]
png(file.path(figdir, "fig6_upset_plot.png"), width = 9, height = 6, units = "in", res = 300, bg = "white")
print(upset(fromList(sets), sets = rev(names(sets)), keep.order = TRUE, order.by = "freq",
            sets.bar.color = DCOL[rev(names(sets))], main.bar.color = "grey25", matrix.color = "grey25",
            text.scale = c(1.5, 1.3, 1.2, 1.1, 1.5, 1.3), mainbar.y.label = "Genes in intersection",
            sets.x.label = "Significant genes"))
dev.off()
# ggsave is not applicable to UpSetR's grid output; set the PNG resolution tag explicitly
invisible(system2("sips", c("-s", "dpiWidth", "300", "-s", "dpiHeight", "300", shQuote(file.path(figdir, "fig6_upset_plot.png"))), stdout = FALSE))
ov3 <- Reduce(intersect, sets)

# =============================== PART D =======================================
pct <- function(i) sprintf("%d/%d (%.0f%%)", PW$n_same_direction[i], PW$n_sig_both[i], PW$pct_same_direction[i])
dc  <- P[!is.na(P$direction_consistent) & P$direction_consistent, ]
ndc <- nrow(dc); nmeas <- sum(P$n_datasets_measured == 3)
res_par <- sprintf(paste(
  "Differential expression results from the three datasets were merged on gene symbol, yielding %d genes measured in all three.",
  "Genome-wide log2 fold changes (AD − control) were positively correlated in every pairwise comparison (microarray vs RNA-seq r = %.2f; microarray vs SEA-AD r = %.2f; RNA-seq vs SEA-AD r = %.2f).",
  "Among genes significant in both datasets of a pair (padj < 0.05, |log2FC| > 0.5), the direction of change agreed for %s of microarray/RNA-seq genes, %s of microarray/SEA-AD genes and %s of RNA-seq/SEA-AD genes.",
  "%d genes were significant in all three datasets.",
  "Of the 33-gene candidate panel, %d of the %d genes measured in all three datasets changed in the same direction in every dataset, with neuronal and neuropeptide genes (e.g. VGF, SST, NPTX2, SNAP25) decreased and glial genes (e.g. GFAP, GJA1, VIM) increased in AD.",
  "Panel genes reached significance most often in the microarray (%d genes), which had the largest sample, and least often in SEA-AD (%d genes).",
  "Overall, the direction of the AD transcriptomic signal was reproducible across microarray, bulk RNA-seq and single-nucleus RNA-seq platforms and across temporal-cortex regions, although effect sizes differed between datasets."),
  nrow(M), PW$r_all_genes[1], PW$r_all_genes[2], PW$r_all_genes[3], pct(1), pct(2), pct(3), length(ov3),
  ndc, nmeas, sum(!is.na(dirs[, "micro"]) & dirs[, "micro"] != "NS"), sum(!is.na(dirs[, "seaad"]) & dirs[, "seaad"] != "NS"))
lim_par <- paste(
  "These comparisons have several limitations.",
  "The bulk RNA-seq dataset sampled the lateral temporal lobe whereas the microarray and SEA-AD datasets sampled the middle temporal gyrus, so regional differences cannot be fully separated from platform differences.",
  "Sample sizes differed tenfold (195 vs 22 vs 20), so absence of significance in the smaller datasets largely reflects power, and RNA-seq/SEA-AD fold changes were apeglm-shrunk while microarray fold changes were not, which compresses their magnitudes.",
  "Finally, the SEA-AD whole-tissue pseudobulk is confounded by cell-type composition: dementia donors contributed fewer nuclei and relatively more glia, so part of the apparent glial up-regulation and neuronal down-regulation reflects differences in which cells were sampled rather than changes in expression per cell.")
fmt <- function(df) capture.output(print(df, row.names = FALSE))
out <- c("CROSS-PLATFORM COMPARISON — SUMMARY", strrep("=", 70), "",
         "1. RESULTS", "", strwrap(res_par, 90), "",
         "2. PAIRWISE CORRELATIONS", "", fmt(PW), "",
         sprintf("3. 33-GENE PANEL — DIRECTION-CONSISTENT GENES (%d of %d measured in all three datasets)", ndc, nmeas), "",
         fmt(transform(dc[, c("gene", "log2FC_micro", "padj_micro", "log2FC_rnaseq", "padj_rnaseq",
                              "log2FC_seaad", "padj_seaad", "n_datasets_sig")],
                       log2FC_micro = round(log2FC_micro, 2), log2FC_rnaseq = round(log2FC_rnaseq, 2),
                       log2FC_seaad = round(log2FC_seaad, 3), padj_micro = signif(padj_micro, 2),
                       padj_rnaseq = signif(padj_rnaseq, 2), padj_seaad = signif(padj_seaad, 2))), "",
         "4. LIMITATIONS", "", strwrap(lim_par, 90))
writeLines(out, file.path(resdir, "comparison_summary.txt"))

cat("\n================ PAIRWISE CORRELATIONS ================\n"); print(PW, row.names = FALSE)
cat(sprintf("\nGenes significant in all three datasets: %d (%s)\n", length(ov3), paste(head(sort(ov3), 40), collapse = ", ")))
cat(sprintf("33-gene panel: measured in all three %d; direction-consistent %d; inconsistent: %s; not measured in all: %s\n",
            nmeas, ndc, paste(P$gene[!is.na(P$direction_consistent) & !P$direction_consistent], collapse = ", "),
            paste(P$gene[P$n_datasets_measured < 3], collapse = ", ")))
cat("\nSaved results/: pairwise_correlations.csv, panel_33gene_summary.csv, comparison_summary.txt, 04_run_log.txt\n")
cat("Saved figures/:", paste(list.files(figdir), collapse = ", "), "\n\n")
print(sessioninfo::session_info())
sink()
