# =============================================================================
# 03b_figures.R — Dataset 3 (SEA-AD MTG snRNA-seq pseudobulk) figure set
# Reads results/seaad_DEGs.csv, results/seaad_pseudobulk.rds and
# results/seaad_celltype_counts_per_donor.csv (all from 03_pseudobulk_deseq2.R).
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
proj   <- ad_root()
outdir <- file.path(proj, "analysis", "dataset3_seaad"); figdir <- file.path(outdir, "figures")
suppressPackageStartupMessages(library(DESeq2))

deg <- read.csv(file.path(outdir, "results", "seaad_DEGs.csv"))
pb  <- readRDS(file.path(outdir, "results", "seaad_pseudobulk.rds"))
cond <- setNames(pb$coldata$condition, pb$coldata$donor_id)
counts <- pb$counts[deg$gene, names(cond)]
dds <- DESeqDataSetFromMatrix(counts, data.frame(condition = cond, row.names = names(cond)), ~ condition)
dds <- estimateSizeFactors(dds)
meanexpr <- log10(rowMeans(counts(dds, normalized = TRUE)) + 1)
mat <- assay(vst(dds, blind = TRUE))

ds_title  <- "SEA-AD MTG snRNA-seq pseudobulk"
ds_sub    <- "raw UMIs summed per donor; DESeq2 + apeglm; 10 dementia / 10 no dementia"
ds_method <- "DESeq2 with apeglm-shrunk fold changes on donor-level sums of raw UMI counts, 10 dementia vs 10 no-dementia donors, middle temporal gyrus"
xlab_ma   <- "log10 mean normalised count (baseMean + 1)"
ma_desc   <- "log10 mean DESeq2-normalised pseudobulk count"
heat_desc <- "variance-stabilised (VST) pseudobulk expression"
extra_legends <- "Figure 5. Cell-type composition of each donor's single-nucleus sample, shown as the percentage of nuclei assigned to SEA-AD subclasses (grouped into nine classes), with donors ordered control (left) then dementia (right) and total nuclei shown above each bar."
# =============================================================================
# SHARED PLOTTING BLOCK (identical in all three figure scripts)
# Colour scheme fixed across datasets: AD #E63946, control #457B9D,
# Up #E63946, Down #457B9D, NS #AAAAAA. 300 dpi, white background.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(pheatmap) })
COL_COND <- c(control = "#457B9D", AD = "#E63946")
COL_DIR  <- c(Up = "#E63946", Down = "#457B9D", NS = "#AAAAAA")
PADJ <- 0.05; LFC <- 0.5
PANEL33 <- c("VGF","SST","SCG2","CARTPT","CHGB","VSNL1","TAC1","SNAP25","NPTX2","NEFL","GABRA1","NRGN",
             "GABRG2","SYT4","CCK","NEFM","SHANK2","GRIN2A","NECAB1","LRRTM1","CBLN2","CAMK2A","SYN2",
             "HOMER1","NPY","SCN1A","PENK","NLGN1","CADM1","VIM","TREM2","GJA1","GFAP")
theme_fig <- theme_bw(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "top", plot.title = element_text(face = "bold"))
dir_labels <- function(d) c(Up = sprintf("Up in AD (%d)", sum(d == "Up")),
                            Down = sprintf("Down in AD (%d)", sum(d == "Down")), NS = "Not significant")
save_fig <- function(p, name, w, h) ggsave(file.path(figdir, name), p, width = w, height = h, dpi = 300, bg = "white")

make_figures <- function() {
  deg <- deg[!is.na(deg$padj), ]
  deg$direction <- factor(deg$direction, levels = c("Up", "Down", "NS"))
  top10 <- head(deg[order(deg$padj, deg$pvalue), "gene"], 10)

  # ---- fig1: volcano, top 10 genes by padj labelled -------------------------
  deg$label <- ifelse(deg$gene %in% top10, deg$gene, NA)
  p1 <- ggplot(deg, aes(log2FC, -log10(padj))) +
    geom_point(aes(colour = direction), size = 0.9, alpha = 0.6) +
    scale_colour_manual(values = COL_DIR, labels = dir_labels(deg$direction), drop = FALSE) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(PADJ), linetype = "dashed", colour = "grey50") +
    geom_text_repel(aes(label = label), size = 3.2, max.overlaps = Inf, na.rm = TRUE, seed = 42,
                    box.padding = 0.5, min.segment.length = 0) +
    labs(title = paste(ds_title, "— volcano"), subtitle = ds_sub,
         x = "log2 fold change (AD − control)", y = "-log10(adjusted p)", colour = NULL) + theme_fig
  save_fig(p1, "fig1_volcano.png", 8, 6.5)

  # ---- fig2: MA plot -----------------------------------------------------------
  deg$A <- meanexpr[deg$gene]
  p2 <- ggplot(deg[order(deg$direction == "NS", decreasing = TRUE), ], aes(A, log2FC)) +
    geom_point(aes(colour = direction), size = 0.8, alpha = 0.6) +
    scale_colour_manual(values = COL_DIR, labels = dir_labels(deg$direction), drop = FALSE) +
    geom_hline(yintercept = 0, colour = "grey30") +
    geom_hline(yintercept = c(-LFC, LFC), linetype = "dashed", colour = "grey50") +
    geom_text_repel(aes(label = label), size = 3.2, max.overlaps = Inf, na.rm = TRUE, seed = 42,
                    min.segment.length = 0) +
    labs(title = paste(ds_title, "— MA plot"), subtitle = ds_sub, x = xlab_ma,
         y = "log2 fold change (AD − control)", colour = NULL) + theme_fig
  save_fig(p2, "fig2_ma_plot.png", 8, 6)

  # ---- fig3: top 10 up + top 10 down by log2FC among padj < 0.05 genes --------
  sig <- deg[deg$padj < PADJ, ]
  up <- head(sig[sig$log2FC > 0, ][order(-sig$log2FC[sig$log2FC > 0]), ], 10)
  dn <- head(sig[sig$log2FC < 0, ][order(sig$log2FC[sig$log2FC < 0]), ], 10)
  bar <- rbind(up, dn); bar$gene <- factor(bar$gene, levels = bar$gene[order(bar$log2FC)])
  bar$dir2 <- ifelse(bar$log2FC > 0, "Up", "Down")
  p3 <- ggplot(bar, aes(log2FC, gene, fill = dir2)) +
    geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = "grey30") +
    scale_fill_manual(values = COL_DIR[c("Up", "Down")], labels = c(Up = "Up in AD", Down = "Down in AD")) +
    labs(title = paste(ds_title, "— largest changes"),
         subtitle = sprintf("Top %d up / %d down by log2FC among genes with padj < 0.05", nrow(up), nrow(dn)),
         x = "log2 fold change (AD − control)", y = NULL, fill = NULL) + theme_fig
  save_fig(p3, "fig3_top20_barplot.png", 7, 7)

  # ---- fig4: 33-gene panel heatmap (z-score per gene, rows clustered) ---------
  g <- intersect(PANEL33, rownames(mat))
  missing <- setdiff(PANEL33, g)
  m <- mat[g, names(cond)[order(cond, names(cond))], drop = FALSE]
  m <- m[apply(m, 1, sd) > 0, , drop = FALSE]
  z <- t(scale(t(m)))
  pj <- deg$padj[match(rownames(z), deg$gene)]
  rownames(z) <- ifelse(!is.na(pj) & pj < PADJ, paste0(rownames(z), " *"), rownames(z))
  ann <- data.frame(condition = cond[colnames(z)], row.names = colnames(z))
  gap <- sum(cond == "control")
  hm <- pheatmap(pmin(pmax(z, -3), 3), cluster_rows = TRUE, cluster_cols = FALSE,
                 color = colorRampPalette(c("#457B9D", "white", "#E63946"))(101), breaks = seq(-3, 3, length.out = 102),
                 annotation_col = ann, annotation_colors = list(condition = COL_COND),
                 show_colnames = ncol(z) <= 40, fontsize_row = 8, fontsize_col = 7, gaps_col = gap,
                 border_color = NA, silent = TRUE,
                 main = sprintf("%s\n33-gene panel (z-score per gene; * padj < 0.05)%s", ds_title,
                                if (length(missing)) paste0("\nnot measured: ", paste(missing, collapse = ", ")) else ""))
  ggsave(file.path(figdir, "fig4_33gene_heatmap.png"), hm$gtable, width = if (ncol(z) > 40) 10 else 8,
         height = 8.5, dpi = 300, bg = "white")

  # ---- figure legends ---------------------------------------------------------
  n_up <- sum(deg$direction == "Up"); n_dn <- sum(deg$direction == "Down")
  leg <- c(
    sprintf("Figure 1. Volcano plot of differential expression in %s (%s); each point is a gene, coloured red (up, n = %d) or blue (down, n = %d) when padj < 0.05 and |log2FC| > 0.5, with dashed lines marking these thresholds and the 10 genes with the smallest padj labelled.",
            ds_title, ds_method, n_up, n_dn),
    sprintf("Figure 2. MA plot of log2 fold change (AD − control) against %s for all tested genes in %s, coloured by significance as in Figure 1 and with the same 10 genes labelled.", ma_desc, ds_title),
    sprintf("Figure 3. The 10 most up-regulated and 10 most down-regulated genes by log2 fold change among genes with padj < 0.05 in %s (red, higher in AD; blue, lower in AD).", ds_title),
    sprintf("Figure 4. Heatmap of the 33-gene direction-consistent candidate panel in %s, showing per-gene z-scored %s across samples (columns, grouped control then AD as indicated by the top bar) with genes hierarchically clustered (Euclidean distance, complete linkage) and asterisks marking genes with padj < 0.05%s.",
            ds_title, heat_desc, if (length(missing)) paste0("; ", paste(missing, collapse = " and "), " had no measurement in this dataset") else ""))
  if (exists("extra_legends")) leg <- c(leg, extra_legends)
  writeLines(leg, file.path(figdir, "figure_legends.txt"))
  cat("Figures written to", figdir, "\n"); print(list.files(figdir))
}
make_figures()

# ---- fig5: cell-type composition per donor (SEA-AD only) ----------------------
# 24 SEA-AD subclasses grouped into 9 classes so colours stay distinguishable.
cc <- read.csv(file.path(outdir, "results", "seaad_celltype_counts_per_donor.csv"), row.names = 1, check.names = FALSE)
grp <- c("L2/3 IT" = "Excitatory IT", "L4 IT" = "Excitatory IT", "L5 IT" = "Excitatory IT", "L6 IT" = "Excitatory IT",
         "L6 IT Car3" = "Excitatory IT", "L5 ET" = "Excitatory deep (ET/NP/CT/L6b)", "L5/6 NP" = "Excitatory deep (ET/NP/CT/L6b)",
         "L6 CT" = "Excitatory deep (ET/NP/CT/L6b)", "L6b" = "Excitatory deep (ET/NP/CT/L6b)",
         "Pvalb" = "Pvalb / Chandelier", "Chandelier" = "Pvalb / Chandelier", "Sst" = "Sst / Sst Chodl", "Sst Chodl" = "Sst / Sst Chodl",
         "Vip" = "Vip / Lamp5 / Sncg / Pax6", "Lamp5" = "Vip / Lamp5 / Sncg / Pax6", "Lamp5 Lhx6" = "Vip / Lamp5 / Sncg / Pax6",
         "Sncg" = "Vip / Lamp5 / Sncg / Pax6", "Pax6" = "Vip / Lamp5 / Sncg / Pax6", "Astrocyte" = "Astrocyte",
         "Immune" = "Microglia / immune", "Oligodendrocyte" = "Oligodendrocyte / OPC", "OPC" = "Oligodendrocyte / OPC",
         "Endothelial" = "Vascular", "VLMC & Perivascular" = "Vascular")
stopifnot(all(colnames(cc) %in% names(grp)))
lv <- c("Excitatory IT", "Excitatory deep (ET/NP/CT/L6b)", "Pvalb / Chandelier", "Sst / Sst Chodl",
        "Vip / Lamp5 / Sncg / Pax6", "Astrocyte", "Microglia / immune", "Oligodendrocyte / OPC", "Vascular")
long <- do.call(rbind, lapply(rownames(cc), function(d) {
  v <- tapply(unlist(cc[d, ]), grp[colnames(cc)], sum)
  data.frame(donor = d, class = names(v), n = as.numeric(v), pct = 100 * as.numeric(v) / sum(v))
}))
long$class <- factor(long$class, levels = rev(lv))
ord <- names(cond)[order(cond, names(cond))]
long$donor <- factor(long$donor, levels = ord)
tot <- data.frame(donor = factor(ord, levels = ord), n = rowSums(cc)[ord], condition = cond[ord])
pal <- c("#2a78d6", "#6250d6", "#1baf7a", "#eb6834", "#eda100", "#e87ba4", "#008300", "#8b5a2b", "#AAAAAA")
names(pal) <- lv
p5 <- ggplot(long, aes(donor, pct, fill = class)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  geom_text(data = tot, aes(donor, 103, label = format(n, big.mark = ",")), inherit.aes = FALSE, size = 2.6, angle = 90, hjust = 0) +
  geom_tile(data = tot, aes(donor, -2.5, fill = NULL, colour = NULL), inherit.aes = FALSE, height = 2.5,
            fill = COL_COND[as.character(tot$condition)]) +
  annotate("text", x = sum(cond == "control") / 2 + 0.5, y = -8, label = "control (no dementia)", colour = COL_COND[["control"]], fontface = "bold", size = 3.5) +
  annotate("text", x = sum(cond == "control") + sum(cond == "AD") / 2 + 0.5, y = -8, label = "AD (dementia)", colour = COL_COND[["AD"]], fontface = "bold", size = 3.5) +
  geom_vline(xintercept = sum(cond == "control") + 0.5, colour = "grey20", linewidth = 0.6) +
  scale_fill_manual(values = pal, breaks = lv) +
  scale_y_continuous(limits = c(-11, 118), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(title = "SEA-AD MTG — cell-type composition per donor",
       subtitle = "% of nuclei per class; number above bar = total nuclei", x = NULL, y = "% of nuclei", fill = NULL) +
  theme_fig + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), legend.position = "right")
save_fig(p5, "fig5_celltype_composition.png", 11, 6.5)
print(list.files(figdir))
