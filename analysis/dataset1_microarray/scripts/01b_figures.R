# =============================================================================
# 01b_figures.R — Dataset 1 (GSE132903 microarray) figure set
# Reads results/microarray_DEGs.csv (from 01_limma_analysis.R) and the raw GEO
# series matrix (for per-sample expression). Writes figures/fig1-4 + legends.
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
outdir <- file.path(proj, "analysis", "dataset1_microarray"); figdir <- file.path(outdir, "figures")
suppressPackageStartupMessages(library(illuminaHumanv4.db))

# ---- DE results ----
deg <- read.csv(file.path(outdir, "results", "microarray_DEGs.csv"))

# ---- per-sample log2 intensities, same probe rule as script 1 (highest-mean probe per gene) ----
infile <- file.path(proj, "data", "GSE132903", "GSE132903_matrix.txt.gz")
con <- gzfile(infile, "rt"); header <- readLines(con, n = 200); close(con)
fld <- function(tag) gsub('"', '', strsplit(header[grep(paste0("^", tag), header)][1], "\t")[[1]][-1])
titles <- fld("!Sample_title"); gsm <- fld("!Sample_geo_accession")
expr <- as.matrix(read.delim(infile, comment.char = "!", row.names = 1, check.names = FALSE)[, gsm])
sym  <- unname(mapIds(illuminaHumanv4.db, keys = rownames(expr), column = "SYMBOL", keytype = "PROBEID", multiVals = "first"))
ok   <- !is.na(sym) & sym != ""; expr <- expr[ok, ]; sym <- sym[ok]
o    <- order(sym, -rowMeans(expr)); expr <- expr[o, ]; sym <- sym[o]
mat  <- expr[!duplicated(sym), ]; rownames(mat) <- sym[!duplicated(sym)]
cond <- factor(setNames(ifelse(grepl("_AD_", titles), "AD", "control"), gsm), levels = c("control", "AD"))
meanexpr <- rowMeans(mat)

ds_title  <- "GSE132903 microarray (MTG)"
ds_sub    <- "limma; 97 AD / 98 control; positive log2FC = higher in AD"
ds_method <- "limma, 97 AD vs 98 control, middle temporal gyrus"
xlab_ma   <- "Mean log2 intensity (A)"
ma_desc   <- "mean log2 probe intensity"
heat_desc <- "log2 intensities"
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
