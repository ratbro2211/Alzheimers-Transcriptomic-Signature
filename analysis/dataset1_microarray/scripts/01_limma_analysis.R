# =============================================================================
# 01_limma_analysis.R  —  Dataset 1: GSE132903 microarray
# -----------------------------------------------------------------------------
# Data   : GSE132903, Illumina HumanHT-12 v4, middle temporal gyrus (MTG)
#          97 Alzheimer's disease (AD) + 98 non-demented controls (ND)
# Model  : limma (lmFit + eBayes) on log2 intensities
# Contrast: AD − ND  (positive log2FC = higher in AD)
#
# Why limma and not limma-voom:
#   voom models the mean-variance trend of COUNT data (RNA-seq). The GEO series
#   matrix already contains normalised log2 intensities, for which the standard
#   limma linear model + empirical-Bayes moderation is the correct method.
#   Applying voom to log-intensities would be a modelling error.
#
# Self-contained: reads only the raw GEO series matrix, writes only into
# analysis/dataset1_microarray/. Run from any working directory.
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

set.seed(42)   # no randomness below, set for reproducibility per project rules

suppressPackageStartupMessages({
  library(limma)               # linear models + empirical Bayes for microarrays
  library(illuminaHumanv4.db)  # Illumina HT-12 v4 probe -> gene symbol annotation
  library(ggplot2)             # volcano plot
  library(ggrepel)             # non-overlapping gene labels
})

# ---- Paths -------------------------------------------------------------------
proj   <- ad_root()
infile <- file.path(proj, "data", "GSE132903", "GSE132903_matrix.txt.gz")
outdir <- file.path(proj, "analysis", "dataset1_microarray")
stopifnot(file.exists(infile))

# ---- Step 1: sample labels from the series-matrix header ---------------------
# Header lines start with "!". Sample titles look like
# "middle_temporal_gyrus_AD_96_02-41" / "..._ND_...", which encode the diagnosis.
con <- gzfile(infile, "rt"); header <- readLines(con, n = 200); close(con)
get_field <- function(tag) gsub('"', '', strsplit(header[grep(paste0("^", tag), header)][1], "\t")[[1]][-1])
titles <- get_field("!Sample_title")
gsm    <- get_field("!Sample_geo_accession")
stopifnot(length(titles) == length(gsm))
diag_raw <- ifelse(grepl("_AD_", titles), "AD", ifelse(grepl("_ND_", titles), "ND", NA))
stopifnot(!anyNA(diag_raw))                       # every sample must be labelled AD or ND
names(diag_raw) <- gsm

# ---- Step 2: expression matrix ----------------------------------------------
# The table between !series_matrix_table_begin/end: probes x samples (log2 intensities).
expr <- read.delim(infile, comment.char = "!", row.names = 1, check.names = FALSE)
expr <- as.matrix(expr[, gsm])                    # align columns to the label vector
cat(sprintf("Expression matrix: %d probes x %d samples; value range %.1f-%.1f (log2 scale)\n",
            nrow(expr), ncol(expr), min(expr, na.rm = TRUE), max(expr, na.rm = TRUE)))

# ---- Step 3: explicit factor levels (THE critical line) ----------------------
# ND is the reference level, so the AD coefficient = mean(AD) − mean(ND).
condition <- factor(diag_raw[colnames(expr)], levels = c("ND", "AD"))
print(table(condition))

# ---- Step 4: limma model -----------------------------------------------------
design <- model.matrix(~ condition)               # columns: (Intercept), conditionAD
stopifnot(identical(colnames(design), c("(Intercept)", "conditionAD")))  # guard against level flips
fit <- eBayes(lmFit(expr, design))
tt  <- topTable(fit, coef = "conditionAD", number = Inf, sort.by = "none")

# ---- Step 5: probe -> gene, one probe per gene --------------------------------
# Multiple probes can map to one gene. We keep the probe with the HIGHEST MEAN
# EXPRESSION — a choice made without looking at the AD/ND result, so it does not
# bias p-values (picking the most significant probe would).
tt$probe  <- rownames(tt)
tt$gene   <- unname(mapIds(illuminaHumanv4.db, keys = tt$probe, column = "SYMBOL",
                           keytype = "PROBEID", multiVals = "first"))
tt <- tt[!is.na(tt$gene) & tt$gene != "", ]
tt <- tt[order(tt$gene, -tt$AveExpr), ]
tt <- tt[!duplicated(tt$gene), ]
# Benjamini-Hochberg FDR recomputed on the gene-level table
tt$padj <- p.adjust(tt$P.Value, method = "BH")

# ---- Step 6: direction call --------------------------------------------------
PADJ <- 0.05; LFC <- 0.5
tt$direction <- ifelse(tt$padj < PADJ & tt$logFC >  LFC, "Up",
                ifelse(tt$padj < PADJ & tt$logFC < -LFC, "Down", "NS"))

res <- data.frame(gene = tt$gene, log2FC = tt$logFC, pvalue = tt$P.Value,
                  padj = tt$padj, direction = tt$direction)
res <- res[order(res$padj, res$pvalue), ]
write.csv(res, file.path(outdir, "results", "microarray_DEGs.csv"), row.names = FALSE)
cat(sprintf("\nGenes tested: %d | significant (padj<%.2f, |log2FC|>%.1f): %d (Up %d / Down %d)\n",
            nrow(res), PADJ, LFC, sum(res$direction != "NS"),
            sum(res$direction == "Up"), sum(res$direction == "Down")))

# ---- Step 7: sanity check on known AD markers --------------------------------
# Expected from the literature: GFAP, CD44, SERPINA3 up; VGF, SST, NPTX2, BDNF down.
markers <- c("GFAP", "CD44", "SERPINA3", "VGF", "SST", "NPTX2", "BDNF")
cat("\nSanity check — canonical AD markers:\n")
print(res[res$gene %in% markers, c("gene", "log2FC", "padj", "direction")], row.names = FALSE)

# ---- Step 8: volcano plot ----------------------------------------------------
res$label <- ifelse(res$gene %in% c(head(res$gene[res$direction != "NS"], 20), markers), res$gene, NA)
p <- ggplot(res, aes(log2FC, -log10(padj))) +
  geom_point(aes(colour = direction), size = 0.9, alpha = 0.6) +
  scale_colour_manual(values = c(Up = "#c42b2b", Down = "#1f5fb0", NS = "grey75"),
                      labels = c(Up = sprintf("Up in AD (%d)", sum(res$direction == "Up")),
                                 Down = sprintf("Down in AD (%d)", sum(res$direction == "Down")),
                                 NS = "Not significant")) +
  geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(PADJ), linetype = "dashed", colour = "grey50") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 40, na.rm = TRUE, seed = 42) +
  labs(title = "GSE132903 microarray — AD vs control (MTG)",
       subtitle = "limma; 97 AD / 98 ND; positive log2FC = higher in AD",
       x = "log2 fold change (AD − ND)", y = "-log10(adjusted p)", colour = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(outdir, "figures", "fig1_volcano.png"), p, width = 8, height = 6.5, dpi = 200)
cat("\nSaved: results/microarray_DEGs.csv, figures/volcano.png\n\n")

# ---- Step 9: session info ----------------------------------------------------
sessioninfo::session_info()
