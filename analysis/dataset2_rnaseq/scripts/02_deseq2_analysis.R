# =============================================================================
# 02_deseq2_analysis.R  —  Dataset 2: GSE104704 bulk RNA-seq
# -----------------------------------------------------------------------------
# Data   : GSE104704 (Nativio et al. 2018), lateral temporal lobe, Illumina RNA-seq
#          Groups in GEO: Young (n=8), Old (n=10, age-matched controls), AD (n=12)
#          Analysis: AD vs Old ("control"); Young excluded (not age-matched)
# Model  : DESeq2, design = ~ condition, reference = control
#          log2FC = AD − control, shrunk with apeglm
#
# ABOUT THE INPUT COUNTS (important)
#   GEO provides no per-sample raw count files for this series (Sample_supplementary_file
#   = NONE). The only table, GSE104704_RNA-Seq_Table.txt.gz, contains "DESeq library
#   size-adjusted tag counts" = raw HTSeq count / sample size factor (non-integers).
#   Because every value in a sample is (integer) / (same size factor), dividing a sample
#   by its smallest non-zero value (= 1 read / size factor) recovers the raw integer
#   counts EXACTLY. The script verifies this (every value must land on an integer to
#   within 1e-6) and stops if it does not. DESeq2 is therefore run on raw counts,
#   never on normalised values.
#
# Self-contained: reads GEO files only, writes only to analysis/dataset2_rnaseq/.
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
gdir   <- file.path(proj, "data", "GSE104704")
outdir <- file.path(proj, "analysis", "dataset2_rnaseq")
logf   <- file.path(outdir, "results", "02_run_log.txt")
sink(logf, split = TRUE)                     # everything printed below also goes to the run log

suppressPackageStartupMessages({
  library(DESeq2)    # negative-binomial GLM for RNA-seq counts
  library(apeglm)    # adaptive shrinkage of log2 fold changes
  library(ggplot2); library(ggrepel)
})

# ---- Step 1: load the GEO table ----------------------------------------------
tab <- read.delim(file.path(gdir, "GSE104704_RNA-Seq_Table.txt.gz"), row.names = 1, check.names = FALSE)
stat_cols <- c("Young_vs_Old", "Old_vs_AD", "Young_vs_AD")     # authors' DE p-values, not data
norm <- as.matrix(tab[, setdiff(colnames(tab), stat_cols)])
cat(sprintf("GEO table: %d genes x %d samples\n", nrow(norm), ncol(norm)))
cat(sprintf("Fraction of entries that are integers as supplied: %.3f  -> values are normalised, not raw\n",
            mean(abs(norm - round(norm)) < 1e-6)))

# ---- Step 2: recover raw integer counts --------------------------------------
min_nz <- apply(norm, 2, function(v) min(v[v > 0]))            # = 1 read / size factor
raw    <- sweep(norm, 2, min_nz, "/")
dev    <- max(abs(raw - round(raw)))
cat(sprintf("Max deviation from integer after rescaling: %.2e\n", dev))
stopifnot(dev < 1e-6)                                          # exact recovery or stop
counts <- round(raw); storage.mode(counts) <- "integer"
cat("Implied DESeq size factors (1 / min non-zero):\n"); print(round(1 / min_nz, 3))
cat(sprintf("Recovered library sizes: %.1f - %.1f million reads\n",
            min(colSums(counts)) / 1e6, max(colSums(counts)) / 1e6))

# ---- Step 3: one row per gene symbol -----------------------------------------
# The authors already kept one representative transcript per gene; apply the
# dataset-1 rule anyway (keep the highest-mean row) so the pipeline is uniform.
n_dup <- sum(duplicated(rownames(counts)))
if (n_dup > 0) {
  o <- order(rownames(counts), -rowMeans(counts))
  counts <- counts[o, ][!duplicated(rownames(counts)[o]), ]
}
cat(sprintf("Duplicate gene rows collapsed: %d\n", n_dup))

# ---- Step 4: sample labels from the GEO series matrix ------------------------
# Sample titles look like "21-1A.RNA.Old"; the part after the last "." is the group.
con <- gzfile(file.path(gdir, "GSE104704_series_matrix.txt.gz"), "rt"); hdr <- readLines(con); close(con)
field <- function(tag) gsub('"', '', strsplit(hdr[grep(paste0("^", tag), hdr)][1], "\t")[[1]][-1])
titles <- field("!Sample_title"); gsm <- field("!Sample_geo_accession")
chars  <- hdr[grep("^!Sample_characteristics_ch1", hdr)]
age    <- as.numeric(sub("age \\(years\\): ", "", gsub('"', '', strsplit(chars[grep("age", chars)], "\t")[[1]][-1])))
meta <- data.frame(sample = sub("\\.[^.]+$", "", titles), geo_group = sub(".*\\.", "", titles),
                   gsm = gsm, age = age, stringsAsFactors = FALSE)
meta <- meta[meta$geo_group %in% c("Old", "AD"), ]              # drop Young
meta$condition <- factor(ifelse(meta$geo_group == "AD", "AD", "control"), levels = c("control", "AD"))
stopifnot(all(meta$sample %in% colnames(counts)))
rownames(meta) <- meta$sample
meta <- meta[order(meta$condition, meta$sample), ]
counts <- counts[, rownames(meta)]

cat("\n================ COLDATA (verify before DESeq2) ================\n")
print(meta[, c("sample", "gsm", "geo_group", "age", "condition")], row.names = FALSE)
print(table(meta$condition))
cat("=================================================================\n\n")

# ---- Step 5: DESeq2 ----------------------------------------------------------
dds <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "control")       # explicit reference
keep <- rowSums(counts(dds) >= 10) >= min(table(dds$condition)) # expressed in >= smallest group
dds  <- dds[keep, ]
cat(sprintf("Genes after low-count filter: %d\n", nrow(dds)))
dds <- DESeq(dds, quiet = TRUE)
stopifnot("condition_AD_vs_control" %in% resultsNames(dds))     # guard: AD vs control, not reverse
res <- results(dds, name = "condition_AD_vs_control", alpha = 0.05)
shr <- lfcShrink(dds, coef = "condition_AD_vs_control", type = "apeglm", quiet = TRUE)

# ---- Step 6: results table ---------------------------------------------------
PADJ <- 0.05; LFC <- 0.5
out <- data.frame(gene = rownames(res), log2FC = shr$log2FoldChange,
                  pvalue = res$pvalue, padj = res$padj)
out$direction <- ifelse(!is.na(out$padj) & out$padj < PADJ & out$log2FC >  LFC, "Up",
                 ifelse(!is.na(out$padj) & out$padj < PADJ & out$log2FC < -LFC, "Down", "NS"))
out <- out[order(out$padj, out$pvalue, na.last = TRUE), ]
write.csv(out, file.path(outdir, "results", "rnaseq_DEGs.csv"), row.names = FALSE)
cat(sprintf("Genes tested: %d | padj NA (independent filtering): %d | significant: %d (Up %d / Down %d)\n",
            nrow(out), sum(is.na(out$padj)), sum(out$direction != "NS"),
            sum(out$direction == "Up"), sum(out$direction == "Down")))

cat("\nTop 20 genes by padj:\n")
print(head(out, 20), row.names = FALSE, digits = 3)

# ---- Step 7: validation on canonical AD markers ------------------------------
expect <- c(VGF = "Down", NPTX2 = "Down", SST = "Down", GFAP = "Up",
            BDNF = "Down", CD44 = "Up", SERPINA3 = "Up")
v <- out[match(names(expect), out$gene), c("gene", "log2FC", "padj", "direction")]
v$expected  <- expect
v$sign_ok   <- ifelse(sign(v$log2FC) == ifelse(expect == "Up", 1, -1), "yes", "NO")
cat("\nValidation — canonical AD markers:\n"); print(v, row.names = FALSE, digits = 3)

# ---- Step 8: volcano plot ----------------------------------------------------
pd <- out[!is.na(out$padj), ]
pd$label <- ifelse(pd$gene %in% c(head(pd$gene[pd$direction != "NS"], 20), names(expect)), pd$gene, NA)
p <- ggplot(pd, aes(log2FC, -log10(padj))) +
  geom_point(aes(colour = direction), size = 0.9, alpha = 0.6) +
  scale_colour_manual(values = c(Up = "#c42b2b", Down = "#1f5fb0", NS = "grey75"),
                      labels = c(Up = sprintf("Up in AD (%d)", sum(pd$direction == "Up")),
                                 Down = sprintf("Down in AD (%d)", sum(pd$direction == "Down")),
                                 NS = "Not significant")) +
  geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(PADJ), linetype = "dashed", colour = "grey50") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 40, na.rm = TRUE, seed = 42) +
  labs(title = "GSE104704 bulk RNA-seq — AD vs age-matched control (LTL)",
       subtitle = sprintf("DESeq2 + apeglm; %d AD / %d control; positive log2FC = higher in AD",
                          sum(meta$condition == "AD"), sum(meta$condition == "control")),
       x = "log2 fold change (AD − control, apeglm-shrunk)", y = "-log10(adjusted p)", colour = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(outdir, "figures", "fig1_volcano.png"), p, width = 8, height = 6.5, dpi = 200)
cat("\nSaved: results/rnaseq_DEGs.csv, figures/volcano.png, results/02_run_log.txt\n\n")

# ---- Step 9: session info ----------------------------------------------------
print(sessioninfo::session_info())
sink()
