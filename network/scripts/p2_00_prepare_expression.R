# =============================================================================
# p2_00_prepare_expression.R — expression matrices for the network analysis
# -----------------------------------------------------------------------------
# Written in R rather than Python: microarray probe -> gene mapping needs
# illuminaHumanv4.db and RNA-seq/pseudobulk normalisation needs DESeq2, neither
# of which exists in Python. All later network steps are Python.
#
# Rebuilds per-sample, all-gene matrices exactly as analysis scripts 01-03 do;
# the analysis results/ folders contain DE tables,
# not expression matrices. No downloads: all inputs are local GEO files.
#   ARR  GSE132903 log2 intensity, highest-mean probe per gene (as script 01)
#   RNA  GSE104704 log2(DESeq2-normalised count + 1), raw counts recovered
#        exactly from the GEO table, AD vs Old (Young excluded)  (as script 02)
#   SEA  SEA-AD donor pseudobulk log2(DESeq2-normalised UMI + 1)  (script 03 rds)
# Writes data/expr_{ARR,RNA,SEA}.csv.gz (genes x samples), data/labels_*.csv
# (group = AD / ND), and results/panel_features_*.csv (samples x 31 panel genes).
# Run from network/.
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
suppressPackageStartupMessages({ library(DESeq2); library(illuminaHumanv4.db) })
root <- ad_root()
ana <- file.path(root, "analysis")
fld <- function(h, tag) gsub('"', '', strsplit(h[grep(paste0("^", tag), h)][1], "\t")[[1]][-1])
save_ds <- function(name, X, grp) {
  stopifnot(identical(colnames(X), names(grp)))
  write.csv(round(X, 4), gzfile(sprintf("data/expr_%s.csv.gz", name)))
  write.csv(data.frame(group = grp, row.names = names(grp)), sprintf("data/labels_%s.csv", name))
  cat(sprintf("%s: %d genes x %d samples (AD %d / ND %d)\n", name, nrow(X), ncol(X), sum(grp == "AD"), sum(grp == "ND")))
}

# ---- ARR ---------------------------------------------------------------------
mx <- file.path(root, "data", "GSE132903", "GSE132903_matrix.txt.gz")
con <- gzfile(mx, "rt"); h <- readLines(con, n = 200); close(con)
gsm <- fld(h, "!Sample_geo_accession"); tit <- fld(h, "!Sample_title")
E <- as.matrix(read.delim(mx, comment.char = "!", row.names = 1, check.names = FALSE)[, gsm])
sym <- unname(mapIds(illuminaHumanv4.db, keys = rownames(E), column = "SYMBOL", keytype = "PROBEID", multiVals = "first"))
ok <- !is.na(sym) & sym != ""; o <- order(sym[ok], -rowMeans(E[ok, ])); E <- E[ok, ][o, ]; s <- sym[ok][o]
E <- E[!duplicated(s), ]; rownames(E) <- s[!duplicated(s)]
stopifnot(nrow(E) == nrow(read.csv(file.path(ana, "dataset1_microarray/results/microarray_DEGs.csv"))))  # same gene set as script 01
save_ds("ARR", E, setNames(ifelse(grepl("_AD_", tit), "AD", "ND"), gsm))

# ---- RNA ---------------------------------------------------------------------
tab <- read.delim(file.path(root, "data/GSE104704/GSE104704_RNA-Seq_Table.txt.gz"), row.names = 1, check.names = FALSE)
nm <- as.matrix(tab[, setdiff(colnames(tab), c("Young_vs_Old", "Old_vs_AD", "Young_vs_AD"))])
raw <- sweep(nm, 2, apply(nm, 2, function(v) min(v[v > 0])), "/"); stopifnot(max(abs(raw - round(raw))) < 1e-6)
cnt <- round(raw); storage.mode(cnt) <- "integer"
con <- gzfile(file.path(root, "data/GSE104704/GSE104704_series_matrix.txt.gz"), "rt"); h2 <- readLines(con); close(con)
t2 <- fld(h2, "!Sample_title"); grp <- setNames(sub(".*\\.", "", t2), sub("\\.[^.]+$", "", t2))
keep <- names(grp)[grp %in% c("Old", "AD")]; cnt <- cnt[, keep]
cnt <- cnt[rowSums(cnt >= 10) >= 10, ]                         # same expression filter as script 02
dds <- estimateSizeFactors(DESeqDataSetFromMatrix(cnt, data.frame(row.names = keep, g = factor(grp[keep])), ~ 1))
save_ds("RNA", log2(counts(dds, normalized = TRUE) + 1), setNames(ifelse(grp[keep] == "AD", "AD", "ND"), keep))

# ---- SEA ---------------------------------------------------------------------
pb <- readRDS(file.path(ana, "dataset3_seaad/results/seaad_pseudobulk.rds"))
cn <- pb$counts[rowSums(pb$counts >= 10) >= 5, ]                # same filter as script 03
d3 <- estimateSizeFactors(DESeqDataSetFromMatrix(cn, pb$coldata, ~ 1))
save_ds("SEA", log2(counts(d3, normalized = TRUE) + 1),
        setNames(ifelse(pb$coldata$condition == "AD", "AD", "ND"), pb$coldata$donor_id))

# ---- 31-gene panel feature matrices --------------------------------------------
panel <- read.csv(file.path(ana, "comparison/results/panel_33gene_summary.csv"))$gene
for (k in c("ARR", "RNA", "SEA")) {
  X <- read.csv(gzfile(sprintf("data/expr_%s.csv.gz", k)), row.names = 1, check.names = FALSE)
  L <- read.csv(sprintf("data/labels_%s.csv", k), row.names = 1)
  g <- intersect(panel, rownames(X)); miss <- setdiff(panel, rownames(X))
  write.csv(data.frame(t(X[g, ]), condition = as.integer(L[colnames(X), "group"] == "AD"), check.names = FALSE),
            sprintf("results/panel_features_%s.csv", k))
  cat(sprintf("  panel genes in %s: %d/33%s\n", k, length(g), if (length(miss)) paste0(" (missing: ", paste(miss, collapse = ", "), ")") else ""))
}
