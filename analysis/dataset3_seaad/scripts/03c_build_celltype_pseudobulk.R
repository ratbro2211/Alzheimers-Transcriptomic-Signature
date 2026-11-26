# =============================================================================
# 03c_build_celltype_pseudobulk.R — SEA-AD donor x cell-type pseudobulk (raw UMIs)
# -----------------------------------------------------------------------------
# Why: whole-tissue pseudobulk (script 03) mixes cell types, and donors differ in
# which cells were sampled, so its DEGs partly reflect composition (validation
# V6-V9, FAIL on the permutation test). Summing UMIs within each cell type per
# donor compares like with like: astrocytes with astrocytes, Sst neurons with
# Sst neurons.
# Cell types: 24 SEA-AD subclasses grouped into 10 classes (same grouping as
# fig5, with oligodendrocytes and OPCs kept separate; vascular cells kept but
# usually too few to analyse).
# Output: results/celltype/seaad_celltype_pseudobulk.rds
#   list(counts = genes x samples (donor__celltype), coldata = sample table)
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
proj  <- ad_root()
h5dir <- file.path(proj, "data", "SEA-AD-MTG", "donors_objects")
out   <- file.path(proj, "analysis", "dataset3_seaad", "results", "celltype")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(hdf5r))

GROUP <- c("L2/3 IT" = "Exc_IT", "L4 IT" = "Exc_IT", "L5 IT" = "Exc_IT", "L6 IT" = "Exc_IT", "L6 IT Car3" = "Exc_IT",
           "L5 ET" = "Exc_deep", "L5/6 NP" = "Exc_deep", "L6 CT" = "Exc_deep", "L6b" = "Exc_deep",
           "Pvalb" = "Inh_Pvalb", "Chandelier" = "Inh_Pvalb", "Sst" = "Inh_Sst", "Sst Chodl" = "Inh_Sst",
           "Vip" = "Inh_CGE", "Lamp5" = "Inh_CGE", "Lamp5 Lhx6" = "Inh_CGE", "Sncg" = "Inh_CGE", "Pax6" = "Inh_CGE",
           "Astrocyte" = "Astrocyte", "Immune" = "Microglia", "Oligodendrocyte" = "Oligodendrocyte", "OPC" = "OPC",
           "Endothelial" = "Vascular", "VLMC & Perivascular" = "Vascular")
CT <- unique(GROUP)

read_obs <- function(f, col) { x <- f[["obs"]][[col]]
  if (inherits(x, "H5Group")) { cats <- x[["categories"]][]; cats[x[["codes"]][] + 1] } else x[] }

files <- sort(list.files(h5dir, pattern = "\\.h5ad$", full.names = TRUE))
blocks <- list(); meta <- list(); genes <- NULL
for (fp in files) {
  f <- H5File$new(fp, mode = "r")
  g <- f[["var/index"]][]; if (is.null(genes)) genes <- g else stopifnot(identical(genes, g))
  sub <- read_obs(f, "Subclass"); stopifnot(all(sub %in% names(GROUP)))
  ct <- match(GROUP[sub], CT)                                  # cell-type index per nucleus
  donor <- read_obs(f, "Donor ID")[1]
  ip <- f[["layers/UMIs/indptr"]][]                             # CSR: row = nucleus
  dat <- f[["layers/UMIs/data"]]; ind <- f[["layers/UMIs/indices"]]
  nnz <- dat$dims; ng <- length(genes)
  acc <- matrix(0, ng, length(CT))
  row_of <- rep.int(seq_along(ct), diff(ip))                    # nucleus index of each non-zero
  for (st in seq(1, nnz, by = 2e7)) {
    en <- min(st + 2e7 - 1, nnz)
    v <- dat[st:en]; stopifnot(all(v == round(v)))              # raw integer UMIs only
    key <- (ct[row_of[st:en]] - 1L) * ng + ind[st:en] + 1L      # (cell type, gene) -> one index
    rs <- rowsum(v, key, reorder = FALSE); k <- as.integer(rownames(rs))
    acc[k] <- acc[k] + rs[, 1]
  }
  ncell <- tabulate(ct, length(CT))
  colnames(acc) <- paste(donor, CT, sep = "__")
  blocks[[donor]] <- acc[, ncell > 0, drop = FALSE]
  meta[[donor]] <- data.frame(sample = colnames(acc), donor_id = donor, celltype = CT, n_nuclei = ncell,
                              cognitive_status = read_obs(f, "Cognitive Status")[1],
                              ADNC = read_obs(f, "Overall AD neuropathological Change")[1],
                              age = read_obs(f, "Age at Death")[1], sex = read_obs(f, "Sex")[1])[ncell > 0, ]
  cat(sprintf("%-11s %6d nuclei  %s\n", donor, length(ct), paste(sprintf("%s=%d", CT, ncell), collapse = " ")))
  f$close_all(); rm(acc, row_of); gc(verbose = FALSE)
}
counts <- do.call(cbind, blocks); rownames(counts) <- genes; storage.mode(counts) <- "integer"
coldata <- do.call(rbind, meta); rownames(coldata) <- coldata$sample
coldata$condition <- factor(ifelse(coldata$cognitive_status == "Dementia", "AD", "control"), levels = c("control", "AD"))
stopifnot(identical(colnames(counts), rownames(coldata)))
saveRDS(list(counts = counts, coldata = coldata), file.path(out, "seaad_celltype_pseudobulk.rds"))
# sanity: summing cell types back must equal the whole-tissue pseudobulk from script 03
wt <- readRDS(file.path(proj, "analysis/dataset3_seaad/results/seaad_pseudobulk.rds"))$counts
back <- sapply(colnames(wt), function(d) rowSums(counts[, coldata$donor_id == d, drop = FALSE]))
cat(sprintf("\nCheck: cell-type sums reproduce whole-tissue pseudobulk exactly: %s\n", identical(unname(back), unname(wt[rownames(counts), ] * 1.0)) || max(abs(back - wt)) == 0))
cat(sprintf("Saved %d samples (donor x cell type), %d genes\n", ncol(counts), nrow(counts)))
