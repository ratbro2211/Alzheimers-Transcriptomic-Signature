# =============================================================================
# 03_pseudobulk_deseq2.R  —  Dataset 3: SEA-AD MTG snRNA-seq, donor pseudobulk
# -----------------------------------------------------------------------------
# Data   : SEA-AD middle temporal gyrus, one h5ad per donor (20 donors) in
#          data/SEA-AD-MTG/donors_objects/. (Reference_MTG_...2022 h5ad is the
#          neurotypical reference atlas, not this cohort, and is not used.)
# Input  : RAW UMI counts = adata.layers["UMIs"].
#          adata.X is log-normalised, ln(UP10K+1) (adata.uns["X_normalization"]),
#          and is NEVER used as DESeq2 input. There is no adata.raw in these files.
# Model  : sum UMIs over all nuclei per donor -> DESeq2, ~ condition,
#          reference = control, apeglm-shrunk log2FC (AD − control)
#
# Usage  : Rscript 03_pseudobulk_deseq2.R --build-only   # steps 1-3, then STOP
#          Rscript 03_pseudobulk_deseq2.R                # full run (steps 1-8)
# The h5ad files are read natively with hdf5r (see note below).
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
build_only <- "--build-only" %in% commandArgs(trailingOnly = TRUE)

proj   <- ad_root()
h5dir  <- file.path(proj, "data", "SEA-AD-MTG", "donors_objects")
outdir <- file.path(proj, "analysis", "dataset3_seaad")
logf   <- file.path(outdir, "results", if (build_only) "03_build_log.txt" else "03_run_log.txt")
sink(logf, split = TRUE)

suppressPackageStartupMessages(library(hdf5r))   # native HDF5 reader (h5ad = HDF5)

# NOTE ON THE READER: the plan was reticulate + Python anndata, but the miniconda
# Python on this machine is an x86_64 build while R is arm64, so reticulate cannot
# load it into the R process. hdf5r reads the h5ad (an HDF5 file) directly instead,
# which also avoids loading the per-cell matrix into memory.

# Read one obs column (anndata stores strings as categorical codes + categories)
read_obs <- function(f, col) {
  x <- f[["obs"]][[col]]
  if (inherits(x, "H5Group")) { cats <- x[["categories"]][]; cats[x[["codes"]][] + 1] } else x[]
}

describe <- function(path) {
  f <- H5File$new(path, mode = "r"); on.exit(f$close_all())
  shp <- h5attr(f[["layers/UMIs"]], "shape")
  list(n_obs = shp[1], n_vars = shp[2], layers = names(f[["layers"]]),
       has_raw = "raw" %in% names(f), X_normalization = f[["uns/X_normalization"]][],
       obs_keys = names(f[["obs"]]))
}

# Sum raw UMIs per gene over all nuclei of one donor. layers/UMIs is CSR
# (rows = nuclei, `indices` = gene columns); summing `data` grouped by `indices`
# gives the pseudobulk. Read in chunks of 20M non-zeros to bound memory.
donor_pseudobulk <- function(path, chunk = 2e7) {
  f <- H5File$new(path, mode = "r"); on.exit(f$close_all())
  genes <- f[["var/index"]][]
  dat <- f[["layers/UMIs/data"]]; ind <- f[["layers/UMIs/indices"]]
  nnz <- dat$dims; sums <- numeric(length(genes)); umax <- 0; is_int <- TRUE
  for (st in seq(1, nnz, by = chunk)) {
    en <- min(st + chunk - 1, nnz)
    v <- dat[st:en]; g <- ind[st:en] + 1L                 # 0-based -> 1-based gene index
    umax <- max(umax, max(v)); is_int <- is_int && all(v == round(v))
    rs <- rowsum(v, g, reorder = FALSE)
    k <- as.integer(rownames(rs)); sums[k] <- sums[k] + rs[, 1]
  }
  sub <- read_obs(f, "Subclass")
  one <- function(col) read_obs(f, col)[1]
  list(genes = genes, umi = sums, umi_max = umax, umi_integer = is_int, n_cells = length(sub),
       donor = one("Donor ID"), cognitive = one("Cognitive Status"),
       adnc = one("Overall AD neuropathological Change"), braak = one("Braak"),
       age = as.character(one("Age at Death")), sex = one("Sex"), comp = table(sub))
}

files <- sort(list.files(h5dir, pattern = "\\.h5ad$", full.names = TRUE))
cat(sprintf("Donor h5ad files found: %d\n", length(files)))
stopifnot(length(files) > 0)

# ---- Step 1: inspect the file structure (first donor) ------------------------
d <- describe(files[1])
cat("\nFile structure (first donor):\n")
cat("  n_obs x n_vars      :", d$n_obs, "x", d$n_vars, "\n")
cat("  layers              :", paste(unlist(d$layers), collapse = ", "), "\n")
cat("  adata.raw present   :", d$has_raw, "\n")
cat("  X_normalization     :", d$X_normalization, "  -> .X is NOT used\n")
cat("  obs keys (selected) :", paste(intersect(unlist(d$obs_keys),
    c("Donor ID", "Cognitive Status", "Overall AD neuropathological Change", "Braak",
      "Class", "Subclass", "Supertype")), collapse = ", "), "\n")
stopifnot("UMIs" %in% unlist(d$layers))

# ---- Step 2-3: per-donor raw-UMI pseudobulk (all nuclei, no cell-type subset) -
pb <- list(); meta <- list(); comp <- list(); genes <- NULL
for (f in files) {
  r <- donor_pseudobulk(f)
  # CRITICAL CHECKS: raw integer counts with a realistic dynamic range
  if (r$umi_max < 50) stop("ERROR: These look like normalised values, not raw counts.")
  if (!r$umi_integer) stop("ERROR: UMI layer contains non-integer values: ", basename(f))
  if (is.null(genes)) genes <- r$genes else stopifnot(identical(genes, r$genes))
  pb[[r$donor]] <- r$umi
  meta[[r$donor]] <- data.frame(donor_id = r$donor, cognitive_status = r$cognitive,
                                n_cells = r$n_cells, ADNC = r$adnc, braak = r$braak,
                                age = r$age, sex = r$sex, max_UMI_in_a_cell = r$umi_max)
  comp[[r$donor]] <- setNames(as.integer(r$comp), names(r$comp))
  cat(sprintf("  %-11s %6d nuclei  max UMI/cell-gene = %5.0f  integer = %s\n",
              r$donor, r$n_cells, r$umi_max, r$umi_integer))
}
counts <- do.call(cbind, pb); rownames(counts) <- genes
storage.mode(counts) <- "integer"
coldata <- do.call(rbind, meta); rownames(coldata) <- coldata$donor_id

# Condition = cognitive status (SEA-AD clinical diagnosis); ADNC kept as information
coldata$condition <- factor(ifelse(coldata$cognitive_status == "Dementia", "AD",
                            ifelse(coldata$cognitive_status == "No dementia", "control", NA)),
                            levels = c("control", "AD"))
stopifnot(!anyNA(coldata$condition))
coldata <- coldata[order(coldata$condition, coldata$donor_id), ]
counts  <- counts[, rownames(coldata)]

cat(sprintf("\nPseudobulk matrix: %d genes x %d donors; max = %d; total UMIs %.2e\n",
            nrow(counts), ncol(counts), max(counts), sum(as.numeric(counts))))
if (max(counts) < 50) stop("ERROR: These look like normalised values, not raw counts.")

# Cell-type composition (n nuclei per Subclass per donor)
ct <- sort(unique(unlist(lapply(comp, names))))
compm <- t(sapply(rownames(coldata), function(dn) { v <- comp[[dn]][ct]; v[is.na(v)] <- 0; v }))
colnames(compm) <- ct
write.csv(compm, file.path(outdir, "results", "seaad_celltype_counts_per_donor.csv"))
cat("\nCell-type composition (n nuclei per Subclass per donor):\n")
print(compm)

# Save the pseudobulk so the DESeq2 stage reproduces exactly what was reviewed
saveRDS(list(counts = counts, coldata = coldata), file.path(outdir, "results", "seaad_pseudobulk.rds"))
write.csv(coldata, file.path(outdir, "results", "seaad_pseudobulk_coldata.csv"), row.names = FALSE)

cat("\n================ PSEUDOBULK COLDATA (verify before DESeq2) ================\n")
print(coldata[, c("donor_id", "condition", "cognitive_status", "n_cells", "ADNC", "braak", "age", "sex")],
      row.names = FALSE)
print(table(coldata$condition))
cat("===========================================================================\n")

if (build_only) {
  cat("\n--build-only: stopping before DESeq2. Re-run without the flag to continue.\n")
  sink(); quit(save = "no", status = 0)
}

# ---- Step 4: DESeq2 ----------------------------------------------------------
suppressPackageStartupMessages({ library(DESeq2); library(apeglm); library(ggplot2); library(ggrepel) })
dds <- DESeqDataSetFromMatrix(counts, colData = coldata, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "control")
# "Expressed" = at least 10 UMIs in the donor pseudobulk; keep genes expressed in >= 5 donors
keep <- rowSums(counts(dds) >= 10) >= 5
dds  <- dds[keep, ]
cat(sprintf("\nGenes expressed (>=10 UMIs) in >=5 donors: %d\n", nrow(dds)))
dds <- DESeq(dds, quiet = TRUE)
stopifnot("condition_AD_vs_control" %in% resultsNames(dds))
res <- results(dds, name = "condition_AD_vs_control", alpha = 0.05)
shr <- lfcShrink(dds, coef = "condition_AD_vs_control", type = "apeglm", quiet = TRUE)

# ---- Step 5: results table ---------------------------------------------------
PADJ <- 0.05; LFC <- 0.5
out <- data.frame(gene = rownames(res), log2FC = shr$log2FoldChange, pvalue = res$pvalue, padj = res$padj)
out$direction <- ifelse(!is.na(out$padj) & out$padj < PADJ & out$log2FC >  LFC, "Up",
                 ifelse(!is.na(out$padj) & out$padj < PADJ & out$log2FC < -LFC, "Down", "NS"))
out <- out[order(out$padj, out$pvalue, na.last = TRUE), ]
write.csv(out, file.path(outdir, "results", "seaad_DEGs.csv"), row.names = FALSE)
cat(sprintf("Genes tested: %d | padj NA: %d | significant: %d (Up %d / Down %d)\n",
            nrow(out), sum(is.na(out$padj)), sum(out$direction != "NS"),
            sum(out$direction == "Up"), sum(out$direction == "Down")))
cat("\nTop 20 genes by padj:\n"); print(head(out, 20), row.names = FALSE, digits = 3)

# ---- Step 7: validation ------------------------------------------------------
expect <- c(VGF = "Down", NPTX2 = "Down", SST = "Down", GFAP = "Up",
            BDNF = "Down", CD44 = "Up", SERPINA3 = "Up")
v <- out[match(names(expect), out$gene), c("gene", "log2FC", "padj", "direction")]
v$expected <- expect
v$sign_ok  <- ifelse(sign(v$log2FC) == ifelse(expect == "Up", 1, -1), "yes", "NO")
cat("\nValidation — canonical AD markers:\n"); print(v, row.names = FALSE, digits = 3)

# ---- Step 6: volcano ---------------------------------------------------------
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
  labs(title = "SEA-AD MTG snRNA-seq pseudobulk — dementia vs no dementia",
       subtitle = sprintf("raw UMIs summed per donor; DESeq2 + apeglm; %d AD / %d control",
                          sum(coldata$condition == "AD"), sum(coldata$condition == "control")),
       x = "log2 fold change (AD − control, apeglm-shrunk)", y = "-log10(adjusted p)", colour = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(outdir, "figures", "fig1_volcano.png"), p, width = 8, height = 6.5, dpi = 200)
cat("\nSaved: results/seaad_DEGs.csv, figures/volcano.png\n")

# ---- Step 7b: second model with nuclei count as a covariate -------------------
# Donors differ ~10-fold in nuclei recovered (2,918-27,431), and AD donors have
# fewer on average. n_cells enters as a centred, scaled log2 value so the GLM
# coefficient is well-conditioned; condition is still last in the design, so the
# reported effect is AD vs control adjusted for nuclei count.
coldata$log2_ncells_z <- as.numeric(scale(log2(coldata$n_cells)))
cat(sprintf("\nn_cells by group (median): control %d, AD %d; Wilcoxon p = %.3f\n",
            as.integer(median(coldata$n_cells[coldata$condition == "control"])),
            as.integer(median(coldata$n_cells[coldata$condition == "AD"])),
            wilcox.test(n_cells ~ condition, data = coldata)$p.value))
dds2 <- DESeqDataSetFromMatrix(counts[rownames(dds), ], colData = coldata, design = ~ log2_ncells_z + condition)
dds2$condition <- relevel(dds2$condition, ref = "control")
dds2 <- DESeq(dds2, quiet = TRUE)
stopifnot("condition_AD_vs_control" %in% resultsNames(dds2))
res2 <- results(dds2, name = "condition_AD_vs_control", alpha = 0.05)
shr2 <- lfcShrink(dds2, coef = "condition_AD_vs_control", type = "apeglm", quiet = TRUE)
out2 <- data.frame(gene = rownames(res2), log2FC = shr2$log2FoldChange, pvalue = res2$pvalue, padj = res2$padj)
out2$direction <- ifelse(!is.na(out2$padj) & out2$padj < PADJ & out2$log2FC >  LFC, "Up",
                  ifelse(!is.na(out2$padj) & out2$padj < PADJ & out2$log2FC < -LFC, "Down", "NS"))
out2 <- out2[order(out2$padj, out2$pvalue, na.last = TRUE), ]
write.csv(out2, file.path(outdir, "results", "seaad_DEGs_ncells_adjusted.csv"), row.names = FALSE)
cat(sprintf("n_cells-adjusted model: significant %d (Up %d / Down %d)\n",
            sum(out2$direction != "NS"), sum(out2$direction == "Up"), sum(out2$direction == "Down")))

# Validation, both models side by side
v2 <- out2[match(names(expect), out2$gene), c("log2FC", "padj", "direction")]
vt <- data.frame(gene = names(expect), expected = expect,
                 main_log2FC = v$log2FC, main_padj = v$padj, main_dir = v$direction,
                 ncells_log2FC = v2$log2FC, ncells_padj = v2$padj, ncells_dir = v2$direction)
cat("\nValidation — main model vs n_cells-adjusted model:\n"); print(vt, row.names = FALSE, digits = 3)

# 33-gene direction-consistent candidate panel
panel33 <- c("VGF","SST","SCG2","CARTPT","CHGB","VSNL1","TAC1","SNAP25","NPTX2","NEFL","GABRA1","NRGN",
             "GABRG2","SYT4","CCK","NEFM","SHANK2","GRIN2A","NECAB1","LRRTM1","CBLN2","CAMK2A","SYN2",
             "HOMER1","NPY","SCN1A","PENK","NLGN1","CADM1","VIM","TREM2","GJA1","GFAP")
pt <- data.frame(gene = panel33,
                 main_log2FC = out$log2FC[match(panel33, out$gene)], main_padj = out$padj[match(panel33, out$gene)],
                 main_dir = out$direction[match(panel33, out$gene)],
                 ncells_log2FC = out2$log2FC[match(panel33, out2$gene)], ncells_padj = out2$padj[match(panel33, out2$gene)],
                 ncells_dir = out2$direction[match(panel33, out2$gene)])
write.csv(pt, file.path(outdir, "results", "seaad_panel33.csv"), row.names = FALSE)
cat("\n33-gene panel in SEA-AD:\n"); print(pt[order(pt$main_padj), ], row.names = FALSE, digits = 3)
for (m in c("main", "ncells")) {
  d <- pt[[paste0(m, "_dir")]]; q <- pt[[paste0(m, "_padj")]]
  cat(sprintf("%-6s model: tested %d/33 | padj<0.05: %d | padj<0.05 & |log2FC|>0.5: %d (%s)\n", m,
              sum(!is.na(q)), sum(q < 0.05, na.rm = TRUE), sum(d != "NS", na.rm = TRUE),
              paste(pt$gene[!is.na(d) & d != "NS"], collapse = ", ")))
}
cat("\nSaved: results/seaad_DEGs_ncells_adjusted.csv, results/seaad_panel33.csv, results/03_run_log.txt\n\n")

# ---- Step 8: session info ----------------------------------------------------
print(sessioninfo::session_info())
sink()
