# =============================================================================
# 05d_validate_comparison.R — validation of the cross-platform comparison +
# assembly of the final report (validation/VALIDATION_REPORT.md).
# Tests whether the concordance could arise by chance, by gene-matching
# artefacts, by panel selection (circularity) or by single influential samples.
# Run after 05a-05c.
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
vdir <- file.path(proj, "validation"); rdir <- file.path(vdir, "results"); fdir <- file.path(vdir, "figures")
suppressPackageStartupMessages({ library(ggplot2) })
sink(file.path(rdir, "05d_log.txt"), split = TRUE)
checks <- data.frame()
add <- function(id, check, value, verdict, note = "") {
  checks <<- rbind(checks, data.frame(dataset = "Comparison", id = id, check = check, value = value, verdict = verdict, note = note))
  cat(sprintf("[%s] %-4s %-50s %s  %s\n", verdict, id, check, value, note))
}
PADJ <- 0.05; LFC <- 0.5; DS <- c(micro = "Microarray", rnaseq = "RNA-seq", seaad = "SEA-AD")
PANEL33 <- c("VGF","SST","SCG2","CARTPT","CHGB","VSNL1","TAC1","SNAP25","NPTX2","NEFL","GABRA1","NRGN","GABRG2","SYT4","CCK","NEFM",
             "SHANK2","GRIN2A","NECAB1","LRRTM1","CBLN2","CAMK2A","SYN2","HOMER1","NPY","SCN1A","PENK","NLGN1","CADM1","VIM","TREM2","GJA1","GFAP")
rd <- function(f) { d <- read.csv(f); d <- d[!duplicated(d$gene) & !is.na(d$padj), ]; d$sig <- d$direction != "NS"; d }
L <- list(micro = rd(file.path(proj, "dataset1_microarray/results/microarray_DEGs.csv")),
          rnaseq = rd(file.path(proj, "dataset2_rnaseq/results/rnaseq_DEGs.csv")),
          seaad = rd(file.path(proj, "dataset3_seaad/results/seaad_DEGs.csv")))
pw <- function(A, B) {                                # pair stats from two DE tables
  m <- merge(A[, c("gene", "log2FC", "sig")], B[, c("gene", "log2FC", "sig")], by = "gene")
  both <- m$sig.x & m$sig.y
  list(m = m, n = nrow(m), r = cor(m$log2FC.x, m$log2FC.y), rho = cor(m$log2FC.x, m$log2FC.y, method = "spearman"),
       nb = sum(both), same = sum(both & sign(m$log2FC.x) == sign(m$log2FC.y)), sa = sum(m$sig.x), sb = sum(m$sig.y))
}
pairs <- list(c("micro", "rnaseq"), c("micro", "seaad"), c("rnaseq", "seaad"))

# ---- C1 gene-symbol matching ----------------------------------------------------
all3 <- Reduce(intersect, lapply(L, `[[`, "gene"))
fr <- sapply(L, function(d) mean(d$gene %in% all3))
add("C1", "Gene-symbol matching across platforms", sprintf("%d genes in all three; matched %s", length(all3),
    paste(sprintf("%s %.0f%%", DS, 100 * fr), collapse = ", ")), ifelse(min(fr) > 0.4, "PASS", "WARN"),
    "unmatched genes are mostly platform-specific (lncRNAs, unannotated probes)")

null_rows <- list()
for (p in pairs) {
  s <- pw(L[[p[1]]], L[[p[2]]]); lab <- paste(DS[p[1]], "vs", DS[p[2]])
  # C2 overlap of significant sets vs chance
  exp_ov <- s$sa * s$sb / s$n; ov <- sum(s$m$sig.x & s$m$sig.y)
  ph <- phyper(ov - 1, s$sa, s$n - s$sa, s$sb, lower.tail = FALSE)
  add("C2", paste("Overlap above chance:", lab), sprintf("%d shared vs %.1f expected (%.1fx), hypergeometric p = %.1e", ov, exp_ov, ov / exp_ov, ph),
      ifelse(ph < 0.001, "PASS", "WARN"))
  # C3 direction agreement vs 50/50
  pb <- binom.test(s$same, s$nb, 0.5)$p.value
  add("C3", paste("Direction agreement vs coin-flip:", lab), sprintf("%d/%d (%.0f%%), binomial p = %.1e", s$same, s$nb, 100 * s$same / s$nb, pb),
      ifelse(pb < 0.001 & s$same / s$nb > 0.8, "PASS", "WARN"))
  # C4 correlation CI (bootstrap over genes) + rank correlation
  bs <- replicate(1000, { i <- sample(s$n, replace = TRUE); cor(s$m$log2FC.x[i], s$m$log2FC.y[i]) })
  add("C4", paste("Correlation robust (bootstrap CI, Spearman):", lab), sprintf("r = %.2f [95%% CI %.2f-%.2f]; Spearman rho = %.2f",
      s$r, quantile(bs, .025), quantile(bs, .975), s$rho), ifelse(quantile(bs, .025) > 0, "PASS", "FAIL"),
      "Spearman is insensitive to apeglm shrinkage of magnitudes")
  # C5 gene-label permutation null
  nl <- t(replicate(1000, { y <- sample(s$m$log2FC.y); sy <- sample(s$m$sig.y); b <- s$m$sig.x & sy
    c(r = cor(s$m$log2FC.x, y), agree = if (sum(b)) mean(sign(s$m$log2FC.x[b]) == sign(y[b])) else NA) }))
  null_rows[[lab]] <- data.frame(pair = lab, null_r = nl[, "r"], null_agree = nl[, "agree"], obs_r = s$r, obs_agree = s$same / s$nb)
  add("C5", paste("Gene-shuffle null:", lab), sprintf("null r 95%% range %.3f-%.3f (obs %.2f); null agreement %.0f%% (obs %.0f%%)",
      quantile(nl[, "r"], .025), quantile(nl[, "r"], .975), s$r, 100 * mean(nl[, "agree"], na.rm = TRUE), 100 * s$same / s$nb),
      ifelse(s$r > quantile(nl[, "r"], .999), "PASS", "FAIL"))
}
write.csv(do.call(rbind, null_rows), file.path(rdir, "comparison_gene_shuffle_null.csv"), row.names = FALSE)

# ---- C6 panel circularity ------------------------------------------------------
# The 33-gene panel was picked from earlier results, so its consistency is not an
# independent test. Compare with (a) random genes and (b) random genes drawn the
# way a panel would be selected: significant in the microarray.
lf <- sapply(L, function(d) d$log2FC[match(all3, d$gene)]); rownames(lf) <- all3
cons <- apply(sign(lf), 1, function(s) length(unique(s)) == 1)
pan <- intersect(PANEL33, all3); obs <- sum(cons[pan])
rand_any <- replicate(5000, sum(cons[sample(all3, length(pan))]))
msig <- intersect(L$micro$gene[L$micro$sig], all3)
rand_sig <- replicate(5000, sum(cons[sample(msig, length(pan))]))
add("C6a", "33-gene panel vs random genes", sprintf("panel %d/%d consistent; random genes mean %.1f, P(>= obs) = %.4f", obs, length(pan),
    mean(rand_any), mean(rand_any >= obs)), "INFO", "random genes are consistent ~25% of the time by chance (2 of 8 sign patterns)")
add("C6b", "Panel vs random microarray-significant genes", sprintf("random sig genes mean %.1f/%d consistent, P(>= %d) = %.3f", mean(rand_sig), length(pan), obs,
    mean(rand_sig >= obs)), ifelse(mean(rand_sig >= obs) < 0.05, "PASS", "WARN"),
    "random microarray-significant genes are already ~79% consistent: most of the panel's consistency reflects replication of significant genes in general")
write.csv(data.frame(random_any = rand_any, random_microarray_sig = rand_sig, observed = obs), file.path(rdir, "panel_circularity_null.csv"), row.names = FALSE)

# ---- C7 robustness to the most influential sample in each small dataset ---------
fi <- function(k) { f <- list.files(rdir, pattern = sprintf("^%s_DEGs_drop_.*csv$", k), full.names = TRUE); if (length(f)) rd(f[1]) else NULL }
Ld <- L; dr <- fi("rnaseq"); ds <- fi("seaad")
if (!is.null(dr)) Ld$rnaseq <- dr
if (!is.null(ds)) Ld$seaad <- ds
for (p in pairs) {
  a <- pw(L[[p[1]]], L[[p[2]]]); b <- pw(Ld[[p[1]]], Ld[[p[2]]])
  add("C7", paste("Drop most-influential sample:", DS[p[1]], "vs", DS[p[2]]), sprintf("r %.2f -> %.2f; agreement %d/%d -> %d/%d",
      a$r, b$r, a$same, a$nb, b$same, b$nb), ifelse(abs(a$r - b$r) < 0.05 & (b$same / max(b$nb, 1)) > 0.85, "PASS", "WARN"),
      sprintf("dropped: %s", paste(c(readLines(file.path(rdir, "rnaseq_most_influential.txt")),
                                     if (file.exists(file.path(rdir, "seaad_most_influential.txt"))) readLines(file.path(rdir, "seaad_most_influential.txt"))), collapse = ", ")))
}

# ---- C8 effect-size attenuation (informational) --------------------------------
for (p in pairs) {
  s <- pw(L[[p[1]]], L[[p[2]]]); b <- s$m$sig.x & s$m$sig.y
  add("C8", paste("Effect-size ratio (|y|/|x|, sig in both):", DS[p[1]], "vs", DS[p[2]]),
      sprintf("median %.2f", median(abs(s$m$log2FC.y[b]) / abs(s$m$log2FC.x[b]))), "INFO", "<1 = second dataset shows smaller fold changes (shrinkage, power, composition)")
}
write.csv(checks, file.path(rdir, "checks_comparison.csv"), row.names = FALSE)

# =================== FIGURES ====================================================
th <- theme_bw(base_size = 11) + theme(plot.background = element_rect(fill = "white", colour = NA), plot.title = element_text(face = "bold"))
pv <- do.call(rbind, lapply(c("microarray", "rnaseq", "seaad"), function(k) { f <- file.path(rdir, sprintf("pvals_%s.csv", k)); if (file.exists(f)) read.csv(f) }))
g <- ggplot(pv, aes(p)) + geom_histogram(breaks = seq(0, 1, 0.025), fill = "grey40", colour = "white") + facet_wrap(~dataset, scales = "free_y") +
  labs(title = "p-value distributions", subtitle = "Expected: flat right side (null genes) + spike near 0 (true effects)", x = "p-value", y = "genes") + th
ggsave(file.path(fdir, "pvalue_histograms.png"), g, width = 10, height = 3.6, dpi = 300, bg = "white")
pn <- do.call(rbind, lapply(c("microarray", "rnaseq", "seaad"), function(k) { f <- file.path(rdir, sprintf("permnull_%s.csv", k)); if (file.exists(f)) read.csv(f) }))
g <- ggplot(pn, aes(null_DEGs)) + geom_histogram(bins = 25, fill = "grey55", colour = "white") +
  geom_vline(aes(xintercept = observed), colour = "#E63946", linewidth = 1) + facet_wrap(~dataset, scales = "free") +
  labs(title = "Label-permutation null", subtitle = "grey = DEGs found with shuffled AD/control labels; red line = observed DEGs", x = "number of DEGs", y = "permutations") + th
ggsave(file.path(fdir, "permutation_null.png"), g, width = 10, height = 3.6, dpi = 300, bg = "white")
lo <- do.call(rbind, lapply(c("rnaseq", "seaad"), function(k) { f <- file.path(rdir, sprintf("loo_%s.csv", k)); if (file.exists(f)) cbind(dataset = DS[k], read.csv(f)[, c("left_out", "condition", "n_DEGs", "r_lfc_vs_full", "pct_change_DEGs")]) }))
g <- ggplot(lo, aes(r_lfc_vs_full, pct_change_DEGs, colour = condition, label = left_out)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(size = 2.4) + ggrepel::geom_text_repel(size = 2.5, seed = 42, max.overlaps = 15) + facet_wrap(~dataset, scales = "free") +
  scale_colour_manual(values = c(control = "#457B9D", AD = "#E63946")) +
  labs(title = "Leave-one-sample-out sensitivity", x = "r of log2FC with full analysis", y = "% change in DEG count") + th
ggsave(file.path(fdir, "leave_one_out.png"), g, width = 11, height = 4.5, dpi = 300, bg = "white")
cn <- rbind(data.frame(null = "random genes", n = rand_any), data.frame(null = "random microarray-significant genes", n = rand_sig))
g <- ggplot(cn, aes(n)) + geom_histogram(binwidth = 1, fill = "grey55", colour = "white") + geom_vline(xintercept = obs, colour = "#E63946", linewidth = 1) +
  facet_wrap(~null, scales = "free_y") + labs(title = "Is the 33-gene panel special? (selection-bias check)",
  subtitle = sprintf("genes (of %d) with the same sign in all three datasets; red = panel (%d)", length(pan), obs), x = "direction-consistent genes", y = "random sets") + th
ggsave(file.path(fdir, "panel_circularity.png"), g, width = 10, height = 3.8, dpi = 300, bg = "white")

# =================== REPORT =====================================================
C <- do.call(rbind, lapply(c("microarray", "rnaseq", "seaad", "seaad_celltype", "comparison"), function(k) {
  f <- file.path(rdir, sprintf("checks_%s.csv", k)); if (file.exists(f)) read.csv(f) }))
write.csv(C, file.path(rdir, "all_checks.csv"), row.names = FALSE)
tally <- table(factor(C$dataset, levels = c("Microarray", "RNA-seq", "SEA-AD", "SEA-AD cell-type", "Comparison")), factor(C$verdict, levels = c("PASS", "WARN", "FAIL", "INFO")))
esc <- function(x) gsub("\\|", "/", x)
tbl <- function(d) c("| ID | Check | Result | Verdict | Note |", "|---|---|---|---|---|",
                     sprintf("| %s | %s | %s | **%s** | %s |", d$id, esc(d$check), esc(d$value), d$verdict, esc(d$note)))
rep <- c("# Validation report — analysis", "",
         sprintf("Generated %s by `validation/scripts/05a-05d`. Verdicts: PASS = no problem; WARN = a caveat to report; FAIL = result not trustworthy; INFO = descriptive.", Sys.Date()), "",
         "## Summary", "", "| Dataset | PASS | WARN | FAIL | INFO |", "|---|---|---|---|---|",
         sprintf("| %s | %d | %d | %d | %d |", rownames(tally), tally[, 1], tally[, 2], tally[, 3], tally[, 4]), "")
for (d in rownames(tally)) rep <- c(rep, sprintf("## %s", d), "", tbl(C[C$dataset == d, ]), "")
rep <- c(rep, "## Figures", "", "- `figures/pvalue_histograms.png`, `figures/permutation_null.png`, `figures/leave_one_out.png`,",
         "  `figures/panel_circularity.png`, `figures/pca_microarray.png`, `figures/pca_rnaseq.png`, `figures/pca_seaad.png`")
ip <- file.path(vdir, "INTERPRETATION.md")                 # reviewer-written interpretation, kept in its own file
if (file.exists(ip)) rep <- c(rep[1:(grep("^## Summary", rep) - 1)], readLines(ip), "", rep[grep("^## Summary", rep):length(rep)])
writeLines(rep, file.path(vdir, "VALIDATION_REPORT.md"))
print(tally)
sink()
