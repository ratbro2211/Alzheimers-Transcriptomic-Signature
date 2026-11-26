# =============================================================================
# 06_ml_classification.R — AD classifier on the cross-platform gene panel
# -----------------------------------------------------------------------------
# Train on GSE132903 microarray (97 AD / 98 ND), evaluate by LOO-CV, then apply
# the frozen model to two independent datasets (GSE104704 RNA-seq, 12 AD / 10
# control; SEA-AD pseudobulk, 10 / 10) with no retraining or tuning.
#
# Features: the panel genes measured in all three datasets (31 of 33; CAMK2A and
# PENK have no annotated microarray probe).
# Expression: log2 scale in every dataset
#   microarray  log2 intensity (GEO series matrix; highest-mean probe per gene)
#   RNA-seq     log2(DESeq2-normalised count + 1) from raw counts recovered
#               exactly from the GEO table (see script 02)
#   SEA-AD      log2(DESeq2-normalised pseudobulk UMI + 1) (script 03 output)
# Cross-platform scaling: each gene is z-scored WITHIN each dataset (label-free).
#   Absolute values are not comparable across platforms (intensity vs counts);
#   within-dataset standardisation is the standard way to transfer a model.
#   In LOO-CV the scaling is fitted on the training fold only (no leakage).
# Models: L2 logistic regression (glmnet, alpha = 0, lambda by inner CV) and
#   random forest (500 trees). Decision threshold chosen on microarray LOO-CV
#   predictions to give sensitivity >= 0.85, then frozen for the test sets.
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
suppressPackageStartupMessages({ library(glmnet); library(randomForest); library(pROC); library(ggplot2)
  library(patchwork); library(parallel); library(DESeq2); library(illuminaHumanv4.db) })
root <- ad_root()
proj <- file.path(root, "analysis"); mdir <- file.path(proj, "ml_classification")
rdir <- file.path(mdir, "results"); fdir <- file.path(mdir, "figures")
sink(file.path(rdir, "06_run_log.txt"), split = TRUE)
CORES <- 4; SENS_TARGET <- 0.85; NTREE <- 500; LC_REPEATS <- 10
MCOL <- c("Logistic (L2)" = "#2a78d6", "Random forest" = "#eb6834")

# ============================ STEP 1: feature matrices ==========================
panel <- read.csv(file.path(proj, "comparison/results/panel_33gene_summary.csv"))
genes <- panel$gene[panel$n_datasets_measured == 3]
cat(sprintf("Panel: %d genes listed, %d measured in all three datasets; excluded: %s\n",
            nrow(panel), length(genes), paste(setdiff(panel$gene, genes), collapse = ", ")))

# -- microarray (local copy of the GEO series matrix) --
mx <- file.path(root, "data", "GSE132903", "GSE132903_matrix.txt.gz")
con <- gzfile(mx, "rt"); hdr <- readLines(con, n = 200); close(con)
fld <- function(h, tag) gsub('"', '', strsplit(h[grep(paste0("^", tag), h)][1], "\t")[[1]][-1])
gsm <- fld(hdr, "!Sample_geo_accession"); tit <- fld(hdr, "!Sample_title")
E <- as.matrix(read.delim(mx, comment.char = "!", row.names = 1, check.names = FALSE)[, gsm])
sym <- unname(mapIds(illuminaHumanv4.db, keys = rownames(E), column = "SYMBOL", keytype = "PROBEID", multiVals = "first"))
ok <- !is.na(sym) & sym != ""; o <- order(sym[ok], -rowMeans(E[ok, ])); E <- E[ok, ][o, ]; s <- sym[ok][o]
E <- E[!duplicated(s), ]; rownames(E) <- s[!duplicated(s)]
Xm <- t(E[intersect(genes, rownames(E)), ]); ym <- as.integer(grepl("_AD_", tit))

# -- RNA-seq (local copy of the GEO supplementary table; raw counts recovered exactly) --
tab <- read.delim(file.path(root, "data/GSE104704/GSE104704_RNA-Seq_Table.txt.gz"), row.names = 1, check.names = FALSE)
nm <- as.matrix(tab[, setdiff(colnames(tab), c("Young_vs_Old", "Old_vs_AD", "Young_vs_AD"))])
raw <- sweep(nm, 2, apply(nm, 2, function(v) min(v[v > 0])), "/"); stopifnot(max(abs(raw - round(raw))) < 1e-6)
cnt <- round(raw); storage.mode(cnt) <- "integer"
con <- gzfile(file.path(root, "data/GSE104704/GSE104704_series_matrix.txt.gz"), "rt"); h2 <- readLines(con); close(con)
t2 <- fld(h2, "!Sample_title"); grp <- setNames(sub(".*\\.", "", t2), sub("\\.[^.]+$", "", t2))
keep <- names(grp)[grp %in% c("Old", "AD")]; cnt <- cnt[, keep]
sf <- estimateSizeFactors(DESeqDataSetFromMatrix(cnt, data.frame(row.names = keep, g = factor(grp[keep])), ~ 1))
Lr <- log2(counts(sf, normalized = TRUE) + 1)
Xr <- t(Lr[intersect(genes, rownames(Lr)), ]); yr <- as.integer(grp[keep] == "AD")

# -- SEA-AD pseudobulk (raw UMIs from script 03) --
pb <- readRDS(file.path(proj, "dataset3_seaad/results/seaad_pseudobulk.rds"))
sf3 <- estimateSizeFactors(DESeqDataSetFromMatrix(pb$counts, pb$coldata, ~ 1))
Ls <- log2(counts(sf3, normalized = TRUE) + 1)
Xs <- t(Ls[intersect(genes, rownames(Ls)), ]); ys <- as.integer(pb$coldata$condition == "AD")

feat <- Reduce(intersect, list(colnames(Xm), colnames(Xr), colnames(Xs)))
for (nmx in c("Xm", "Xr", "Xs")) { miss <- setdiff(genes, colnames(get(nmx))); if (length(miss)) cat(nmx, "missing:", miss, "\n") }
Xm <- Xm[, feat]; Xr <- Xr[, feat]; Xs <- Xs[, feat]
cat(sprintf("Features used: %d genes\nMicroarray %d (AD %d / ND %d) | RNA-seq %d (AD %d / ctrl %d) | SEA-AD %d (AD %d / ctrl %d)\n",
            length(feat), nrow(Xm), sum(ym), sum(!ym), nrow(Xr), sum(yr), sum(!yr), nrow(Xs), sum(ys), sum(!ys)))
wr <- function(X, y, f) write.csv(data.frame(sample = rownames(X), round(X, 4), condition = y, check.names = FALSE),
                                  file.path(rdir, "feature_matrices", f), row.names = FALSE)
wr(Xm, ym, "microarray_features.csv"); wr(Xr, yr, "rnaseq_features.csv"); wr(Xs, ys, "seaad_features.csv")

# ============================ model helpers =====================================
zfit <- function(X) list(m = colMeans(X), s = apply(X, 2, sd))
zapp <- function(X, z) sweep(sweep(X, 2, z$m, "-"), 2, z$s, "/")
fit_lr <- function(X, y) { nf <- if (nrow(X) >= 60) 10 else 5
  cv.glmnet(X, y, family = "binomial", alpha = 0, nfolds = nf, foldid = sample(rep(seq_len(nf), length.out = nrow(X)))) }
pred_lr <- function(m, X) as.numeric(predict(m, X, s = "lambda.min", type = "response"))
fit_rf <- function(X, y) randomForest(X, factor(y, levels = 0:1), ntree = NTREE)
pred_rf <- function(m, X) as.numeric(predict(m, X, type = "prob")[, "1"])
MODELS <- list("Logistic (L2)" = list(fit = fit_lr, pred = pred_lr), "Random forest" = list(fit = fit_rf, pred = pred_rf))

loo <- function(X, y, mod, seed = 1) unlist(mclapply(seq_len(nrow(X)), function(i) {
  set.seed(seed * 1000 + i); z <- zfit(X[-i, , drop = FALSE])
  m <- mod$fit(zapp(X[-i, , drop = FALSE], z), y[-i]); mod$pred(m, zapp(X[i, , drop = FALSE], z)) }, mc.cores = CORES))

metrics <- function(y, p, thr) { pr <- as.integer(p >= thr)
  tp <- sum(pr == 1 & y == 1); tn <- sum(pr == 0 & y == 0); fp <- sum(pr == 1 & y == 0); fn <- sum(pr == 0 & y == 1)
  r <- roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE); ci <- as.numeric(ci.auc(r, method = "delong"))
  data.frame(AUC = round(as.numeric(auc(r)), 3), AUC_CI_low = round(ci[1], 3), AUC_CI_high = round(ci[3], 3),
             AUC_wilcox_p = signif(wilcox.test(p[y == 1], p[y == 0])$p.value, 2), threshold = round(thr, 3),
             sensitivity = round(tp / (tp + fn), 3), specificity = round(tn / (tn + fp), 3),
             balanced_accuracy = round((tp / (tp + fn) + tn / (tn + fp)) / 2, 3), TP = tp, FN = fn, TN = tn, FP = fp) }
sens_threshold <- function(y, p) {      # highest threshold that still gives sensitivity >= target
  cand <- sort(unique(p), decreasing = TRUE); for (t in cand) if (mean(p[y == 1] >= t) >= SENS_TARGET) return(t); min(p) }

# ============================ STEP 2: LOO-CV on microarray ======================
cat("\n== Step 2: LOO-CV on microarray ==\n")
loo_pred <- lapply(MODELS, function(mod) loo(Xm, ym, mod))
thr <- sapply(names(MODELS), function(k) sens_threshold(ym, loo_pred[[k]]))
L <- do.call(rbind, lapply(names(MODELS), function(k) rbind(
  cbind(model = k, rule = sprintf("sensitivity >= %.2f", SENS_TARGET), metrics(ym, loo_pred[[k]], thr[k])),
  cbind(model = k, rule = "threshold 0.5", metrics(ym, loo_pred[[k]], 0.5)))))
write.csv(L, file.path(rdir, "loocv_results.csv"), row.names = FALSE); print(L, row.names = FALSE)
write.csv(data.frame(sample = rownames(Xm), condition = ym, lr = loo_pred[[1]], rf = loo_pred[[2]]), file.path(rdir, "loocv_predictions.csv"), row.names = FALSE)

# ============================ STEP 3: cross-dataset ==============================
cat("\n== Step 3: train on all microarray, predict RNA-seq and SEA-AD (frozen model + threshold) ==\n")
zm <- zfit(Xm); set.seed(42)
final <- lapply(MODELS, function(mod) mod$fit(zapp(Xm, zm), ym))
TEST <- list("RNA-seq (GSE104704)" = list(X = Xr, y = yr), "SEA-AD pseudobulk" = list(X = Xs, y = ys))
test_pred <- list()
Cx <- do.call(rbind, lapply(names(TEST), function(tn) { T <- TEST[[tn]]; Z <- zapp(T$X, zfit(T$X))   # scaled with its own (label-free) stats
  do.call(rbind, lapply(names(MODELS), function(k) { p <- MODELS[[k]]$pred(final[[k]], Z); test_pred[[paste(tn, k)]] <<- data.frame(sample = rownames(Z), condition = T$y, prob = p)
    rbind(cbind(dataset = tn, model = k, rule = "frozen microarray threshold (sens >= 0.85)", metrics(T$y, p, thr[k])),
          cbind(dataset = tn, model = k, rule = "threshold 0.5", metrics(T$y, p, 0.5))) })) }))
write.csv(Cx, file.path(rdir, "crossdataset_results.csv"), row.names = FALSE); print(Cx, row.names = FALSE)
write.csv(do.call(rbind, lapply(names(test_pred), function(k) cbind(set = k, test_pred[[k]]))), file.path(rdir, "crossdataset_predictions.csv"), row.names = FALSE)

# ============================ STEP 4: learning curve =============================
cat("\n== Step 4: learning curve (LOO-CV AUC on stratified subsamples) ==\n")
LC <- do.call(rbind, lapply(c(0.2, 0.4, 0.6, 0.8, 1.0), function(fr) {
  reps <- if (fr == 1) 1 else LC_REPEATS
  do.call(rbind, lapply(seq_len(reps), function(r) { set.seed(100 * fr + r)
    idx <- if (fr == 1) seq_along(ym) else unlist(lapply(split(seq_along(ym), ym), function(i) sample(i, round(fr * length(i)))))
    do.call(rbind, lapply(names(MODELS), function(k) {
      p <- if (fr == 1) loo_pred[[k]] else loo(Xm[idx, ], ym[idx], MODELS[[k]], seed = r)
      data.frame(fraction = fr, n = length(idx), repeat_id = r, model = k,
                 AUC = as.numeric(auc(roc(ym[idx], p, levels = c(0, 1), direction = "<", quiet = TRUE)))) })) })) }))
LCs <- aggregate(AUC ~ fraction + n + model, LC, function(a) c(mean = mean(a), sd = sd(a)))
LCs <- data.frame(LCs[, 1:3], AUC_mean = round(LCs$AUC[, "mean"], 3), AUC_sd = round(LCs$AUC[, "sd"], 3))
write.csv(LCs, file.path(rdir, "learning_curve.csv"), row.names = FALSE); write.csv(LC, file.path(rdir, "learning_curve_all_repeats.csv"), row.names = FALSE)
print(LCs, row.names = FALSE)

# ============================ STEP 5: feature importance =========================
cat("\n== Step 5: feature importance (full microarray models) ==\n")
gini <- importance(final[["Random forest"]], type = 2)[, 1]
coefs <- as.matrix(coef(final[["Logistic (L2)"]], s = "lambda.min"))[feat, 1]
FI <- data.frame(gene = feat, rf_gini = round(gini[feat], 3), lr_coef = round(coefs, 4))
FI$rf_rank <- rank(-FI$rf_gini); FI$lr_rank <- rank(-abs(FI$lr_coef)); FI$mean_rank <- (FI$rf_rank + FI$lr_rank) / 2
FI$lr_coef_pushes_toward <- ifelse(FI$lr_coef > 0, "AD", "control")      # conditional on all other genes (collinear panel)
FI$marginal_direction_in_AD <- ifelse(colMeans(Xm[ym == 1, feat]) > colMeans(Xm[ym == 0, feat]), "higher", "lower")[match(FI$gene, feat)]
FI$expected_marker <- FI$gene %in% c("VGF", "GFAP", "SST", "SNAP25")
FI <- FI[order(FI$mean_rank), ]; write.csv(FI, file.path(rdir, "feature_importance.csv"), row.names = FALSE)
print(head(FI, 15), row.names = FALSE)
cat("Expected markers (rank RF / LR):", paste(sprintf("%s %d/%d", FI$gene[FI$expected_marker], FI$rf_rank[FI$expected_marker], FI$lr_rank[FI$expected_marker]), collapse = "; "), "\n")

# ============================ FIGURES ===========================================
th <- theme_bw(base_size = 11) + theme(plot.background = element_rect(fill = "white", colour = NA), plot.title = element_text(face = "bold"))
rocdf <- function(y, p, set, k) { r <- roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
  data.frame(set = set, model = k, fpr = 1 - r$specificities, tpr = r$sensitivities, auc = as.numeric(auc(r))) }
R <- rbind(do.call(rbind, lapply(names(MODELS), function(k) rocdf(ym, loo_pred[[k]], "Microarray LOO-CV (n = 195)", k))),
           do.call(rbind, lapply(names(test_pred), function(k) { d <- test_pred[[k]]; tn <- sub(" (Logistic.*|Random.*)$", "", k)
             rocdf(d$condition, d$prob, ifelse(grepl("RNA", tn), "RNA-seq test (n = 22)", "SEA-AD test (n = 20)"), sub(".*\\) |.*pseudobulk ", "", k)) })))
R$set <- factor(R$set, levels = c("Microarray LOO-CV (n = 195)", "RNA-seq test (n = 22)", "SEA-AD test (n = 20)"))
R <- R[order(R$set, R$model, R$fpr, R$tpr), ]
lab <- unique(R[, c("set", "model", "auc")]); lab$y <- ifelse(lab$model == "Logistic (L2)", 0.16, 0.06)
g1 <- ggplot(R, aes(fpr, tpr, colour = model)) + geom_abline(linetype = "dashed", colour = "grey60") + geom_step(linewidth = 0.9, direction = "vh") +
  geom_text(data = lab, aes(x = 0.98, y = y, label = sprintf("%s AUC = %.2f", model, auc)), hjust = 1, size = 3.2, show.legend = FALSE) +
  facet_wrap(~set) + coord_equal() + scale_colour_manual(values = MCOL) +
  labs(title = "AD classification with the 31-gene panel", subtitle = "Trained on microarray; RNA-seq and SEA-AD predicted by the frozen model",
       x = "1 − specificity", y = "Sensitivity", colour = NULL) + th + theme(legend.position = "top")
ggsave(file.path(fdir, "fig1_roc_curves.pdf"), g1, width = 11, height = 4.6)
g2 <- ggplot(LCs, aes(n, AUC_mean, colour = model)) + geom_ribbon(aes(ymin = AUC_mean - AUC_sd, ymax = AUC_mean + AUC_sd, fill = model), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2) + scale_colour_manual(values = MCOL) + scale_fill_manual(values = MCOL, guide = "none") +
  scale_x_continuous(breaks = unique(LCs$n)) + coord_cartesian(ylim = c(0.5, 1)) +
  labs(title = "Learning curve (microarray)", subtitle = sprintf("LOO-CV AUC on stratified subsamples; mean ± SD over %d repeats", LC_REPEATS),
       x = "training samples", y = "LOO-CV AUC", colour = NULL) + th + theme(legend.position = "top")
ggsave(file.path(fdir, "fig2_learning_curve.pdf"), g2, width = 6.5, height = 4.8)
top <- head(FI, 15); top$gene <- factor(top$gene, levels = rev(top$gene))
g3a <- ggplot(top, aes(rf_gini, gene, fill = expected_marker)) + geom_col(width = 0.7) +
  scale_fill_manual(values = c(`TRUE` = "#eb6834", `FALSE` = "grey55"), labels = c(`TRUE` = "expected AD marker", `FALSE` = "other")) +
  labs(title = "Random forest", x = "Mean decrease in Gini", y = NULL, fill = NULL) + th + theme(legend.position = "bottom")
g3b <- ggplot(top, aes(lr_coef, gene, fill = lr_coef_pushes_toward)) + geom_col(width = 0.7) + geom_vline(xintercept = 0) +
  scale_fill_manual(values = c(AD = "#E63946", control = "#457B9D"), labels = c(AD = "pushes prediction toward AD", control = "toward control")) +
  labs(title = "L2 logistic regression", subtitle = "coefficients are conditional on the other (correlated) genes", x = "Coefficient (z-scored genes)", y = NULL, fill = NULL) + th + theme(legend.position = "bottom")
ggsave(file.path(fdir, "fig3_feature_importance.pdf"), (g3a | g3b) + plot_annotation(title = "Top 15 genes (ranked by mean of RF and LR ranks)",
       theme = theme(plot.title = element_text(face = "bold"))), width = 11, height = 5.5)
cm <- rbind(do.call(rbind, lapply(names(MODELS), function(k) { m <- L[L$model == k & grepl("sens", L$rule), ]; cbind(set = "Microarray LOO-CV", model = k, m[, c("TP", "FN", "TN", "FP")]) })),
            do.call(rbind, lapply(seq_len(nrow(Cx)), function(i) if (grepl("frozen", Cx$rule[i])) cbind(set = Cx$dataset[i], model = Cx$model[i], Cx[i, c("TP", "FN", "TN", "FP")]))))
cml <- do.call(rbind, lapply(seq_len(nrow(cm)), function(i) data.frame(set = cm$set[i], model = cm$model[i],
       truth = c("AD", "AD", "control", "control"), predicted = c("AD", "control", "control", "AD"), n = unlist(cm[i, c("TP", "FN", "TN", "FP")]))))
cml$correct <- cml$truth == cml$predicted; cml$set <- factor(cml$set, levels = unique(cml$set))
g4 <- ggplot(cml, aes(predicted, truth, fill = correct)) + geom_tile(colour = "white", linewidth = 1.5) + geom_text(aes(label = n), size = 5) +
  facet_grid(model ~ set) + scale_fill_manual(values = c(`TRUE` = "#cfe3f7", `FALSE` = "#f7d4d4"), guide = "none") +
  labs(title = "Confusion matrices at the frozen threshold (sensitivity >= 0.85 on microarray LOO-CV)", x = "Predicted", y = "True") + th
ggsave(file.path(fdir, "fig4_confusion_matrices.pdf"), g4, width = 10, height = 5.5)
cat("\nThresholds (from microarray LOO-CV):", paste(sprintf("%s %.3f", names(thr), thr), collapse = "; "), "\n\n")
print(sessioninfo::session_info())
sink()
