# =============================================================================
# 05a_validate_microarray.R — validation of Dataset 1 (GSE132903, limma)
# Re-derives the analysis from the raw GEO matrix and tests whether the result
# is correct and unbiased. Every check is written to results/checks_microarray.csv
# as PASS / WARN / FAIL with the number behind it.
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
proj <- ad_root()
vdir <- file.path(proj, "analysis", "validation")
suppressPackageStartupMessages({ library(limma); library(illuminaHumanv4.db); library(ggplot2) })
sink(file.path(vdir, "results", "05a_log.txt"), split = TRUE)

checks <- data.frame()
add <- function(id, check, value, verdict, note = "") {
  checks <<- rbind(checks, data.frame(dataset = "Microarray", id = id, check = check, value = value, verdict = verdict, note = note))
  cat(sprintf("[%s] %-4s %-48s %s  %s\n", verdict, id, check, value, note))
}
PADJ <- 0.05; LFC <- 0.5
MARK_DN <- c("VGF","SST","NPTX2","BDNF","SCG2","CHGB","PCSK1","CORT","CRH","RGS4")
MARK_UP <- c("GFAP","CD44","SERPINA3","CHI3L1","C4B","AQP4")

# ---- reload data exactly as script 01 ---------------------------------------
infile <- file.path(proj, "data", "GSE132903", "GSE132903_matrix.txt.gz")
con <- gzfile(infile, "rt"); hdr <- readLines(con, n = 200); close(con)
fld <- function(tag, n = 1) gsub('"', '', strsplit(hdr[grep(paste0("^", tag), hdr)][n], "\t")[[1]][-1])
gsm <- fld("!Sample_geo_accession"); titles <- fld("!Sample_title")
ch <- hdr[grep("^!Sample_characteristics_ch1", hdr)]
getch <- function(key) { l <- ch[grep(key, ch)][1]; sub(".*: ", "", gsub('"', '', strsplit(l, "\t")[[1]][-1])) }
meta <- data.frame(gsm = gsm, title_dx = ifelse(grepl("_AD_", titles), "AD", "ND"), geo_dx = getch("diagnosis"),
                   sex = getch("Sex"), age = as.numeric(sub("\\+", "", getch("expired_age"))))
expr <- as.matrix(read.delim(infile, comment.char = "!", row.names = 1, check.names = FALSE)[, gsm])
meta$condition <- factor(meta$title_dx, levels = c("ND", "AD"))
deg <- read.csv(file.path(proj, "analysis/dataset1_microarray/results/microarray_DEGs.csv"))

# probe -> gene (highest-mean probe), as in script 01
sym <- unname(mapIds(illuminaHumanv4.db, keys = rownames(expr), column = "SYMBOL", keytype = "PROBEID", multiVals = "first"))
ok <- !is.na(sym) & sym != ""
o <- order(sym[ok], -rowMeans(expr[ok, ])); E <- expr[ok, ][o, ]; s2 <- sym[ok][o]
G <- E[!duplicated(s2), ]; rownames(G) <- s2[!duplicated(s2)]
fitde <- function(X, design, coef) {
  tt <- topTable(eBayes(lmFit(X, design)), coef = coef, number = Inf, sort.by = "none")
  tt$padj <- p.adjust(tt$P.Value, "BH"); tt
}
nsig <- function(tt) sum(tt$padj < PADJ & abs(tt$logFC) > LFC)

# ---- V1 labels ----------------------------------------------------------------
agree <- mean(meta$title_dx == meta$geo_dx)
add("V1", "Title-derived labels match GEO 'diagnosis' field", sprintf("%.0f%% (AD %d / ND %d)", 100 * agree,
    sum(meta$condition == "AD"), sum(meta$condition == "ND")), ifelse(agree == 1, "PASS", "FAIL"))

# ---- V2 contrast direction: canonical AD markers ------------------------------
m <- deg[deg$gene %in% c(MARK_DN, MARK_UP), ]
m$expected <- ifelse(m$gene %in% MARK_UP, 1, -1); ok2 <- sum(sign(m$log2FC) == m$expected)
bp <- binom.test(ok2, nrow(m), 0.5)$p.value
add("V2", "Canonical AD markers in expected direction", sprintf("%d/%d (binomial p = %.1e)", ok2, nrow(m), bp),
    ifelse(ok2 / nrow(m) >= 0.8, "PASS", "FAIL"), "tests that log2FC = AD - control, not reversed")

# ---- V3 covariate balance -----------------------------------------------------
pa <- wilcox.test(age ~ condition, meta)$p.value; ps <- fisher.test(table(meta$sex, meta$condition))$p.value
add("V3a", "Age balanced between AD and ND", sprintf("median AD %.0f vs ND %.0f, p = %.2f",
    median(meta$age[meta$condition == "AD"]), median(meta$age[meta$condition == "ND"]), pa),
    ifelse(pa < 0.05, "WARN", "PASS"), "ages 90+ coded as 90")
add("V3b", "Sex balanced between AD and ND", sprintf("female AD %.0f%% vs ND %.0f%%, p = %.2f",
    100 * mean(meta$sex[meta$condition == "AD"] == "female"), 100 * mean(meta$sex[meta$condition == "ND"] == "female"), ps),
    ifelse(ps < 0.05, "WARN", "PASS"))

# ---- V4 covariate-adjusted model ---------------------------------------------
t0 <- fitde(G, model.matrix(~ condition, meta), "conditionAD")
t1 <- fitde(G, model.matrix(~ age + sex + condition, meta), "conditionAD")
r01 <- cor(t0$logFC, t1$logFC); j <- function(a, b) length(intersect(a, b)) / length(union(a, b))
s0 <- rownames(t0)[t0$padj < PADJ & abs(t0$logFC) > LFC]; s1 <- rownames(t1)[t1$padj < PADJ & abs(t1$logFC) > LFC]
add("V4", "Age + sex adjustment changes little", sprintf("r(log2FC) = %.3f; DEGs %d -> %d; Jaccard %.2f", r01, length(s0), length(s1), j(s0, s1)),
    ifelse(r01 > 0.95, "PASS", ifelse(r01 > 0.85, "WARN", "FAIL")))
stopifnot(abs(nsig(t0) - sum(deg$direction != "NS")) == 0)   # reproduces script 01 exactly

# ---- V5 p-value distribution ---------------------------------------------------
p <- t0$P.Value; pi0 <- min(1, mean(p > 0.5) / 0.5)
tail_ratio <- mean(p > 0.9) / mean(p > 0.5 & p <= 0.6)
add("V5", "p-value histogram: uniform tail + low-p spike", sprintf("pi0 = %.2f; tail ratio (p>0.9 vs 0.5-0.6) = %.2f", pi0, tail_ratio),
    ifelse(tail_ratio > 1.3 | tail_ratio < 0.7, "WARN", "PASS"), "pi0 = estimated fraction of unchanged genes")
write.csv(data.frame(dataset = "Microarray", p = p), file.path(vdir, "results", "pvals_microarray.csv"), row.names = FALSE)

# ---- V6 permutation null (FDR calibration) -------------------------------------
NPERM <- 200
nullsig <- replicate(NPERM, { mp <- meta; mp$condition <- sample(mp$condition)
  nsig(fitde(G, model.matrix(~ condition, mp), "conditionAD")) })
pe <- (sum(nullsig >= nsig(t0)) + 1) / (length(nullsig) + 1)
add("V6", "Label-permutation null: DEGs under shuffled labels",
    sprintf("observed %d; null median %.0f, max %d; %d/%d perms >= observed; empirical p = %.3f",
            nsig(t0), median(nullsig), max(nullsig), sum(nullsig >= nsig(t0)), length(nullsig), pe),
    if (pe > 0.05) "FAIL" else if (any(nullsig >= 0.5 * nsig(t0))) "WARN" else "PASS",
    "FAIL if empirical p > 0.05; WARN if any permutation reaches 50% of observed")
write.csv(data.frame(dataset = "Microarray", null_DEGs = nullsig, observed = nsig(t0)), file.path(vdir, "results", "permnull_microarray.csv"), row.names = FALSE)

# ---- V7 PCA: structure and outliers --------------------------------------------
pc <- prcomp(t(G[order(-apply(G, 1, var))[1:2000], ]), scale. = FALSE)
ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
assoc <- sapply(1:5, function(k) c(cond = wilcox.test(pc$x[, k] ~ meta$condition)$p.value,
                                    sex = wilcox.test(pc$x[, k] ~ meta$sex)$p.value,
                                    age = cor.test(pc$x[, k], meta$age, method = "spearman")$p.value))
d <- sqrt(rowSums(scale(pc$x[, 1:5])^2)); out <- names(d)[d > sqrt(qchisq(0.999, 5))]
add("V7a", "PCA: condition captured on leading PCs", sprintf("PC1 %.0f%% var; condition p on PC1-5: %s", ve[1],
    paste(sprintf("%.0e", assoc["cond", ]), collapse = ", ")), "INFO")
add("V7b", "PCA: sex/age not dominating PC1-2", sprintf("sex p PC1/PC2 = %.2g/%.2g; age p = %.2g/%.2g",
    assoc["sex", 1], assoc["sex", 2], assoc["age", 1], assoc["age", 2]),
    ifelse(min(assoc[c("sex", "age"), 1:2]) < 1e-3, "WARN", "PASS"))
add("V7c", "Outlier samples (Mahalanobis on PC1-5, p<0.001)", sprintf("%d: %s", length(out), paste(out, collapse = ", ")),
    ifelse(length(out) <= 2, "PASS", "WARN"))
pcd <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], condition = meta$condition, outlier = gsm %in% out)
g <- ggplot(pcd, aes(PC1, PC2, colour = condition)) + geom_point(size = 2, alpha = .8) +
  geom_point(data = pcd[pcd$outlier, ], shape = 21, size = 4, colour = "black") +
  scale_colour_manual(values = c(ND = "#457B9D", AD = "#E63946"), labels = c(ND = "control", AD = "AD")) +
  labs(title = "Microarray PCA (top 2,000 variable genes)", x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]),
       subtitle = "circled = Mahalanobis outlier (PC1-5, p < 0.001)") + theme_bw() + theme(plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(vdir, "figures", "pca_microarray.png"), g, width = 6.5, height = 5, dpi = 300, bg = "white")

# ---- V8 split-half reproducibility ---------------------------------------------
sh <- replicate(50, {
  idx <- unlist(lapply(split(seq_len(nrow(meta)), meta$condition), function(i) sample(i, length(i) %/% 2)))
  a <- fitde(G[, idx], model.matrix(~ condition, meta[idx, ]), "conditionAD")
  b <- fitde(G[, -idx], model.matrix(~ condition, meta[-idx, ]), "conditionAD")
  cor(a$logFC, b$logFC) })
add("V8", "Split-half reproducibility of log2FC (50 splits)", sprintf("median r = %.2f (range %.2f-%.2f)", median(sh), min(sh), max(sh)),
    ifelse(median(sh) > 0.5, "PASS", "WARN"), "internal consistency of the AD signal")

# ---- V9 probe-selection sensitivity --------------------------------------------
tp <- fitde(expr[ok, ], model.matrix(~ condition, meta), "conditionAD"); tp$gene <- sym[ok]
best <- tp[order(tp$P.Value), ]; best <- best[!duplicated(best$gene), ]
cmp <- merge(data.frame(gene = rownames(t0), a = t0$logFC), data.frame(gene = best$gene, b = best$logFC))
nb <- sum(p.adjust(best$P.Value, "BH") < PADJ & abs(best$logFC) > LFC)
add("V9", "Probe choice: highest-mean vs most-significant probe", sprintf("r = %.3f; DEGs %d vs %d (+%.0f%% with best-p probe)",
    cor(cmp$a, cmp$b), nsig(t0), nb, 100 * (nb / nsig(t0) - 1)), ifelse(cor(cmp$a, cmp$b) > 0.9, "PASS", "WARN"),
    "most-significant-probe selection inflates DEG counts; highest-mean is unbiased")

write.csv(checks, file.path(vdir, "results", "checks_microarray.csv"), row.names = FALSE)
sink()
