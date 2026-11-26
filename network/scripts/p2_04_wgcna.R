# =============================================================================
# p2_04_wgcna.R
# Unsupervised WGCNA on the microarray (GSE132903, n = 195, best-powered dataset):
#   - genes: 5,000 most variable genes  UNION  the 33-gene set
#   - signed network, soft-threshold chosen by scale-free fit (R^2 >= 0.8)
#   - question: do the 3 hand-defined modules (Glial / Neuropeptide / Synaptic)
#     re-emerge as separate data-driven co-expression modules?
#   - module eigengene association with AD, and module membership (kME) of the 33 genes
# Run from network/ (microarray matrix from p2_00)
# =============================================================================
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE); enableWGCNAThreads(4)

g33 <- strsplit("VGF SST SCG2 CARTPT CHGB VSNL1 TAC1 SNAP25 NPTX2 NEFL GABRA1 NRGN GABRG2 SYT4 CCK NEFM SHANK2 GRIN2A NECAB1 LRRTM1 CBLN2 CAMK2A SYN2 HOMER1 NPY SCN1A PENK NLGN1 CADM1 VIM TREM2 GJA1 GFAP", " ")[[1]]
mod_pre <- setNames(rep("Synaptic", length(g33)), g33)
mod_pre[c("GFAP","VIM","GJA1","TREM2")] <- "Glial"
mod_pre[c("VGF","SST","NPY","CARTPT","TAC1","PENK","CHGB","SCG2","VSNL1")] <- "Neuropeptide"

X   <- read.csv(gzfile("data/expr_ARR.csv.gz"), row.names = 1, check.names = FALSE)
lab <- read.csv("data/labels_ARR.csv", row.names = 1)
lab <- lab[colnames(X), "group"]
X   <- X[rowMeans(X) > quantile(rowMeans(X), 0.25), ]                 # drop lowest-expressed quartile
v   <- apply(X, 1, var)
keep <- union(names(sort(v, decreasing = TRUE))[1:5000], intersect(g33, rownames(X)))
datExpr <- t(X[keep, ])
cat("WGCNA input:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")
stopifnot(goodSamplesGenes(datExpr)$allOK)

sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 20, 2)), networkType = "signed", verbose = 0)
fit <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
power <- sft$fitIndices$Power[which(fit >= 0.80)[1]]
if (is.na(power)) power <- 12
cat("Soft-threshold power:", power, " (scale-free R^2 =", round(fit[sft$fitIndices$Power == power], 2), ")\n")
write.csv(cbind(sft$fitIndices, signed_R2 = fit), "results/wgcna_soft_threshold.csv", row.names = FALSE)

net <- blockwiseModules(datExpr, power = power, networkType = "signed", TOMType = "signed",
                        minModuleSize = 30, mergeCutHeight = 0.25, maxBlockSize = 6000,
                        numericLabels = FALSE, pamRespectsDendro = FALSE, randomSeed = 1, verbose = 0)
cols <- net$colors
cat("Modules:", length(unique(cols)), "\n"); print(sort(table(cols), decreasing = TRUE))

ME  <- orderMEs(net$MEs)
ad  <- as.numeric(lab == "AD")
mt  <- data.frame(module = sub("^ME", "", colnames(ME)),
                  n_genes = as.integer(table(cols)[sub("^ME", "", colnames(ME))]),
                  r_AD = as.numeric(cor(ME, ad)),
                  p = corPvalueStudent(as.numeric(cor(ME, ad)), nrow(datExpr)))
mt$fdr <- p.adjust(mt$p, "BH")
mt <- mt[order(mt$p), ]
write.csv(mt, "results/wgcna_module_trait.csv", row.names = FALSE)
cat("\nModule eigengene vs AD:\n"); print(mt, row.names = FALSE, digits = 3)

kME <- signedKME(datExpr, ME)
gs  <- intersect(g33, colnames(datExpr))
tab <- data.frame(gene = gs, predefined_module = mod_pre[gs], wgcna_module = cols[match(gs, colnames(datExpr))])
tab$kME_own <- mapply(function(g, m) kME[g, paste0("kME", m)], tab$gene, tab$wgcna_module)
tab <- tab[order(tab$predefined_module, tab$wgcna_module), ]
write.csv(tab, "results/wgcna_33gene_membership.csv", row.names = FALSE)
cat("\nPredefined module x WGCNA module (33 genes):\n"); print(table(tab$predefined_module, tab$wgcna_module))
print(tab, row.names = FALSE, digits = 2)

# Top hub genes (highest kME) of the modules that contain the 33 genes, for context
hubs <- do.call(rbind, lapply(unique(tab$wgcna_module), function(m) {
  k <- kME[cols == m, paste0("kME", m)]; names(k) <- colnames(datExpr)[cols == m]
  k <- sort(k, decreasing = TRUE)[1:15]
  data.frame(module = m, gene = names(k), kME = round(k, 3), in_33 = names(k) %in% g33)
}))
write.csv(hubs, "results/wgcna_module_hubs.csv", row.names = FALSE)
write.csv(data.frame(gene = colnames(datExpr), module = cols), "results/wgcna_all_gene_modules.csv", row.names = FALSE)
saveRDS(list(net = net, power = power), "results/wgcna_net.rds")

pdf("figures/wgcna_dendrogram.pdf", width = 10, height = 4.5)
plotDendroAndColors(net$dendrograms[[1]], cols[net$blockGenes[[1]]], "Module",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = sprintf("WGCNA (signed, power %d) - GSE132903 MTG, n = 195", power))
dev.off()
