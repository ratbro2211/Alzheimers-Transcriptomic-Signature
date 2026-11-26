# Install the R packages used by this repository.
cran <- c("WGCNA", "hdf5r", "glmnet", "randomForest", "pROC", "ggplot2", "ggrepel", "pheatmap", "patchwork", "UpSetR", "sessioninfo")
bioc <- c("limma", "illuminaHumanv4.db", "DESeq2", "apeglm")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE)
