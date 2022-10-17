# References
## https://wikis.utexas.edu/display/bioiteam/Testing+for+Differential+Expression
## https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto

################################################################################
################################################################################
################################################################################
library(org.Mm.eg.db)

txdb <- org.Mm.eg.db

# keytypes(txdb)  ## show list of available keytypes
key <- AnnotationDbi::keys(txdb, keytype = "ENSEMBLTRANS")  # key : ENSMUST
key2 <- AnnotationDbi::keys(txdb, keytype = "ENSEMBL")  # key : ENSMUSG

# columns(txdb)  ## show list of retriveable data using 'select'
tx2gene <- AnnotationDbi::select(txdb, 
                                 keys = key, 
                                 columns = "ENSEMBL", 
                                 keytype = "ENSEMBLTRANS")

ens2sym <- AnnotationDbi::select(txdb,
                                 keys = key2,
                                 columns = "SYMBOL",
                                 keytype = "ENSEMBL")
################################################################################
################################################################################
################################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

################################################################################
################################################################################
################################################################################
library(tximport)

# remember output of read.table is a dataframe --> you need to specify column
file_paths <- read.table('../../data/kallisto/kallisto_file_list',
                         sep = '\t',
                         header = FALSE)$V1

meta_data <- read.table('../../data/kallisto/SampleMetaData',
                        sep = '\t',
                        header = TRUE,
                        stringsAsFactors = TRUE)

# Backup
file_paths0 <- file_paths
meta_data0 <- meta_data

# Sample.Treatment --> Treatment (problem with read.table)
colnames(meta_data) <- c("Treatment", "Diet", "Tumor", "Group")

# relevel meta_data's columns
levels(meta_data$Treatment) <- c("saline", "ccl4_1x", "ccl4_3x")
levels(meta_data$Diet) <- c("ncd", "nashd")
levels(meta_data$Tumor) <- c("non_tumor", "adj_tumor", "tumor")
levels(meta_data$Group) <- c(
    "Saline_NCD_Non_Tumor", "Saline_NASHD_Non_Tumor",
    "CCL4_1x_NCD_Non_Tumor", "CCL4_3x_NCD_Non_Tumor",
    "CCL4_1x_NASHD_Adj_Tumor", "CCL4_3x_NASHD_Adj_Tumor",
    "CCL4_1x_NASHD_Tumor", "CCL4_3x_NASHD_Tumor"
)

# set names of file_paths (character vector)
# make sure order of column in meta_data is matched to order of file_paths
names(file_paths) <- rownames(meta_data)

txi <- tximport::tximport(file_paths,
                          type = "kallisto",
                          tx2gene = tx2gene,
                          ignoreTxVersion = TRUE)
################################################################################
################################################################################
################################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

################################################################################
################################################################################
################################################################################
library(DESeq2)
library(apeglm)
library(ashr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(readxl)
library(writexl)

# below function automatically convert string to factors
dds <- DESeq2::DESeqDataSetFromTximport(txi, 
                                        meta_data,
                                        ~Treatment + Diet + Tumor)

ddsX <- DESeq2::DESeq(dds)

rld <- DESeq2::rlog(ddsX, blind = TRUE)
rld_mat <- SummarizedExperiment::assay(rld)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta_data)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Group)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Treatment)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Diet)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Tumor)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)


rld_cor <- stats::cor(rld_mat)

heat_colors <- rev(RColorBrewer::brewer.pal(11, 'RdBu'))
pheatmap::pheatmap(mat = rld_cor, annotation = meta_data, color = heat_colors)

group_mdata <- meta_data[, 4, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = group_mdata, color = heat_colors)

treatment_mdata <- meta_data[, 1, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = treatment_mdata, color = heat_colors)

diet_mdata <- meta_data[, 2, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = diet_mdata, color = heat_colors)

tumor_mdata <- meta_data[, 3, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = tumor_mdata, color = heat_colors)
################################################################################
################################################################################
################################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

################################################################################
####################### Normal vs. HCC 1x Tumor ################################
################################################################################
file_paths <- file_paths0
meta_data <- meta_data0

file_paths <- file_paths[c(11, 12, 13, 14, 26, 27, 28)]
meta_data <- meta_data[c(11, 12, 13, 14, 26, 27, 28), ]

# Sample.Treatment --> Treatment (problem with read.table)
colnames(meta_data) <- c("Treatment", "Diet", "Tumor", "Group")

# relevel meta_data's columns
meta_data$Treatment <- factor(meta_data$Treatment)
levels(meta_data$Treatment) <- c("saline", "ccl4_1x")
levels(meta_data$Diet) <- c("ncd", "nashd")
meta_data$Tumor <- factor(meta_data$Tumor)
levels(meta_data$Tumor) <- c("non_tumor", "tumor")
meta_data$Group <- factor(meta_data$Group)
levels(meta_data$Group) <- c(
    "Saline_NCD_Non_Tumor", "CCL4_1x_NASHD_Tumor"
)

# set names of file_paths (character vector)
# make sure order of column in meta_data is matched to order of file_paths
names(file_paths) <- rownames(meta_data)

txi <- tximport::tximport(file_paths,
                          type = "kallisto",
                          tx2gene = tx2gene,
                          ignoreTxVersion = TRUE)

dds <- DESeq2::DESeqDataSetFromTximport(txi, 
                                        meta_data,
                                        ~Tumor)

ddsX <- DESeq2::DESeq(dds)

rld <- DESeq2::rlog(ddsX, blind = TRUE)
rld_mat <- SummarizedExperiment::assay(rld)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta_data)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Group)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Treatment)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Diet)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Tumor)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)


rld_cor <- stats::cor(rld_mat)

pheatmap::pheatmap(mat = rld_cor, annotation = meta_data, color = heat_colors)

group_mdata <- meta_data[, 4, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = group_mdata, color = heat_colors)

treatment_mdata <- meta_data[, 1, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = treatment_mdata, color = heat_colors)

diet_mdata <- meta_data[, 2, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = diet_mdata, color = heat_colors)

tumor_mdata <- meta_data[, 3, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = tumor_mdata, color = heat_colors)

lfcR <- DESeq2::lfcShrink(ddsX,
                          coef = resultsNames(ddsX)[2],
                          type = "apeglm")

sig_lfcR <- lfcR %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")

sig_lfcR$Gene <- rep('', nrow(sig_lfcR))
for (i in 1:nrow(sig_lfcR)) {
    sig_lfcR$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_lfcR$ENSG[i]]
}

sig_lfcR <- sig_lfcR[, c(7, 1, 2, 3, 4, 5, 6)]

hcc1x_con_up <- sig_lfcR[sig_lfcR$log2FoldChange > 0, ]
hcc1x_con_down <- sig_lfcR[sig_lfcR$log2FoldChange < 0, ]

writexl::write_xlsx(sig_lfcR,
                    path = "../../results/con_vs_hcc1x/siglfc.xlsx",
                    col_names = TRUE)

norm_cnt <- counts(ddsX, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ENSG")

sig_norm_cnt <- norm_cnt %>%
    dplyr::filter(ENSG %in% sig_lfcR$ENSG)

sig_norm_cnt$Gene <- rep('', nrow(sig_norm_cnt))
for (i in 1:nrow(sig_norm_cnt)) {
    sig_norm_cnt$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_norm_cnt$ENSG[i]]
}

sig_norm_cnt <- sig_norm_cnt[, -1]
sig_norm_cnt <- sig_norm_cnt %>%
    column_to_rownames(var = "Gene")

pheatmap::pheatmap(
    mat = sig_norm_cnt,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_colnames = TRUE,
    show_rownames = FALSE,
    scale = "row",
    color = heat_colors
)


# Top10 DEGs
sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = TRUE), ]
top10_df_up <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                          log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                          padj = sig_lfcR$padj[1:10])

sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = FALSE), ]
top10_df_down <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                            log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                            padj = sig_lfcR$padj[1:10])

top10_df <- rbind(top10_df_up, top10_df_down)
top10_df <- top10_df[order(top10_df$log2FoldChange), ]

gene_lev <- top10_df$Gene

top10_df$Gene <- factor(top10_df$Gene, levels = gene_lev)

ggplot(top10_df, aes(Gene, log2FoldChange)) +
    geom_col(aes(fill=-log(padj))) +
    coord_flip() + 
    scale_fill_continuous(low = 'blue', high = 'red')
################################################################################
################################################################################
################################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

################################################################################
####################### Normal vs. HCC 3x Tumor ################################
################################################################################
file_paths <- file_paths0
meta_data <- meta_data0

file_paths <- file_paths[c(19, 20, 21, 22, 26, 27, 28)]
meta_data <- meta_data[c(19, 20, 21, 22, 26, 27, 28), ]

# Sample.Treatment --> Treatment (problem with read.table)
colnames(meta_data) <- c("Treatment", "Diet", "Tumor", "Group")

# relevel meta_data's columns
meta_data$Treatment <- factor(meta_data$Treatment)
levels(meta_data$Treatment) <- c("saline", "ccl4_3x")
levels(meta_data$Diet) <- c("ncd", "nashd")
meta_data$Tumor <- factor(meta_data$Tumor)
levels(meta_data$Tumor) <- c("non_tumor", "tumor")
meta_data$Group <- factor(meta_data$Group)
levels(meta_data$Group) <- c(
    "Saline_NCD_Non_Tumor", "CCL4_3x_NASHD_Tumor"
)

# set names of file_paths (character vector)
# make sure order of column in meta_data is matched to order of file_paths
names(file_paths) <- rownames(meta_data)

txi <- tximport::tximport(file_paths,
                          type = "kallisto",
                          tx2gene = tx2gene,
                          ignoreTxVersion = TRUE)

dds <- DESeq2::DESeqDataSetFromTximport(txi, 
                                        meta_data,
                                        ~Tumor)

ddsX <- DESeq2::DESeq(dds)

rld <- DESeq2::rlog(ddsX, blind = TRUE)
rld_mat <- SummarizedExperiment::assay(rld)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta_data)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Group)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Treatment)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Diet)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Tumor)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)


rld_cor <- stats::cor(rld_mat)

pheatmap::pheatmap(mat = rld_cor, annotation = meta_data, color = heat_colors)

group_mdata <- meta_data[, 4, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = group_mdata, color = heat_colors)

treatment_mdata <- meta_data[, 1, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = treatment_mdata, color = heat_colors)

diet_mdata <- meta_data[, 2, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = diet_mdata, color = heat_colors)

tumor_mdata <- meta_data[, 3, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = tumor_mdata, color = heat_colors)

lfcR <- DESeq2::lfcShrink(ddsX,
                          coef = resultsNames(ddsX)[2],
                          type = "apeglm")

sig_lfcR <- lfcR %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")

sig_lfcR$Gene <- rep('', nrow(sig_lfcR))
for (i in 1:nrow(sig_lfcR)) {
    sig_lfcR$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_lfcR$ENSG[i]]
}

sig_lfcR <- sig_lfcR[, c(7, 1, 2, 3, 4, 5, 6)]

hcc3x_con_up <- sig_lfcR[sig_lfcR$log2FoldChange > 0, ]
hcc3x_con_down <- sig_lfcR[sig_lfcR$log2FoldChange < 0, ]

writexl::write_xlsx(sig_lfcR,
                    path = "../../results/con_vs_hcc3x/siglfc.xlsx",
                    col_names = TRUE)

norm_cnt <- counts(ddsX, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ENSG")

sig_norm_cnt <- norm_cnt %>%
    dplyr::filter(ENSG %in% sig_lfcR$ENSG)

sig_norm_cnt$Gene <- rep('', nrow(sig_norm_cnt))
for (i in 1:nrow(sig_norm_cnt)) {
    sig_norm_cnt$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_norm_cnt$ENSG[i]]
}

sig_norm_cnt <- sig_norm_cnt[, -1]
sig_norm_cnt <- sig_norm_cnt %>%
    column_to_rownames(var = "Gene")

pheatmap::pheatmap(
    mat = sig_norm_cnt,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_colnames = TRUE,
    show_rownames = FALSE,
    scale = "row",
    color = heat_colors
)

# Top10 DEGs
sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = TRUE), ]
top10_df_up <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                          log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                          padj = sig_lfcR$padj[1:10])

sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = FALSE), ]
top10_df_down <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                            log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                            padj = sig_lfcR$padj[1:10])

top10_df <- rbind(top10_df_up, top10_df_down)
top10_df <- top10_df[order(top10_df$log2FoldChange), ]

gene_lev <- top10_df$Gene

top10_df$Gene <- factor(top10_df$Gene, levels = gene_lev)

ggplot(top10_df, aes(Gene, log2FoldChange)) +
    geom_col(aes(fill=-log(padj))) +
    coord_flip() + 
    scale_fill_continuous(low = 'blue', high = 'red')
################################################################################
################################################################################
################################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

################################################################################
##################### Normal vs. HCC 1x Adj.Tumor ##############################
################################################################################
file_paths <- file_paths0
meta_data <- meta_data0

file_paths <- file_paths[c(7, 8, 9, 10, 26, 27, 28)]
meta_data <- meta_data[c(7, 8, 9, 10, 26, 27, 28), ]

# Sample.Treatment --> Treatment (problem with read.table)
colnames(meta_data) <- c("Treatment", "Diet", "Tumor", "Group")

# relevel meta_data's columns
meta_data$Treatment <- factor(meta_data$Treatment)
levels(meta_data$Treatment) <- c("saline", "ccl4_1x")
levels(meta_data$Diet) <- c("ncd", "nashd")
meta_data$Tumor <- factor(meta_data$Tumor)
levels(meta_data$Tumor) <- c("non_tumor", "adj_tumor")
meta_data$Group <- factor(meta_data$Group)
levels(meta_data$Group) <- c(
    "Saline_NCD_Non_Tumor", "CCL4_1x_NASHD_Adj_Tumor"
)

# set names of file_paths (character vector)
# make sure order of column in meta_data is matched to order of file_paths
names(file_paths) <- rownames(meta_data)

txi <- tximport::tximport(file_paths,
                          type = "kallisto",
                          tx2gene = tx2gene,
                          ignoreTxVersion = TRUE)

dds <- DESeq2::DESeqDataSetFromTximport(txi, 
                                        meta_data,
                                        ~Tumor)

ddsX <- DESeq2::DESeq(dds)

rld <- DESeq2::rlog(ddsX, blind = TRUE)
rld_mat <- SummarizedExperiment::assay(rld)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta_data)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Group)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Treatment)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Diet)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Tumor)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)


rld_cor <- stats::cor(rld_mat)

pheatmap::pheatmap(mat = rld_cor, annotation = meta_data, color = heat_colors)

group_mdata <- meta_data[, 4, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = group_mdata, color = heat_colors)

treatment_mdata <- meta_data[, 1, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = treatment_mdata, color = heat_colors)

diet_mdata <- meta_data[, 2, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = diet_mdata, color = heat_colors)

tumor_mdata <- meta_data[, 3, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = tumor_mdata, color = heat_colors)

lfcR <- DESeq2::lfcShrink(ddsX,
                          coef = resultsNames(ddsX)[2],
                          type = "apeglm")

sig_lfcR <- lfcR %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")

sig_lfcR$Gene <- rep('', nrow(sig_lfcR))
for (i in 1:nrow(sig_lfcR)) {
    sig_lfcR$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_lfcR$ENSG[i]]
}

sig_lfcR <- sig_lfcR[, c(7, 1, 2, 3, 4, 5, 6)]

writexl::write_xlsx(sig_lfcR,
                    path = "../../results/con_vs_adjhcc1x//siglfc.xlsx",
                    col_names = TRUE)

norm_cnt <- counts(ddsX, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ENSG")

sig_norm_cnt <- norm_cnt %>%
    dplyr::filter(ENSG %in% sig_lfcR$ENSG)

sig_norm_cnt$Gene <- rep('', nrow(sig_norm_cnt))
for (i in 1:nrow(sig_norm_cnt)) {
    sig_norm_cnt$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_norm_cnt$ENSG[i]]
}

sig_norm_cnt <- sig_norm_cnt[, -1]
sig_norm_cnt <- sig_norm_cnt %>%
    column_to_rownames(var = "Gene")

pheatmap::pheatmap(
    mat = sig_norm_cnt,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_colnames = TRUE,
    show_rownames = FALSE,
    scale = "row",
    color = heat_colors
)

# Top10 DEGs
sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = TRUE), ]
top10_df_up <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                          log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                          padj = sig_lfcR$padj[1:10])

sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = FALSE), ]
top10_df_down <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                            log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                            padj = sig_lfcR$padj[1:10])

top10_df <- rbind(top10_df_up, top10_df_down)
top10_df <- top10_df[order(top10_df$log2FoldChange), ]

gene_lev <- top10_df$Gene

top10_df$Gene <- factor(top10_df$Gene, levels = gene_lev)

ggplot(top10_df, aes(Gene, log2FoldChange)) +
    geom_col(aes(fill=-log(padj))) +
    coord_flip() + 
    scale_fill_continuous(low = 'blue', high = 'red')
################################################################################
################################################################################
################################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

################################################################################
##################### Normal vs. HCC 3x Adj.Tumor ##############################
################################################################################
file_paths <- file_paths0
meta_data <- meta_data0

file_paths <- file_paths[c(16, 17, 18, 26, 27, 28)]
meta_data <- meta_data[c(16, 17, 18, 26, 27, 28), ]

# Sample.Treatment --> Treatment (problem with read.table)
colnames(meta_data) <- c("Treatment", "Diet", "Tumor", "Group")

# relevel meta_data's columns
meta_data$Treatment <- factor(meta_data$Treatment)
levels(meta_data$Treatment) <- c("saline", "ccl4_3x")
levels(meta_data$Diet) <- c("ncd", "nashd")
meta_data$Tumor <- factor(meta_data$Tumor)
levels(meta_data$Tumor) <- c("non_tumor", "adj_tumor")
meta_data$Group <- factor(meta_data$Group)
levels(meta_data$Group) <- c(
    "Saline_NCD_Non_Tumor", "CCL4_3x_NASHD_Adj_Tumor"
)

# set names of file_paths (character vector)
# make sure order of column in meta_data is matched to order of file_paths
names(file_paths) <- rownames(meta_data)

txi <- tximport::tximport(file_paths,
                          type = "kallisto",
                          tx2gene = tx2gene,
                          ignoreTxVersion = TRUE)

dds <- DESeq2::DESeqDataSetFromTximport(txi, 
                                        meta_data,
                                        ~Tumor)

ddsX <- DESeq2::DESeq(dds)

rld <- DESeq2::rlog(ddsX, blind = TRUE)
rld_mat <- SummarizedExperiment::assay(rld)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta_data)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Group)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)
    #ggrepel::geom_text_repel(aes(label = rownames(pca_df)))

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Treatment)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Diet)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Tumor)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)


rld_cor <- stats::cor(rld_mat)

pheatmap::pheatmap(mat = rld_cor, annotation = meta_data, color = heat_colors)

group_mdata <- meta_data[, 4, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = group_mdata, color = heat_colors)

treatment_mdata <- meta_data[, 1, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = treatment_mdata, color = heat_colors)

diet_mdata <- meta_data[, 2, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = diet_mdata, color = heat_colors)

tumor_mdata <- meta_data[, 3, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = tumor_mdata, color = heat_colors)

lfcR <- DESeq2::lfcShrink(ddsX,
                          coef = resultsNames(ddsX)[2],
                          type = "apeglm")

sig_lfcR <- lfcR %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")

sig_lfcR$Gene <- rep('', nrow(sig_lfcR))
for (i in 1:nrow(sig_lfcR)) {
    sig_lfcR$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_lfcR$ENSG[i]]
}

sig_lfcR <- sig_lfcR[, c(7, 1, 2, 3, 4, 5, 6)]

writexl::write_xlsx(sig_lfcR,
                    path = "../../results/con_vs_adjhcc3x/siglfc.xlsx",
                    col_names = TRUE)

norm_cnt <- counts(ddsX, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ENSG")

sig_norm_cnt <- norm_cnt %>%
    dplyr::filter(ENSG %in% sig_lfcR$ENSG)

sig_norm_cnt$Gene <- rep('', nrow(sig_norm_cnt))
for (i in 1:nrow(sig_norm_cnt)) {
    sig_norm_cnt$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_norm_cnt$ENSG[i]]
}

sig_norm_cnt <- sig_norm_cnt[, -1]
sig_norm_cnt <- sig_norm_cnt %>%
    column_to_rownames(var = "Gene")

pheatmap::pheatmap(
    mat = sig_norm_cnt,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_colnames = TRUE,
    show_rownames = FALSE,
    scale = "row",
    color = heat_colors
)

# Top10 DEGs
sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = TRUE), ]
top10_df_up <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                          log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                          padj = sig_lfcR$padj[1:10])

sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = FALSE), ]
top10_df_down <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                            log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                            padj = sig_lfcR$padj[1:10])

top10_df <- rbind(top10_df_up, top10_df_down)
top10_df <- top10_df[order(top10_df$log2FoldChange), ]

gene_lev <- top10_df$Gene

top10_df$Gene <- factor(top10_df$Gene, levels = gene_lev)

ggplot(top10_df, aes(Gene, log2FoldChange)) +
    geom_col(aes(fill=-log(padj))) +
    coord_flip() + 
    scale_fill_continuous(low = 'blue', high = 'red')
################################################################################
################################################################################
################################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

################################################################################
######################### Normal vs. Nash Diet #################################
################################################################################
file_paths <- file_paths0
meta_data <- meta_data0

file_paths <- file_paths[c(23, 24, 25, 26, 27, 28)]
meta_data <- meta_data[c(23, 24, 25, 26, 27, 28), ]

# Sample.Treatment --> Treatment (problem with read.table)
colnames(meta_data) <- c("Treatment", "Diet", "Tumor", "Group")

# relevel meta_data's columns
meta_data$Treatment <- factor(meta_data$Treatment)
levels(meta_data$Treatment) <- c("saline")
levels(meta_data$Diet) <- c("ncd", "nashd")
meta_data$Tumor <- factor(meta_data$Tumor)
levels(meta_data$Tumor) <- c("non_tumor")
meta_data$Group <- factor(meta_data$Group)
levels(meta_data$Group) <- c(
    "Saline_NCD_Non_Tumor", "Saline_NASHD_Non_Tumor"
)

# set names of file_paths (character vector)
# make sure order of column in meta_data is matched to order of file_paths
names(file_paths) <- rownames(meta_data)

txi <- tximport::tximport(file_paths,
                          type = "kallisto",
                          tx2gene = tx2gene,
                          ignoreTxVersion = TRUE)

dds <- DESeq2::DESeqDataSetFromTximport(txi, 
                                        meta_data,
                                        ~Diet)

ddsX <- DESeq2::DESeq(dds)

rld <- DESeq2::rlog(ddsX, blind = TRUE)
rld_mat <- SummarizedExperiment::assay(rld)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta_data)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Group)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)
#ggrepel::geom_text_repel(aes(label = rownames(pca_df)))

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Treatment)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Diet)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Tumor)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)


rld_cor <- stats::cor(rld_mat)

pheatmap::pheatmap(mat = rld_cor, annotation = meta_data, color = heat_colors)

group_mdata <- meta_data[, 4, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = group_mdata, color = heat_colors)

treatment_mdata <- meta_data[, 1, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = treatment_mdata, color = heat_colors)

diet_mdata <- meta_data[, 2, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = diet_mdata, color = heat_colors)

tumor_mdata <- meta_data[, 3, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = tumor_mdata, color = heat_colors)

lfcR <- DESeq2::lfcShrink(ddsX,
                          coef = resultsNames(ddsX)[2],
                          type = "apeglm")

sig_lfcR <- lfcR %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")

sig_lfcR$Gene <- rep('', nrow(sig_lfcR))
for (i in 1:nrow(sig_lfcR)) {
    sig_lfcR$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_lfcR$ENSG[i]]
}

sig_lfcR <- sig_lfcR[, c(7, 1, 2, 3, 4, 5, 6)]

writexl::write_xlsx(sig_lfcR,
                    path = "../../results/con_vs_nash/siglfc.xlsx",
                    col_names = TRUE)

norm_cnt <- counts(ddsX, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ENSG")

sig_norm_cnt <- norm_cnt %>%
    dplyr::filter(ENSG %in% sig_lfcR$ENSG)

sig_norm_cnt$Gene <- rep('', nrow(sig_norm_cnt))
for (i in 1:nrow(sig_norm_cnt)) {
    sig_norm_cnt$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_norm_cnt$ENSG[i]]
}

sig_norm_cnt <- sig_norm_cnt[, -1]
sig_norm_cnt <- sig_norm_cnt %>%
    column_to_rownames(var = "Gene")

pheatmap::pheatmap(
    mat = sig_norm_cnt,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_colnames = TRUE,
    show_rownames = FALSE,
    scale = "row",
    color = heat_colors
)

# Top10 DEGs
sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = TRUE), ]
top10_df_up <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                          log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                          padj = sig_lfcR$padj[1:10])

sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = FALSE), ]
top10_df_down <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                            log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                            padj = sig_lfcR$padj[1:10])

top10_df <- rbind(top10_df_up, top10_df_down)
top10_df <- top10_df[order(top10_df$log2FoldChange), ]

gene_lev <- top10_df$Gene

top10_df$Gene <- factor(top10_df$Gene, levels = gene_lev)

ggplot(top10_df, aes(Gene, log2FoldChange)) +
    geom_col(aes(fill=-log(padj))) +
    coord_flip() + 
    scale_fill_continuous(low = 'blue', high = 'red')
################################################################################
################################################################################
################################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

################################################################################
######################### Normal vs. CCl4 1x ###################################
################################################################################
file_paths <- file_paths0
meta_data <- meta_data0

file_paths <- file_paths[c(1, 2, 3, 26, 27, 28)]
meta_data <- meta_data[c(1, 2, 3, 26, 27, 28), ]

# Sample.Treatment --> Treatment (problem with read.table)
colnames(meta_data) <- c("Treatment", "Diet", "Tumor", "Group")

# relevel meta_data's columns
meta_data$Treatment <- factor(meta_data$Treatment)
levels(meta_data$Treatment) <- c("saline", "ccl4_1x")
meta_data$Diet <- factor(meta_data$Diet)
levels(meta_data$Diet) <- c("ncd")
meta_data$Tumor <- factor(meta_data$Tumor)
levels(meta_data$Tumor) <- c("non_tumor")
meta_data$Group <- factor(meta_data$Group)
levels(meta_data$Group) <- c(
    "Saline_NCD_Non_Tumor", "CCL4_1x_NCD_Non_Tumor"
)

# set names of file_paths (character vector)
# make sure order of column in meta_data is matched to order of file_paths
names(file_paths) <- rownames(meta_data)

txi <- tximport::tximport(file_paths,
                          type = "kallisto",
                          tx2gene = tx2gene,
                          ignoreTxVersion = TRUE)

dds <- DESeq2::DESeqDataSetFromTximport(txi, 
                                        meta_data,
                                        ~Treatment)

ddsX <- DESeq2::DESeq(dds)

rld <- DESeq2::rlog(ddsX, blind = TRUE)
rld_mat <- SummarizedExperiment::assay(rld)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta_data)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Group)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)
#ggrepel::geom_text_repel(aes(label = rownames(pca_df)))

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Treatment)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Diet)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Tumor)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)


rld_cor <- stats::cor(rld_mat)

pheatmap::pheatmap(mat = rld_cor, annotation = meta_data, color = heat_colors)

group_mdata <- meta_data[, 4, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = group_mdata, color = heat_colors)

treatment_mdata <- meta_data[, 1, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = treatment_mdata, color = heat_colors)

diet_mdata <- meta_data[, 2, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = diet_mdata, color = heat_colors)

tumor_mdata <- meta_data[, 3, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = tumor_mdata, color = heat_colors)

lfcR <- DESeq2::lfcShrink(ddsX,
                          coef = resultsNames(ddsX)[2],
                          type = "apeglm")

sig_lfcR <- lfcR %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")

sig_lfcR$Gene <- rep('', nrow(sig_lfcR))
for (i in 1:nrow(sig_lfcR)) {
    sig_lfcR$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_lfcR$ENSG[i]]
}

sig_lfcR <- sig_lfcR[, c(7, 1, 2, 3, 4, 5, 6)]

writexl::write_xlsx(sig_lfcR,
                    path = "../../results/con_vs_ccl41x/siglfc.xlsx",
                    col_names = TRUE)

norm_cnt <- counts(ddsX, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ENSG")

sig_norm_cnt <- norm_cnt %>%
    dplyr::filter(ENSG %in% sig_lfcR$ENSG)

sig_norm_cnt$Gene <- rep('', nrow(sig_norm_cnt))
for (i in 1:nrow(sig_norm_cnt)) {
    sig_norm_cnt$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_norm_cnt$ENSG[i]]
}

sig_norm_cnt <- sig_norm_cnt[, -1]
sig_norm_cnt <- sig_norm_cnt %>%
    column_to_rownames(var = "Gene")

pheatmap::pheatmap(
    mat = sig_norm_cnt,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_colnames = TRUE,
    show_rownames = FALSE,
    scale = "row",
    color = heat_colors
)

# Top10 DEGs
sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = TRUE), ]
top10_df_up <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                          log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                          padj = sig_lfcR$padj[1:10])

sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = FALSE), ]
top10_df_down <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                            log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                            padj = sig_lfcR$padj[1:10])

top10_df <- rbind(top10_df_up, top10_df_down)
top10_df <- top10_df[order(top10_df$log2FoldChange), ]

gene_lev <- top10_df$Gene

top10_df$Gene <- factor(top10_df$Gene, levels = gene_lev)

ggplot(top10_df, aes(Gene, log2FoldChange)) +
    geom_col(aes(fill=-log(padj))) +
    coord_flip() + 
    scale_fill_continuous(low = 'blue', high = 'red')
################################################################################
################################################################################
################################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

################################################################################
######################### Normal vs. CCl4 3x ###################################
################################################################################
file_paths <- file_paths0
meta_data <- meta_data0

file_paths <- file_paths[c(4, 5, 6, 26, 27, 28)]
meta_data <- meta_data[c(4, 5, 6, 26, 27, 28), ]

# Sample.Treatment --> Treatment (problem with read.table)
colnames(meta_data) <- c("Treatment", "Diet", "Tumor", "Group")

# relevel meta_data's columns
meta_data$Treatment <- factor(meta_data$Treatment)
levels(meta_data$Treatment) <- c("saline", "ccl4_3x")
meta_data$Diet <- factor(meta_data$Diet)
levels(meta_data$Diet) <- c("ncd")
meta_data$Tumor <- factor(meta_data$Tumor)
levels(meta_data$Tumor) <- c("non_tumor")
meta_data$Group <- factor(meta_data$Group)
levels(meta_data$Group) <- c(
    "Saline_NCD_Non_Tumor", "CCL4_3x_NCD_Non_Tumor"
)

# set names of file_paths (character vector)
# make sure order of column in meta_data is matched to order of file_paths
names(file_paths) <- rownames(meta_data)

txi <- tximport::tximport(file_paths,
                          type = "kallisto",
                          tx2gene = tx2gene,
                          ignoreTxVersion = TRUE)

dds <- DESeq2::DESeqDataSetFromTximport(txi, 
                                        meta_data,
                                        ~Treatment)

ddsX <- DESeq2::DESeq(dds)

rld <- DESeq2::rlog(ddsX, blind = TRUE)
rld_mat <- SummarizedExperiment::assay(rld)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta_data)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Group)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)
#ggrepel::geom_text_repel(aes(label = rownames(pca_df)))

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Treatment)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Diet)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Tumor)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)


rld_cor <- stats::cor(rld_mat)

pheatmap::pheatmap(mat = rld_cor, annotation = meta_data, color = heat_colors)

group_mdata <- meta_data[, 4, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = group_mdata, color = heat_colors)

treatment_mdata <- meta_data[, 1, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = treatment_mdata, color = heat_colors)

diet_mdata <- meta_data[, 2, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = diet_mdata, color = heat_colors)

tumor_mdata <- meta_data[, 3, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = tumor_mdata, color = heat_colors)

lfcR <- DESeq2::lfcShrink(ddsX,
                          coef = resultsNames(ddsX)[2],
                          type = "apeglm")

sig_lfcR <- lfcR %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")

sig_lfcR$Gene <- rep('', nrow(sig_lfcR))
for (i in 1:nrow(sig_lfcR)) {
    sig_lfcR$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_lfcR$ENSG[i]]
}

sig_lfcR <- sig_lfcR[, c(7, 1, 2, 3, 4, 5, 6)]

writexl::write_xlsx(sig_lfcR,
                    path = "../../results/con_vs_ccl43x/siglfc.xlsx",
                    col_names = TRUE)

norm_cnt <- counts(ddsX, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ENSG")

sig_norm_cnt <- norm_cnt %>%
    dplyr::filter(ENSG %in% sig_lfcR$ENSG)

sig_norm_cnt$Gene <- rep('', nrow(sig_norm_cnt))
for (i in 1:nrow(sig_norm_cnt)) {
    sig_norm_cnt$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_norm_cnt$ENSG[i]]
}

sig_norm_cnt <- sig_norm_cnt[, -1]
sig_norm_cnt <- sig_norm_cnt %>%
    column_to_rownames(var = "Gene")

pheatmap::pheatmap(
    mat = sig_norm_cnt,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_colnames = TRUE,
    show_rownames = FALSE,
    scale = "row",
    color = heat_colors
)

# Top10 DEGs
sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = TRUE), ]
top10_df_up <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                          log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                          padj = sig_lfcR$padj[1:10])

sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = FALSE), ]
top10_df_down <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                            log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                            padj = sig_lfcR$padj[1:10])

top10_df <- rbind(top10_df_up, top10_df_down)
top10_df <- top10_df[order(top10_df$log2FoldChange), ]

gene_lev <- top10_df$Gene

top10_df$Gene <- factor(top10_df$Gene, levels = gene_lev)

ggplot(top10_df, aes(Gene, log2FoldChange)) +
    geom_col(aes(fill=-log(padj))) +
    coord_flip() + 
    scale_fill_continuous(low = 'blue', high = 'red')
################################################################################
################################################################################
################################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

################################################################################
################## HCC 1x Adj.Tumor vs. HCC 1x Tumor ###########################
################################################################################
file_paths <- file_paths0
meta_data <- meta_data0

file_paths <- file_paths[c(7, 8, 9, 10, 11, 12, 13, 14)]
meta_data <- meta_data[c(7, 8, 9, 10, 11, 12, 13, 14), ]

# Sample.Treatment --> Treatment (problem with read.table)
colnames(meta_data) <- c("Treatment", "Diet", "Tumor", "Group")

# relevel meta_data's columns
meta_data$Treatment <- factor(meta_data$Treatment)
levels(meta_data$Treatment) <- c("ccl4_1x")
meta_data$Diet <- factor(meta_data$Diet)
levels(meta_data$Diet) <- c("nashd")
meta_data$Tumor <- factor(meta_data$Tumor)
levels(meta_data$Tumor) <- c("adj_tumor", "tumor")
meta_data$Group <- factor(meta_data$Group)
levels(meta_data$Group) <- c(
    "CCL4_1x_NASHD_Adj_Tumor", "CCL4_1x_NASHD_Tumor"
)

# set names of file_paths (character vector)
# make sure order of column in meta_data is matched to order of file_paths
names(file_paths) <- rownames(meta_data)

txi <- tximport::tximport(file_paths,
                          type = "kallisto",
                          tx2gene = tx2gene,
                          ignoreTxVersion = TRUE)

dds <- DESeq2::DESeqDataSetFromTximport(txi, 
                                        meta_data,
                                        ~Tumor)

ddsX <- DESeq2::DESeq(dds)

rld <- DESeq2::rlog(ddsX, blind = TRUE)
rld_mat <- SummarizedExperiment::assay(rld)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta_data)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Group)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)
#ggrepel::geom_text_repel(aes(label = rownames(pca_df)))

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Treatment)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Diet)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Tumor)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)


rld_cor <- stats::cor(rld_mat)

pheatmap::pheatmap(mat = rld_cor, annotation = meta_data, color = heat_colors)

group_mdata <- meta_data[, 4, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = group_mdata, color = heat_colors)

treatment_mdata <- meta_data[, 1, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = treatment_mdata, color = heat_colors)

diet_mdata <- meta_data[, 2, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = diet_mdata, color = heat_colors)

tumor_mdata <- meta_data[, 3, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = tumor_mdata, color = heat_colors)

lfcR <- DESeq2::lfcShrink(ddsX,
                          coef = resultsNames(ddsX)[2],
                          type = "apeglm")

sig_lfcR <- lfcR %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")

sig_lfcR$Gene <- rep('', nrow(sig_lfcR))
for (i in 1:nrow(sig_lfcR)) {
    sig_lfcR$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_lfcR$ENSG[i]]
}

sig_lfcR <- sig_lfcR[, c(7, 1, 2, 3, 4, 5, 6)]

writexl::write_xlsx(sig_lfcR,
                    path = "../../results/adjhcc1x_vs_hcc1x/siglfc.xlsx",
                    col_names = TRUE)

norm_cnt <- counts(ddsX, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ENSG")

sig_norm_cnt <- norm_cnt %>%
    dplyr::filter(ENSG %in% sig_lfcR$ENSG)

sig_norm_cnt$Gene <- rep('', nrow(sig_norm_cnt))
for (i in 1:nrow(sig_norm_cnt)) {
    sig_norm_cnt$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_norm_cnt$ENSG[i]]
}

sig_norm_cnt <- sig_norm_cnt[, -1]
sig_norm_cnt <- sig_norm_cnt %>%
    column_to_rownames(var = "Gene")

pheatmap::pheatmap(
    mat = sig_norm_cnt,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_colnames = TRUE,
    show_rownames = FALSE,
    scale = "row",
    color = heat_colors
)

# Top10 DEGs
sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = TRUE), ]
top10_df_up <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                          log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                          padj = sig_lfcR$padj[1:10])

sig_lfcR <- sig_lfcR[order(sig_lfcR$log2FoldChange, decreasing = FALSE), ]
top10_df_down <- data.frame(Gene = sig_lfcR$Gene[1:10], 
                            log2FoldChange = sig_lfcR$log2FoldChange[1:10],
                            padj = sig_lfcR$padj[1:10])

top10_df <- rbind(top10_df_up, top10_df_down)
top10_df <- top10_df[order(top10_df$log2FoldChange), ]

gene_lev <- top10_df$Gene

top10_df$Gene <- factor(top10_df$Gene, levels = gene_lev)

ggplot(top10_df, aes(Gene, log2FoldChange)) +
    geom_col(aes(fill=-log(padj))) +
    coord_flip() + 
    scale_fill_continuous(low = 'blue', high = 'red')
################################################################################
################################################################################
################################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

################################################################################
################## HCC 3x Adj.Tumor vs. HCC 3x Tumor ###########################
################################################################################
file_paths <- file_paths0
meta_data <- meta_data0

file_paths <- file_paths[c(16, 17, 18, 20, 21, 22)]
meta_data <- meta_data[c(16, 17, 18, 20, 21, 22), ]

# Sample.Treatment --> Treatment (problem with read.table)
colnames(meta_data) <- c("Treatment", "Diet", "Tumor", "Group")

# relevel meta_data's columns
meta_data$Treatment <- factor(meta_data$Treatment)
levels(meta_data$Treatment) <- c("ccl4_3x")
meta_data$Diet <- factor(meta_data$Diet)
levels(meta_data$Diet) <- c("nashd")
meta_data$Tumor <- factor(meta_data$Tumor)
levels(meta_data$Tumor) <- c("adj_tumor", "tumor")
meta_data$Group <- factor(meta_data$Group)
levels(meta_data$Group) <- c(
    "CCL4_3x_NASHD_Adj_Tumor", "CCL4_3x_NASHD_Tumor"
)

# set names of file_paths (character vector)
# make sure order of column in meta_data is matched to order of file_paths
names(file_paths) <- rownames(meta_data)

txi <- tximport::tximport(file_paths,
                          type = "kallisto",
                          tx2gene = tx2gene,
                          ignoreTxVersion = TRUE)

dds <- DESeq2::DESeqDataSetFromTximport(txi, 
                                        meta_data,
                                        ~Tumor)

ddsX <- DESeq2::DESeq(dds)

rld <- DESeq2::rlog(ddsX, blind = TRUE)
rld_mat <- SummarizedExperiment::assay(rld)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta_data)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Group)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)
#ggrepel::geom_text_repel(aes(label = rownames(pca_df)))

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Treatment)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Diet)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = Tumor)) + 
    xlab(paste0("PC1 : ", percentVar[1], "%")) + 
    ylab(paste0("PC2 : ", percentVar[2], "%")) + 
    coord_fixed(ratio = sd_ratio)


rld_cor <- stats::cor(rld_mat)

pheatmap::pheatmap(mat = rld_cor, annotation = meta_data, color = heat_colors)

group_mdata <- meta_data[, 4, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = group_mdata, color = heat_colors)

treatment_mdata <- meta_data[, 1, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = treatment_mdata, color = heat_colors)

diet_mdata <- meta_data[, 2, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = diet_mdata, color = heat_colors)

tumor_mdata <- meta_data[, 3, drop = FALSE]
pheatmap::pheatmap(mat = rld_cor, annotation = tumor_mdata, color = heat_colors)

lfcR <- DESeq2::lfcShrink(ddsX,
                          coef = resultsNames(ddsX)[2],
                          type = "apeglm")

sig_lfcR <- lfcR %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")

sig_lfcR$Gene <- rep('', nrow(sig_lfcR))
for (i in 1:nrow(sig_lfcR)) {
    sig_lfcR$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_lfcR$ENSG[i]]
}

sig_lfcR <- sig_lfcR[, c(7, 1, 2, 3, 4, 5, 6)]

writexl::write_xlsx(sig_lfcR,
                    path = "../../results/adjhcc3x_vs_hcc3x/siglfc.xlsx",
                    col_names = TRUE)

norm_cnt <- counts(ddsX, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ENSG")

sig_norm_cnt <- norm_cnt %>%
    dplyr::filter(ENSG %in% sig_lfcR$ENSG)

sig_norm_cnt$Gene <- rep('', nrow(sig_norm_cnt))
for (i in 1:nrow(sig_norm_cnt)) {
    sig_norm_cnt$Gene[i] <- ens2sym$SYMBOL[ens2sym$ENSEMBL == sig_norm_cnt$ENSG[i]]
}

sig_norm_cnt <- sig_norm_cnt[, -1]
sig_norm_cnt <- sig_norm_cnt %>%
    column_to_rownames(var = "Gene")

pheatmap::pheatmap(
    mat = sig_norm_cnt,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_colnames = TRUE,
    show_rownames = FALSE,
    scale = "row",
    color = heat_colors
)
