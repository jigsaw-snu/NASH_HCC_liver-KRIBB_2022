# main.R

source("preprocessor.R")
source("deseq.R")

count_data <- BuildCountData()
meta_data <- BuildMetaData()




# Quality Control of Total Sample
dds <- DESeq2::DESeqDataSetFromMatrix(count_data, meta_data, ~Tumor)
ddsX <- DESeq2::DESeq(dds)

rld_T <- DESeq2::rlog(ddsX, blind=TRUE)
rld_mat <- SummarizedExperiment::assay(rld_T)

#plotPCA(rld_T, intgroup = "Group3") + 
#    ggrepel::geom_text_repel(label = rownames(rld_T@colData), size = 3)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percent_var[2] / percent_var[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta_data)

ggplot2::ggplot(pca_df, aes(PC1, PC2)) + 
    geom_point(aes(color = eval(as.name("Group3"))), size=3) + 
    xlab(paste0("PC1 : ", percent_var[1], "%")) + 
    ylab(paste0("PC2 : ", percent_var[2], "%")) + 
    coord_fixed(ratio = sd_ratio) +
    theme_minimal() + 
    ggrepel::geom_text_repel(label = rownames(pca_df), size=2.5)

ggplot2::ggplot(pca_df, aes(PC1, PC2)) + 
    geom_point(aes(color = Tumor), size=3) + 
    xlab(paste0("PC1 : ", percent_var[1], "%")) + 
    ylab(paste0("PC2 : ", percent_var[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

rld_cor <- stats::cor(rld_mat)

meta_selected <- meta_data[, 5, drop = FALSE]
pheatmap::pheatmap(rld_cor, annotation = meta_selected, color = heat_colors)



# Tumor and Adjacent Tumor samples are mixed up let's dig into it
hcc_cdata <- count_data[, c(7:22)]
hcc_mdata <- meta_data[c(7:22), ]

dds <- DESeq2::DESeqDataSetFromMatrix(hcc_cdata, hcc_mdata, ~Tumor)
ddsX <- DESeq2::DESeq(dds)

rld_T <- DESeq2::rlog(ddsX, blind=TRUE)
rld_mat <- SummarizedExperiment::assay(rld_T)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percent_var[2] / percent_var[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], hcc_mdata)

ggplot2::ggplot(pca_df, aes(PC1, PC2)) + 
    geom_point(aes(color = Group3), size=3) + 
    xlab(paste0("PC1 : ", percent_var[1], "%")) + 
    ylab(paste0("PC2 : ", percent_var[2], "%")) + 
    coord_fixed(ratio = sd_ratio)

rld_cor <- stats::cor(rld_mat)

meta_selected <- meta_data[, 5, drop = FALSE]
pheatmap::pheatmap(rld_cor, annotation = hcc_mdata, color = heat_colors)





# 1x and 3x seems separated and Tumor and Adjacent Tumor seems blended
# try with HCC Tumor 1x group
cdata_1x <- count_data[, -c(4:10, 15:22)]
mdata_1x <- meta_data[-c(4:10, 15:22), ]

dds <- DESeq2::DESeqDataSetFromMatrix(cdata_1x, mdata_1x, ~Tumor)
ddsX <- DESeq2::DESeq(dds)

rld_T <- DESeq2::rlog(ddsX, blind=TRUE)
rld_mat <- SummarizedExperiment::assay(rld_T)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percent_var[2] / percent_var[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], mdata_1x)

ggplot2::ggplot(pca_df, aes(PC1, PC2)) + 
    geom_point(aes(color = Group3), size=3) + 
    xlab(paste0("PC1 : ", percent_var[1], "%")) + 
    ylab(paste0("PC2 : ", percent_var[2], "%")) + 
    coord_fixed(ratio = sd_ratio) +
    theme_minimal()

rld_cor <- stats::cor(rld_mat)

meta_selected <- meta_data[, 5, drop = FALSE]
pheatmap::pheatmap(rld_cor, annotation = meta_selected, color = heat_colors)

# Control vs HCC 1x
dds <- DESeq2::DESeqDataSetFromMatrix(
    cdata_1x[, c(5, 6, 7, 11, 12, 13)],
    mdata_1x[c(5, 6, 7, 11, 12, 13), ],
    ~Tumor
)
ddsX <- DESeq2::DESeq(dds)

test_rld <- DESeq2::rlog(ddsX, blind=T)
test_mat <- SummarizedExperiment::assay(test_rld)
test_cor <- stats::cor(test_mat)
pheatmap::pheatmap(test_cor, 
                   color = heat_colors,
                   scale = "row",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)

lfc_res <- DESeq2::lfcShrink(ddsX,
                             coef = DESeq2::resultsNames(ddsX)[2],
                             type = "apeglm")
sig_lfc <- lfc_res %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")

sig_lfc_with_id <- ens2sym(sig_lfc)
writexl::write_xlsx(sig_lfc_with_id, "temporary/con_vs_hcc1x_DEG.xlsx")
con_vs_hcc1x <- sig_lfc_with_id

norm_cnt <- counts(ddsX, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ENSG")

sig_norm_cnt <- norm_cnt %>%
    dplyr::filter(ENSG %in% sig_lfc$ENSG)

ego <- clusterProfiler::enrichGO(
    gene = sig_lfc %>% 
        dplyr::filter(log2FoldChange > 2) %>%
        dplyr::pull(ENSG),
    OrgDb = "org.Mm.eg.db",
    keyType = "ENSEMBL",
    ont = "BP",
    qvalueCutoff = 0.05,
    minGSSize = 20,
    maxGSSize = 1000
)
View(ego@result %>% dplyr::filter(p.adjust < 0.01 & qvalue < 0.01))

gsea_list <- sig_lfc %>% 
    dplyr::filter(log2FoldChange > 0.58) %>%
    dplyr::pull(log2FoldChange)

eg <- clusterProfiler::bitr(
    sig_lfc %>% 
        dplyr::filter(log2FoldChange > 0.58) %>%
        dplyr::pull(ENSG),
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db"
)

names(gsea_list) <- sig_lfc %>% 
    dplyr::filter(log2FoldChange > 0.58) %>%
    dplyr::pull(ENSG)
gsea_list <- sort(gsea_list, decreasing = TRUE)

gse <- clusterProfiler::gseGO(
    geneList = gsea_list,
    ont = "BP",
    OrgDb = "org.Mm.eg.db",
    keyType = "ENSEMBL",
    minGSSize = 10,
    maxGSSize = 1000,
    scoreType = "pos",
    pvalueCutoff = 0.05
)


# Control vs Nash
dds <- DESeq2::DESeqDataSetFromMatrix(
    cdata_1x[, c(8:13)],
    mdata_1x[c(8:13), ],
    ~Diet
)
ddsX <- DESeq2::DESeq(dds)  # no need to check sample cor (only 3 samples each)

lfc_res <- DESeq2::lfcShrink(ddsX,
                             coef = DESeq2::resultsNames(ddsX)[2],
                             type = "apeglm")
sig_lfc <- lfc_res %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")  # nothing

sig_lfc_with_id <- ens2sym(sig_lfc)
writexl::write_xlsx(sig_lfc_with_id, "temporary/con_vs_nash_DEG.xlsx")
con_vs_nash <- sig_lfc_with_id

# Control vs CCL4 1x
dds <- DESeq2::DESeqDataSetFromMatrix(
    cdata_1x[, c(1:3, 11:13)],
    mdata_1x[c(1:3, 11:13), ],
    ~Treatment
)
ddsX <- DESeq2::DESeq(dds)  # no need to check sample cor (only 3 samples each)

lfc_res <- DESeq2::lfcShrink(ddsX,
                             coef = DESeq2::resultsNames(ddsX)[2],
                             type = "apeglm")
sig_lfc <- lfc_res %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")  # nothing

sig_lfc_with_id <- ens2sym(sig_lfc)
writexl::write_xlsx(sig_lfc_with_id, "temporary/con_vs_ccl41x_DEG.xlsx")
con_vs_ccl41x <- sig_lfc_with_id

# Nash vs HCC 1x
dds <- DESeq2::DESeqDataSetFromMatrix(
    cdata_1x[, c(5:10)],
    mdata_1x[c(5:10), ],
    ~Tumor
)
ddsX <- DESeq2::DESeq(dds)

test_rld <- DESeq2::rlog(ddsX, blind=T)
test_mat <- SummarizedExperiment::assay(test_rld)
test_cor <- stats::cor(test_mat)
pheatmap::pheatmap(test_cor, 
                   color = heat_colors,
                   scale = "row",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)

lfc_res <- DESeq2::lfcShrink(ddsX,
                             coef = DESeq2::resultsNames(ddsX)[2],
                             type = "apeglm")
sig_lfc <- lfc_res %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")  # nothing

sig_lfc_with_id <- ens2sym(sig_lfc)
writexl::write_xlsx(sig_lfc_with_id, "temporary/nash_vs_hcc1x_DEG.xlsx")
nash_vs_hcc1x <- sig_lfc_with_id

# CCL4 1x vs HCC 1x
dds <- DESeq2::DESeqDataSetFromMatrix(
    cdata_1x[, c(1:3, 5:7)],
    mdata_1x[c(1:3, 5:7), ],
    ~Tumor
)
ddsX <- DESeq2::DESeq(dds)

test_rld <- DESeq2::rlog(ddsX, blind=T)
test_mat <- SummarizedExperiment::assay(test_rld)
test_cor <- stats::cor(test_mat)
pheatmap::pheatmap(test_cor, 
                   color = heat_colors,
                   scale = "row",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)

lfc_res <- DESeq2::lfcShrink(ddsX,
                             coef = DESeq2::resultsNames(ddsX)[2],
                             type = "apeglm")
sig_lfc <- lfc_res %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")  # nothing

sig_lfc_with_id <- ens2sym(sig_lfc)
writexl::write_xlsx(sig_lfc_with_id, "temporary/ccl41x_vs_hcc1x_DEG.xlsx")
ccl41x_vs_hcc1x <- sig_lfc_with_id

# HCC 1x Tumor vs adj Tumor  -->  literally no significant genes exist
dds <- DESeq2::DESeqDataSetFromMatrix(
    count_data[, c(7:14)],
    meta_data[c(7:14), ],
    ~Tumor
)
ddsX <- DESeq2::DESeq(dds)

test_rld <- DESeq2::rlog(ddsX, blind=T)
test_mat <- SummarizedExperiment::assay(test_rld)
test_cor <- stats::cor(test_mat)
pheatmap::pheatmap(test_cor, 
                   color = heat_colors,
                   scale = "row",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)

lfc_res <- DESeq2::lfcShrink(ddsX,
                             coef = DESeq2::resultsNames(ddsX)[2],
                             type = "apeglm")
sig_lfc <- lfc_res %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")  # nothing

# Control vs HCC 1x AT
dds <- DESeq2::DESeqDataSetFromMatrix(
    count_data[, c(7:10, 26:28)],
    meta_data[c(7:10, 26:28), ],
    ~Group2
)
ddsX <- DESeq2::DESeq(dds)

test_rld <- DESeq2::rlog(ddsX, blind=T)
test_mat <- SummarizedExperiment::assay(test_rld)
test_cor <- stats::cor(test_mat)
pheatmap::pheatmap(test_cor, 
                   color = heat_colors,
                   scale = "row",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)

lfc_res <- DESeq2::lfcShrink(ddsX,
                             coef = DESeq2::resultsNames(ddsX)[2],
                             type = "apeglm")
sig_lfc <- lfc_res %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")  # nothing

sig_lfc_with_id <- ens2sym(sig_lfc)
writexl::write_xlsx(sig_lfc_with_id, "temporary/con_vs_at1x_DEG.xlsx")
con_vs_at1x <- sig_lfc_with_id

# make comparing set
set1 <- setdiff(con_vs_hcc1x$gene_symbol, con_vs_at1x$gene_symbol)  # 651
set2 <- intersect(con_vs_hcc1x$gene_symbol, con_vs_at1x$gene_symbol)  # 2762
set3 <- setdiff(con_vs_at1x$gene_symbol, con_vs_hcc1x$gene_symbol)  # 639

set1 <- setdiff(
    con_vs_hcc1x %>%
        dplyr::filter(abs(log2FoldChange) > 3) %>%
        dplyr::pull(gene_symbol),
    con_vs_at1x %>%
        dplyr::filter(abs(log2FoldChange) > 3) %>%
        dplyr::pull(gene_symbol)
)
set2 <- intersect(
    con_vs_hcc1x %>%
        dplyr::filter(abs(log2FoldChange) > 3) %>%
        dplyr::pull(gene_symbol),
    con_vs_at1x %>%
        dplyr::filter(abs(log2FoldChange) > 3) %>%
        dplyr::pull(gene_symbol)
)
set3 <- setdiff(
    con_vs_at1x %>%
        dplyr::filter(abs(log2FoldChange) > 3) %>%
        dplyr::pull(gene_symbol),
    con_vs_hcc1x %>%
        dplyr::filter(abs(log2FoldChange) > 3) %>%
        dplyr::pull(gene_symbol)
)


set1_df <- con_vs_hcc1x %>% dplyr::filter(gene_symbol %in% set1)
set2_df <- merge(
    con_vs_hcc1x %>% 
        dplyr::filter(gene_symbol %in% set2) %>%
        dplyr::select(gene_symbol, log2FoldChange),
    con_vs_at1x %>%
        dplyr::filter(gene_symbol %in% set2) %>%
        dplyr::select(gene_symbol, log2FoldChange),
    by = "gene_symbol"
)
set2_df$direction <- ifelse(set2_df$log2FoldChange.x * set2_df$log2FoldChange.y > 0, 
                            'same', 
                            'opposite')
set2_df$large_gap <- ifelse(
    abs(set2_df$log2FoldChange.x / set2_df$log2FoldChange.y) > 2 | abs(set2_df$log2FoldChange.x / set2_df$log2FoldChange.y) < 1/2,
    'large',
    'small'
    )
set3_df <- con_vs_at1x %>% dplyr::filter(gene_symbol %in% set3)

writexl::write_xlsx(set1_df, "temporary/set1_1x.xlsx")
writexl::write_xlsx(set2_df, "temporary/set2_1x.xlsx")
writexl::write_xlsx(set3_df, "temporary/set3_1x.xlsx")


# try with HCC Tumor 3x group
cdata_3x <- count_data[, -c(1:3, 7:18)]
mdata_3x <- meta_data[-c(1:3, 7:18), ]

dds <- DESeq2::DESeqDataSetFromMatrix(cdata_3x, mdata_3x, ~Tumor)
ddsX <- DESeq2::DESeq(dds)

rld_T <- DESeq2::rlog(ddsX, blind=TRUE)
rld_mat <- SummarizedExperiment::assay(rld_T)

pca <- stats::prcomp(t(rld_mat), scale = FALSE)
percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percent_var[2] / percent_var[1])

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], mdata_3x)

ggplot2::ggplot(pca_df, aes(PC1, PC2)) + 
    geom_point(aes(color = Group3), size=3) + 
    xlab(paste0("PC1 : ", percent_var[1], "%")) + 
    ylab(paste0("PC2 : ", percent_var[2], "%")) + 
    coord_fixed(ratio = sd_ratio) +
    theme_minimal()

rld_cor <- stats::cor(rld_mat)

meta_selected <- meta_data[, 5, drop = FALSE]
pheatmap::pheatmap(rld_cor, annotation = meta_selected, color = heat_colors)

# Control vs HCC 3x
dds <- DESeq2::DESeqDataSetFromMatrix(
    cdata_3x[, c(4:7, 11:13)],
    mdata_3x[c(4:7, 11:13), ],
    ~Tumor
)
ddsX <- DESeq2::DESeq(dds)

test_rld <- DESeq2::rlog(ddsX, blind=T)
test_mat <- SummarizedExperiment::assay(test_rld)
test_cor <- stats::cor(test_mat)
pheatmap::pheatmap(test_cor, 
                   color = heat_colors,
                   scale = "row",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)

lfc_res <- DESeq2::lfcShrink(ddsX,
                             coef = DESeq2::resultsNames(ddsX)[2],
                             type = "apeglm")
sig_lfc <- lfc_res %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")

sig_lfc_with_id <- ens2sym(sig_lfc)
writexl::write_xlsx(sig_lfc_with_id, "temporary/con_vs_hcc3x_DEG.xlsx")
con_vs_hcc3x <- sig_lfc_with_id

# Control vs CCL4 3x

# Nash vs HCC3x
dds <- DESeq2::DESeqDataSetFromMatrix(
    cdata_3x[, c(4:10)],
    mdata_3x[c(4:10), ],
    ~Tumor
)
ddsX <- DESeq2::DESeq(dds)

test_rld <- DESeq2::rlog(ddsX, blind=T)
test_mat <- SummarizedExperiment::assay(test_rld)
test_cor <- stats::cor(test_mat)
pheatmap::pheatmap(test_cor, 
                   color = heat_colors,
                   scale = "row",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)

lfc_res <- DESeq2::lfcShrink(ddsX,
                             coef = DESeq2::resultsNames(ddsX)[2],
                             type = "apeglm")
sig_lfc <- lfc_res %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")  # nothing

sig_lfc_with_id <- ens2sym(sig_lfc)
writexl::write_xlsx(sig_lfc_with_id, "temporary/ccl41x_vs_hcc1x_DEG.xlsx")
ccl41x_vs_hcc1x <- sig_lfc_with_id

# CCL4 3x vs HCC3x
dds <- DESeq2::DESeqDataSetFromMatrix(
    cdata_3x[, c(1:7)],
    mdata_3x[c(1:7), ],
    ~Tumor
)
ddsX <- DESeq2::DESeq(dds)

test_rld <- DESeq2::rlog(ddsX, blind=T)
test_mat <- SummarizedExperiment::assay(test_rld)
test_cor <- stats::cor(test_mat)
pheatmap::pheatmap(test_cor, 
                   color = heat_colors,
                   scale = "row",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)

lfc_res <- DESeq2::lfcShrink(ddsX,
                             coef = DESeq2::resultsNames(ddsX)[2],
                             type = "apeglm")
sig_lfc <- lfc_res %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")  # nothing

sig_lfc_with_id <- ens2sym(sig_lfc)
writexl::write_xlsx(sig_lfc_with_id, "temporary/ccl41x_vs_hcc1x_DEG.xlsx")
ccl41x_vs_hcc1x <- sig_lfc_with_id

# HCC 3x vs adj Tumor  -->  literally no significant genes exist
dds <- DESeq2::DESeqDataSetFromMatrix(
    count_data[, c(15:22)],
    meta_data[c(15:22), ],
    ~Tumor
)
ddsX <- DESeq2::DESeq(dds)

test_rld <- DESeq2::rlog(ddsX, blind=T)
test_mat <- SummarizedExperiment::assay(test_rld)
test_cor <- stats::cor(test_mat)
pheatmap::pheatmap(test_cor, 
                   color = heat_colors,
                   scale = "row",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)

lfc_res <- DESeq2::lfcShrink(ddsX,
                             coef = DESeq2::resultsNames(ddsX)[2],
                             type = "apeglm")
sig_lfc <- lfc_res %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")  # nothing

# Control vs HCC 3x AT
dds <- DESeq2::DESeqDataSetFromMatrix(
    count_data[, c(15:18, 26:28)],
    meta_data[c(15:18, 26:28), ],
    ~Group2
)
ddsX <- DESeq2::DESeq(dds)

test_rld <- DESeq2::rlog(ddsX, blind=T)
test_mat <- SummarizedExperiment::assay(test_rld)
test_cor <- stats::cor(test_mat)
pheatmap::pheatmap(test_cor, 
                   color = heat_colors,
                   scale = "row",
                   cluster_rows = TRUE,
                   cluster_cols = TRUE)

lfc_res <- DESeq2::lfcShrink(ddsX,
                             coef = DESeq2::resultsNames(ddsX)[2],
                             type = "apeglm")
sig_lfc <- lfc_res %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
    rownames_to_column(var = "ENSG")  # nothing

sig_lfc_with_id <- ens2sym(sig_lfc)
writexl::write_xlsx(sig_lfc_with_id, "temporary/con_vs_at3x_DEG.xlsx")
con_vs_at3x <- sig_lfc_with_id

# make comparing set
set1 <- setdiff(con_vs_hcc3x$gene_symbol, con_vs_at3x$gene_symbol)  # 651
set2 <- intersect(con_vs_hcc3x$gene_symbol, con_vs_at3x$gene_symbol)  # 2762
set3 <- setdiff(con_vs_at3x$gene_symbol, con_vs_hcc3x$gene_symbol)  # 639

set1 <- setdiff(
    con_vs_hcc3x %>%
        dplyr::filter(abs(log2FoldChange) > 3) %>%
        dplyr::pull(gene_symbol),
    con_vs_at3x %>%
        dplyr::filter(abs(log2FoldChange) > 3) %>%
        dplyr::pull(gene_symbol)
)
set2 <- intersect(
    con_vs_hcc3x %>%
        dplyr::filter(abs(log2FoldChange) > 3) %>%
        dplyr::pull(gene_symbol),
    con_vs_at3x %>%
        dplyr::filter(abs(log2FoldChange) > 3) %>%
        dplyr::pull(gene_symbol)
)
set3 <- setdiff(
    con_vs_at3x %>%
        dplyr::filter(abs(log2FoldChange) > 3) %>%
        dplyr::pull(gene_symbol),
    con_vs_hcc3x %>%
        dplyr::filter(abs(log2FoldChange) > 3) %>%
        dplyr::pull(gene_symbol)
)


set1_df <- con_vs_hcc3x %>% dplyr::filter(gene_symbol %in% set1)
set2_df <- merge(
    con_vs_hcc3x %>% 
        dplyr::filter(gene_symbol %in% set2) %>%
        dplyr::select(gene_symbol, log2FoldChange),
    con_vs_at3x %>%
        dplyr::filter(gene_symbol %in% set2) %>%
        dplyr::select(gene_symbol, log2FoldChange),
    by = "gene_symbol"
)
set2_df$direction <- ifelse(set2_df$log2FoldChange.x * set2_df$log2FoldChange.y > 0, 
                            'same', 
                            'opposite')
set2_df$large_gap <- ifelse(
    abs(set2_df$log2FoldChange.x / set2_df$log2FoldChange.y) > 2 | abs(set2_df$log2FoldChange.x / set2_df$log2FoldChange.y) < 1/2,
    'large',
    'small'
)
set3_df <- con_vs_at3x %>% dplyr::filter(gene_symbol %in% set3)

writexl::write_xlsx(set1_df, "temporary/set1_3x.xlsx")
writexl::write_xlsx(set2_df, "temporary/set2_3x.xlsx")
writexl::write_xlsx(set3_df, "temporary/set3_3x.xlsx")



