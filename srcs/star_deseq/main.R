# main.R

source("preprocessor.R")
source("deseq.R")
source("functools.R")



count_data <- BuildCountData()
meta_data <- BuildMetaData()


# Quality Control of Total Sample
deseq_total_samples <- RunDeseq(count_data, meta_data, ~Tumor)


pca_with_total_samples <- GetPCA(
    rld_mat = deseq_total_samples$rld_mat,
    meta_data = meta_data,
    intgroup = "Group3",
    shape = FALSE,
    pt_size = 3,
    label = FALSE,
    txt_size = 3,
    show_legend = TRUE
)

pca_with_total_samples_Tumor <- GetPCA(
    rld_mat = deseq_total_samples$rld_mat,
    meta_data = meta_data,
    intgroup = "Tumor",
    shape = TRUE,
    pt_size = 3,
    label = FALSE,
    txt_size = 3,
    show_legend = TRUE
)


cor_heatmap_with_total_samples <- GetCorHeatmap(
    cor_mat = deseq_total_samples$rld_cor_pearson,
    meta_data = meta_data[, 6, drop = FALSE],
    annotation = FALSE
)

# Tumor and Adjacent Tumor samples are mixed up let's dig into it
hcc_cdata <- count_data[, c(7:22)]
hcc_mdata <- meta_data[c(7:22), ]

deseq_tumor_samples <- RunDeseq(hcc_cdata, hcc_mdata, ~Tumor)

pca_with_tumor_samples <- GetPCA(
    rld_mat = deseq_tumor_samples$rld_mat,
    meta_data = hcc_mdata,
    intgroup = "Group3",
)

cor_heatmap_with_tumor_samples <- GetCorHeatmap(
    cor_mat = deseq_tumor_samples$rld_cor_spearman,
    meta_data = hcc_mdata[, 6, drop = FALSE]
)


# HCC 1x families with AT
deseq_hcc1x_AT_families <- RunDeseq(
    count_data = count_data[, -c(4:6, 15:22)],
    meta_data = meta_data[-c(4:6, 15:22), ],
    ~Tumor
)

pca_with_hcc1x_AT_families <- GetPCA(
    rld_mat = deseq_hcc1x_AT_families$rld_mat,
    meta_data = meta_data[-c(4:6, 15:22), ],
    intgroup = "Group3",
    pt_size = 3,
    show_legend = TRUE
)

cor_heatmap_with_hcc1x_AT_families <- GetCorHeatmap(
    cor_mat = deseq_hcc1x_AT_families$rld_cor_pearson,
    annotation = FALSE
)


# 1x and 3x seems separated and Tumor and Adjacent Tumor seems blended
# try with HCC Tumor 1x group
cdata_1x <- count_data[, -c(4:10, 15:22)]
mdata_1x <- meta_data[-c(4:10, 15:22), ]

deseq_hcc1x_families <- RunDeseq(cdata_1x, mdata_1x, ~Tumor)

pca_with_hcc1x_families <- GetPCA(
    rld_mat = deseq_hcc1x_families$rld_mat,
    meta_data = mdata_1x,
    intgroup = "Group3",
    pt_size = 5,
    show_legend = FALSE
)

cor_heatmap_with_hcc1x_families_pearson <- GetCorHeatmap(
    cor_mat = deseq_hcc1x_families$rld_cor_pearson,
    meta_data = mdata_1x[, 6, drop = FALSE]
)

cor_heatmap_with_hcc1x_families_spearman <- GetCorHeatmap(
    cor_mat = deseq_hcc1x_families$rld_cor_spearman,
    meta_data = mdata_1x[, 6, drop = FALSE]
)


# Control vs HCC 1x
deseq_con_vs_hcc1x <- RunDeseq(
    count_data = cdata_1x[, c(5, 6, 7, 11, 12, 13)],
    meta_data = mdata_1x[c(5, 6, 7, 11, 12, 13), ],
    ~Tumor
)

pca_with_con_vs_hcc1x <- GetPCA(
    rld_mat = deseq_con_vs_hcc1x$rld_mat,
    meta_data = mdata_1x[c(5, 6, 7, 11, 12, 13), ],
    intgroup = "Tumor"
)

cor_heatmap_with_con_vs_hcc1x <- GetCorHeatmap(
    cor_mat = deseq_con_vs_hcc1x$rld_cor,
    meta_data = mdata_1x[c(5, 6, 7, 11, 12, 13), 3, drop = FALSE]
)

ora_con_vs_hcc1x <- ora(deseq_con_vs_hcc1x, 2)
gsea_con_vs_hcc1x <- gsea(deseq_con_vs_hcc1x, 0.58)



# Control vs Nash
deseq_con_vs_nash <- RunDeseq(
    count_data = cdata_1x[, c(8:13)],
    meta_data = mdata_1x[c(8:13), ],
    ~Diet
)

ora_con_vs_nash <- ora(deseq_con_vs_nash, 2)
gsea_con_vs_nash <- gsea(deseq_con_vs_nash, 0.58)


# Control vs CCL4 1x
deseq_con_vs_ccl41x <- RunDeseq(
    count_data = cdata_1x[, c(1:3, 11:13)],
    meta_data = mdata_1x[c(1:3, 11:13), ],
    ~Treatment
)

ora_con_vs_ccl41x <- ora(deseq_con_vs_ccl41x, 2)
gsea_con_vs_ccl41x <- gsea(deseq_con_vs_ccl41x, 0.58)





# Nash vs HCC 1x
dds <- DESeq2::DESeqDataSetFromMatrix(
    cdata_1x[, c(5:10)],
    mdata_1x[c(5:10), ],
    ~Tumor
)


# CCL4 1x vs HCC 1x
dds <- DESeq2::DESeqDataSetFromMatrix(
    cdata_1x[, c(1:3, 5:7)],
    mdata_1x[c(1:3, 5:7), ],
    ~Tumor
)


# HCC 1x Tumor vs adj Tumor  -->  literally no significant genes exist
dds <- DESeq2::DESeqDataSetFromMatrix(
    count_data[, c(7:14)],
    meta_data[c(7:14), ],
    ~Tumor
)


# Control vs HCC 1x AT
deseq_con_vs_hcc1x_AT <- RunDeseq(
    count_data = count_data[, c(7:10, 26:28)],
    meta_data = meta_data[c(7:10, 26:28), ],
    ~Group2
)

ora_con_vs_hcc1x_AT <- ora(deseq_con_vs_hcc1x_AT, 2)


# let us define liver damage genes
liver_damaged_genes_1x <- intersect(
    intersect(
        intersect(deseq_con_vs_nash$sig_lfc_with_sym$gene_symbol,
                  deseq_con_vs_ccl41x$sig_lfc_with_sym$gene_symbol),
        deseq_con_vs_hcc1x$sig_lfc_with_sym$gene_symbol
    ),
    deseq_con_vs_hcc1x_AT$sig_lfc_with_sym$gene_symbol
)

nash_effect_genes_1x <- setdiff(
    deseq_con_vs_nash$sig_lfc_with_sym$gene_symbol,
    liver_damaged_genes_1x
)
nash_effect_1x_df <- deseq_con_vs_nash$sig_lfc_with_sym %>%
    dplyr::filter(gene_symbol %in% nash_effect_genes_1x)

ccl4_effect_genes_1x <- setdiff(
    deseq_con_vs_ccl41x$sig_lfc_with_sym$gene_symbol,
    liver_damaged_genes_1x
)
ccl4_effect_1x_df <- deseq_con_vs_ccl41x$sig_lfc_with_sym %>%
    dplyr::filter(gene_symbol %in% ccl4_effect_genes_1x)

tumor_effect_genes_1x <- setdiff(
    deseq_con_vs_hcc1x$sig_lfc_with_sym$gene_symbol,
    liver_damaged_genes_1x
)
tumor_effect_1x_df <- deseq_con_vs_hcc1x$sig_lfc_with_sym %>%
    dplyr::filter(gene_symbol %in% tumor_effect_genes_1x)
ora_tumor_effect_1x <- ora()


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




# HCC 3x families with AT
deseq_hcc3x_AT_families <- RunDeseq(
    count_data = count_data[, -c(1:3, 7:14, 17)],
    meta_data = meta_data[-c(1:3, 7:14, 17), ],
    ~Tumor
)

pca_with_hcc3x_AT_families <- GetPCA(
    rld_mat = deseq_hcc3x_AT_families$rld_mat,
    meta_data = meta_data[-c(1:3, 7:14, 17), ],
    intgroup = "Group3",
    pt_size = 3,
    show_legend = TRUE
)

cor_heatmap_with_hcc3x_AT_families <- GetCorHeatmap(
    cor_mat = deseq_hcc3x_AT_families$rld_cor_pearson,
    annotation = FALSE
)



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
deseq_con_vs_hcc3x <- RunDeseq(
    count_data = cdata_3x[, c(4, 5, 7, 11:13)],
    meta_data = mdata_3x[c(4, 5, 7, 11:13), ],
    ~Tumor
)

ora_con_vs_hcc3x <- ora(deseq_con_vs_hcc3x, 2)


# Control vs CCL4 3x
deseq_con_vs_ccl43x <- RunDeseq(
    count_data = cdata_3x[, c(1:3, 11:13)],
    meta_data = mdata_3x[c(1:3, 11:13), ],
    ~Treatment
)

ora_con_vs_ccl43x <- ora(deseq_con_vs_ccl43x, 2)


# Control vs HCC 3x AT
deseq_con_vs_hcc3x_AT <- RunDeseq(
    count_data = count_data[, c(15, 16, 18, 26:28)],
    meta_data = meta_data[c(15, 16, 18, 26:28), ],
    ~Group2
)

ora_con_vs_hcc3x_AT <- ora(deseq_con_vs_hcc3x_AT, 2)



# Nash vs HCC3x
dds <- DESeq2::DESeqDataSetFromMatrix(
    cdata_3x[, c(4:10)],
    mdata_3x[c(4:10), ],
    ~Tumor
)


# CCL4 3x vs HCC3x
dds <- DESeq2::DESeqDataSetFromMatrix(
    cdata_3x[, c(1:7)],
    mdata_3x[c(1:7), ],
    ~Tumor
)


# HCC 3x vs adj Tumor  -->  literally no significant genes exist
dds <- DESeq2::DESeqDataSetFromMatrix(
    count_data[, c(15:22)],
    meta_data[c(15:22), ],
    ~Tumor
)


# Control vs HCC 3x AT
dds <- DESeq2::DESeqDataSetFromMatrix(
    count_data[, c(15:18, 26:28)],
    meta_data[c(15:18, 26:28), ],
    ~Group2
)


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



