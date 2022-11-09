# deseq.R

source("idmapper.R")

library(tidyverse)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)


heat_colors <- rev(RColorBrewer::brewer.pal(11, 'RdBu'))


RunDeseq <- function(count_data, meta_data, design) {
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = count_data,
        colData = meta_data,
        design = design
    )
    
    dds_x <- DESeq2::DESeq(dds)
    
    rld_true <- DESeq2::rlog(dds_x, blind = TRUE)
    rld_false <- DESeq2::rlog(dds_x, blind = FALSE)
    
    rld_mat <- SummarizedExperiment::assay(rld_true)
    rld_cor_pearson <- stats::cor(rld_mat, method = "pearson")
    rld_cor_spearman <- stats::cor(rld_mat, method = "spearman")
    
    lfc_res <- DESeq2::lfcShrink(dds_x,
                                 coef = DESeq2::resultsNames(dds_x)[2],
                                 type = "apeglm")
    
    sig_lfc <- lfc_res %>%
        as.data.frame() %>%
        dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
        rownames_to_column(var = "ENSG")
    
    # some mixed samples might have empty sig_lfc
    try(expr = sig_lfc_with_sym <- ens2sym(sig_lfc))
    try(expr = sig_lfc_with_ent <- ens2entrez(sig_lfc))
    
    norm_cnt <- counts(dds_x, normalized = TRUE) %>%
        as.data.frame() %>%
        rownames_to_column(var = "ENSG")
    
    sig_norm_cnt <- norm_cnt %>%
        dplyr::filter(ENSG %in% sig_lfc$ENSG)
    
    final_res <- list(dds_x, rld_true, rld_false, rld_mat, 
                      rld_cor_pearson, rld_cor_spearman, 
                      lfc_res, sig_lfc, 
                      sig_lfc_with_sym, sig_lfc_with_ent,
                      norm_cnt, sig_norm_cnt)
    
    names(final_res) <- c("dds_x", "rld_true", "rld_false", "rld_mat", 
                          "rld_cor_pearson", "rld_cor_spearman",
                          "lfc_res", "sig_lfc", 
                          "sig_lfc_with_sym", "sig_lfc_with_ent",
                          "norm_cnt", "sig_norm_cnt")
    
    return(final_res)
}


GetPCA <- function(rld_mat, meta_data, intgroup,
                   shape = FALSE, pt_size = 3, 
                   label = FALSE, txt_size = 3,
                   show_legend = TRUE) {
    # GetPCA : get PCA result
    # 
    # [params]
    # @rld_mat : "assayed" result of "rlog" trnasformation of DESeq object
    # @meta_data : choose which meta_data to use
    # @intgroup : colname of meta_data interested in
    #
    # [return]
    # ggplot object
    
    # sanity check
    if (!(intgroup %in% colnames(meta_data))) {
        print("intgroup should be one of colname in meta_data!")
        stop()
    }
    
    pca <- stats::prcomp(t(rld_mat), scale = FALSE)
    percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    sd_ratio <- sqrt(percent_var[2] / percent_var[1])
    pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta_data)
    
    
    gg <- ggplot2::ggplot(pca_df, aes(PC1, PC2))
    
    
    if (shape) gg <- gg + 
        ggplot2::geom_point(aes(color = eval(as.name(intgroup)),
                                shape = eval(as.name(intgroup))), 
                            size = pt_size,
                            show.legend = show_legend)
    else gg <- gg + 
        ggplot2::geom_point(aes(color = eval(as.name(intgroup))), 
                            size = pt_size,
                            show.legend = show_legend)
    
    gg <- gg + 
        xlab(paste0("PC1 : ", percent_var[1], "%")) +
        ylab(paste0("PC2 : ", percent_var[2], "%")) +
        coord_fixed(ratio = sd_ratio) +
        theme_minimal() +
        theme(legend.key.width = unit(0.2, 'cm'))
    
    if (label) {
        gg <- gg + ggrepel::geom_text_repel(
            label = rownames(pca_df), size = txt_size
        )
    }

    return(gg)
}


GetCorHeatmap <- function(cor_mat, meta_data = NA, annotation = FALSE) {
    # GetCorHeatmap : get sample correlation heatmap
    #
    # [params]
    # @cor_mat : correlation matrix
    # @meta_data : meta_data containing columns for annotation
    #
    # [return]
    # pheatmap object
    
    if (annotation) annot_data <- meta_data
    else annot_data <- NA
    
    hm <- pheatmap::pheatmap(
        mat = cor_mat,
        #scale = "row",
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        color = heat_colors,
        annotation = annot_data
    )
    
    return(hm)
}
