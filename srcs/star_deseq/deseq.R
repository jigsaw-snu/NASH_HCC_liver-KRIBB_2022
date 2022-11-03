# deseq.R

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
    rld_cor <- stats::cor(rld_mat)
    
    lfc_res <- DESeq2::lfcShrink(dds_x,
                                 coef = DESeq2::resultsNames(dds_x)[2],
                                 type = "apeglm")
    
    sig_lfc <- lfc_res %>%
        as.data.frame() %>%
        dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
        rownames_to_column(var = "ENSG")
    
    norm_cnt <- counts(dds_x, normalized = TRUE) %>%
        as.data.frame() %>%
        rownames_to_column(var = "ENSG")
    
    sig_norm_cnt <- norm_cnt %>%
        dplyr::filter(ENSG %in% sig_lfc$ENSG)
    
    final_res <- list(dds_x, rld_true, rld_false, rld_mat, rld_cor, 
                      lfc_res, sig_lfc, norm_cnt, sig_norm_cnt)
    
    names(final_res) <- c("dds_x", "rld_true", "rld_false", "rld_mat", "rld_cor", 
                          "lfc_res", "sig_lfc", "norm_cnt", "sig_norm_cnt")
    
    return(final_res)
}


GetPCA <- function(rld_mat, meta_data, intgroup, label = FALSE, size = 3) {
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
    
    gg <- ggplot2::ggplot(pca_df, aes(PC1, PC2)) +
        ggplot2::geom_point(aes(color = eval(as.name(intgroup)))) +
        xlab(paste0("PC1 : ", percent_var[1], "%")) +
        ylab(paste0("PC2 : ", percent_var[2], "%")) +
        coord_fixed(ratio = sd_ratio) +
        theme_minimal()
    
    if (label) {
        gg <- gg + ggrepel::geom_text_repel(
            label = rownames(pca_df), size = size
        )
    }

    return(gg)
}


GetCorHeatmap <- function(rld_cor) {
    # GetCorHeatmap : get sample correlation heatmap
    #
    # [params]
    # @rld_cor : 
}