# deseq.R

library(DESeq2)
library(apeglm)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)


heat_colors <- RColorBrewer::brewer.pal(11, 'RdBu')


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


GetPCA <- function() {}


GetHeatmap <- function() {}