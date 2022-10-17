source("annotation.R")
annotation <- BuildAnnotation()


library(tximport)
library(DESeq2)
library(apeglm)
library(ashr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(writexl)


heat_colors <- rev(RColorBrewer::brewer.pal(11, 'RdBu'))
result_path <- "../../results/"

# consider subgroups 
PrepMatrix <- function(subgroup) {
    # remember output of read.table is a dataframe --> you need to specify column
    file_paths <- read.table('../data/kallisto/kallisto_file_list',
                             sep = '\t',
                             header = FALSE)$V1
    
    meta_data <- read.table('../data/kallisto/SampleMetaData',
                            sep = '\t',
                            header = TRUE,
                            stringsAsFactors = TRUE)
    
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
    
    # retrieve subgroup of interest
    file_paths <- file_paths[subgroup]
    meta_data <- meta_data[subgroup, ]
    
    txi <- tximport::tximport(file_paths,
                              type = "kallisto",
                              tx2gene = tx2gene,
                              ignoreTxVersion = TRUE)
    
    res <- list(meta_data, txi)
    names(res) <- c("mData", 'txi')
    
    return(res)
}


RunDeseq <- function(prepared_mat, design) {
    dds <- DESeq2::DESeqDataSetFromTximport(countData = prepared_mat$txi, 
                                            colData = prepared_mat$mData,
                                            design = design)  # custom design
    
    ddsX <- DESeq2::DESeq(dds)
    
    rld <- DESeq2::rlog(ddsX, blind = TRUE)
    rld_mat <- SummarizedExperiment::assay(rld)
    
    pca <- stats::prcomp(t(rld_mat), scale = FALSE)
    percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    sd_ratio <- sqrt(percentVar[2] / percentVar[1])
    
    pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta_data)
    
    rld_cor <- stats::cor(rld_mat)
    
    for (i in 1:length(plt_hl)) {
        ggplot(pca_df, aes(PC1, PC2)) +
            geom_point(aes(color = plt_hl[i])) + 
            xlab(paste0("PC1 : ", percentVar[1], "%")) + 
            ylab(paste0("PC2 : ", percentVar[2], "%")) + 
            coord_fixed(ratio = sd_ratio)
    }
    
    rld_cor <- stats::cor(rld_mat)
    
    
}


BasePlots <- function(ddsX) {

}