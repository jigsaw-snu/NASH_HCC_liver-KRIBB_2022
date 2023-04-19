# workflow.R contains plain analysis codes without functional encapsulations
source("miscellanea.R")


RSRC <- "../../resources"
RAW_DATA <- file.path(RSRC, "rsem_count")
RDATA <- file.path(RSRC, "rdata")
FIGURES <- file.path(RSRC, "figures")




# 1. create count_data and meta_data
library(tidyverse)


# prepare count data
files <- list.files(RAW_DATA)

raw_data <- read.table(file.path(RAW_DATA, files[1]), header = TRUE, sep = '\t')

raw_count <- raw_data[, c("gene_id", "expected_count")]
colnames(raw_count)[ncol(raw_count)] <- strsplit(files[1], '\\.')[[1]][1]

gene_length <- raw_data[, c("gene_id", "length")]
colnames(gene_length)[ncol(gene_length)] <- strsplit(files[1], '\\.')[[1]][1]

tpm <- raw_data[, c("gene_id", "TPM")]
colnames(tpm)[ncol(tpm)] <- strsplit(files[1], '\\.')[[1]][1]

fpkm <- raw_data[, c("gene_id", "FPKM")]
colnames(fpkm)[ncol(fpkm)] <- strsplit(files[1], '\\.')[[1]][1]

for (idx in 2:length(files)) {
    raw_data <- read.table(file.path(RAW_DATA, files[idx]), header = TRUE, sep = '\t')
    
    tmp_count <- raw_data[, c("gene_id", "expected_count")]
    tmp_length <- raw_data[, c("gene_id", "length")]
    tmp_tpm <- raw_data[, c("gene_id", "TPM")]
    tmp_fpkm <- raw_data[, c("gene_id", "FPKM")]
    
    raw_count <- merge(raw_count, tmp_count, by = "gene_id")
    gene_length <- merge(gene_length, tmp_length, by = "gene_id")
    tpm <- merge(tpm, tmp_tpm, by = "gene_id")
    fpkm <- merge(fpkm, tmp_fpkm, by = "gene_id")
    
    colnames(raw_count)[ncol(raw_count)] <- strsplit(files[idx], '\\.')[[1]][1]
    colnames(gene_length)[ncol(gene_length)] <- strsplit(files[idx], '\\.')[[1]][1]
    colnames(tpm)[ncol(tpm)] <- strsplit(files[idx], '\\.')[[1]][1]
    colnames(fpkm)[ncol(fpkm)] <- strsplit(files[idx], '\\.')[[1]][1]
}

# remove gene id version
raw_count$gene_id <- sapply(raw_count$gene_id, function(x) strsplit(x, '\\.')[[1]][1])
gene_length$gene_id <- sapply(gene_length$gene_id, function(x) strsplit(x, '\\.')[[1]][1])
tpm$gene_id <- sapply(tpm$gene_id, function(x) strsplit(x, '\\.')[[1]][1])
fpkm$gene_id <- sapply(fpkm$gene_id, function(x) strsplit(x, '\\.')[[1]][1])

# make gene id as row names
raw_count <- raw_count %>% tibble::column_to_rownames(var = "gene_id")
gene_length <- gene_length %>% tibble::column_to_rownames(var = "gene_id")
tpm <- tpm %>% tibble::column_to_rownames(var = "gene_id")
fpkm <- fpkm %>% tibble::column_to_rownames(var = "gene_id")

# order column names for the ease of making meta data
raw_count <- raw_count[, order(colnames(raw_count))]
gene_length <- gene_length[, order(colnames(gene_length))]
tpm <- tpm[, order(colnames(tpm))]
fpkm <- fpkm[, order(colnames(fpkm))]

# make float to int (round down)
raw_count[] <- lapply(raw_count, function(x) round(x))
gene_length[] <- lapply(gene_length, function(x) round(x))
tpm[] <- lapply(tpm, function(x) round(x))
fpkm[] <- lapply(fpkm, function(x) round(x))

count_data <- raw_count

# use biomaRt to convert ensembl gene id to gene symbol
ensembl <- biomaRt::useMart(biomart = "ensembl",
                            dataset = "mmusculus_gene_ensembl")

bm_res <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
                         filters = "ensembl_gene_id",
                         mart = ensembl,
                         values = rownames(count_data))




# prepare meta data
# HCC : NASH diet + CCl4 treatment
treatment <- grepl("CCI4", colnames(raw_count)) | grepl("HCC", colnames(raw_count))  # logical vector for ccl4
dose <- grepl("3x", colnames(raw_count))  # logical vector for 3x dose
diet <- grepl("NASH", colnames(raw_count)) | grepl("HCC", colnames(raw_count))  # logical vector for NASH diet
tumor <- !grepl("AT", colnames(raw_count)) & grepl("T", colnames(raw_count))  # logical vector for tumor existence

map_treatment <- c("TRUE" = "ccl4", "FALSE" = "saline")
map_dose <- c("TRUE" = "3x", "FALSE" = "1x")
map_diet <- c("TRUE" = "nash_diet", "FALSE" = "normal_chow")
map_tumor <- c("TRUE" = "tumor", "FALSE" = "non-tumor")

stage <- c()

for (sample_name in colnames(count_data)) {
    if ("CCI4" %in% strsplit(sample_name, '_')[[1]]) {
        stage <- c(stage, "Fibrosis")
    }
    
    if ("NCD" %in% strsplit(sample_name, '_')[[1]]) {
        stage <- c(stage, "Normal")
    }
    
    if ("NASH" %in% strsplit(sample_name, '_')[[1]]) {
        stage <- c(stage, "NASH")
    }
    
    if ("HCC" %in% strsplit(sample_name, '_')[[1]] & "AT" %in% strsplit(sample_name, '_')[[1]]) {
        stage <- c(stage, "HCC_AT")
    }
    
    if ("HCC" %in% strsplit(sample_name, '_')[[1]] & !("AT" %in% strsplit(sample_name, '_')[[1]])) {
        stage <- c(stage, "HCC_T")
    }
}

meta_data <- data.frame(Treatment = map_treatment[as.character(treatment)],
                        Dose = map_dose[as.character(dose)],
                        Diet = map_diet[as.character(diet)],
                        Tumor = map_tumor[as.character(tumor)],
                        Stage = stage)

rownames(meta_data) <- colnames(raw_count)

meta_data$Treatment <- factor(meta_data$Treatment, levels = c("saline", "ccl4"))
meta_data$Dose <- factor(meta_data$Dose, levels = c("1x", "3x"))
meta_data$Diet <- factor(meta_data$Diet, levels = c("normal_chow", "nash_diet"))
meta_data$Tumor <- factor(meta_data$Tumor, levels = c("non-tumor", "tumor"))
meta_data$Stage <- factor(meta_data$Stage, levels = c("Normal", "NASH", "Fibrosis", "HCC_AT", "HCC_T"))




# 2. make DESeq object
library(DESeq2)
library(ashr)
library(ggrepel)
library(ggplot2)
library(scales)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
library(fgsea)
library(msigdbr)

hallmarks_gmt <- msigdbr::msigdbr(species = "Mus musculus", category = "H")
hallmarks <- split(x = hallmarks_gmt$gene_symbol, 
                   f = hallmarks_gmt$gs_name)


dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_data,
                                      colData = meta_data,
                                      design = ~Stage)

# Pre-filtering low counts
keep <- rowSums(BiocGenerics::counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds <- DESeq2::DESeq(dds)

# RLE plot for normalization status
vsd <- DESeq2::vst(dds, blind=TRUE)
showRLEplot(SummarizedExperiment::assay(vsd))

# show PCA of total groups
showPCAplot(vsd, meta_data, c("Stage", "Dose"))

# show HCC_T and HCC_AT groups
vsd_sub <- vsd[, vsd$Stage %in% c("HCC_T", "HCC_AT")]
showPCAplot(vsd_sub, meta_data, c("Stage", "Dose"), names = TRUE, text_size = 3)

exclude <- c("HCC_3x_T_5", "HCC_3x_AT_5")
new_count <- count_data[, !(names(count_data) %in% exclude), drop = FALSE]
new_meta <- meta_data[!(rownames(meta_data) %in% exclude), , drop = FALSE]

dds <- DESeq2::DESeqDataSetFromMatrix(countData = new_count,
                                      colData = new_meta,
                                      design = ~Stage)



corr <- stats::cor(SummarizedExperiment::assay(vsd))
ComplexHeatmap::Heatmap(mat = corr,
                        #top_annotation = top_annotation,
                        #left_annotation = left_annotation,
                        rect_gp = grid::gpar(col = "white", lwd = 0.5),
                        heatmap_legend_param = list(title = "corr"))


corr_sub <- stats::cor(SummarizedExperiment::assay(vsd_sub))
ComplexHeatmap::Heatmap(mat = corr_sub,
                        rect_gp = grid::gpar(col = "white", lwd = 0.5),
                        heatmap_legend_param = list(title = "corr"))


################################################################################
################################################################################

keep <- meta_data %>% dplyr::filter(Dose == "1x") %>% rownames()
new_count <- count_data[, keep, drop=FALSE]
new_meta <- meta_data[keep, , drop=FALSE]

dds <- DESeq2::DESeqDataSetFromMatrix(countData = new_count,
                                      colData = new_meta,
                                      design = ~Stage)

dds <- DESeq2::DESeq(dds)

res <- DESeq2::lfcShrink(dds,
                         contrast = c("Stage", "HCC_T", "Normal"),
                         type = "ashr") %>% as.data.frame()

res <- merge(res %>% tibble::rownames_to_column("ensembl_gene_id"), bm_res)

resLFC_HCC_T_vs_Normal <- res %>%
    dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)


EnhancedVolcano::EnhancedVolcano(toptable = res,
                                 lab = res$external_gene_name,
                                 labSize = 0,
                                 x = "log2FoldChange",
                                 y = "padj",
                                 pCutoff = 0.05,
                                 FCcutoff = 2)

EnhancedVolcano::EnhancedVolcano(toptable = res,
                                 lab = res$external_gene_name,
                                 labSize = 3,
                                 labFace = "bold",
                                 drawConnectors = TRUE,
                                 x = "log2FoldChange",
                                 y = "padj",
                                 pCutoff = 10e-5,
                                 FCcutoff = 5)


# use biomaRt to convert ensembl gene id to gene symbol
bm_res_tmp <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
                             filters = "ensembl_gene_id",
                             mart = ensembl,
                             values = rownames(resLFC_HCC_T_vs_Normal))


resLFC_HCC_T_vs_Normal <- merge(resLFC_HCC_T_vs_Normal %>% tibble::rownames_to_column("ensembl_gene_id"), 
                                bm_res_tmp)

resLFC_HCC_T_vs_Normal2 <- resLFC_HCC_T_vs_Normal %>%
    dplyr::filter(abs(log2FoldChange) > 2)

up_reg_HCC_T_vs_Normal <- resLFC_HCC_T_vs_Normal2 %>%
    dplyr::filter(log2FoldChange > 0)

down_reg_HCC_T_vs_Normal <- resLFC_HCC_T_vs_Normal2 %>%
    dplyr::filter(log2FoldChange < 0)


up_kegg_HCC_T_vs_Normal <- clusterProfiler::enrichKEGG(gene = up_reg_HCC_T_vs_Normal$entrezgene_id,
                                                       organism = "mmu",
                                                       keyType = "kegg",
                                                       pvalueCutoff = 0.05,
                                                       qvalueCutoff = 0.05)

up_kegg_HCC_T_vs_Normal_sig <- up_kegg_HCC_T_vs_Normal@result %>%
    dplyr::filter(qvalue < 0.05)


ggplot2::ggplot(data = up_kegg_HCC_T_vs_Normal_sig, 
                aes(x = Count, y = Description))+
    ggplot2::theme_bw() +
    ggplot2::geom_bar(stat = "identity", 
                      aes(fill = -log(qvalue))) +
    ggplot2::scale_fill_gradient(low = "#1E88E5", high = "#E53935")



up_go_HCC_T_vs_Normal <- clusterProfiler::enrichGO(gene = up_reg_HCC_T_vs_Normal$entrezgene_id,
                                                   OrgDb = org.Mm.eg.db,
                                                   keyType = "ENTREZID",
                                                   ont = "BP",
                                                   pvalueCutoff = 0.05,
                                                   qvalueCutoff = 0.05,
                                                   readable = TRUE)

up_go_HCC_T_vs_Normal_sig <- up_go_HCC_T_vs_Normal@result %>%
    dplyr::filter(qvalue < 0.05)


ggplot2::ggplot(data = up_go_HCC_T_vs_Normal_sig %>% filter(qvalue < 3*10e-9), 
                aes(x = Count, y = Description))+
    ggplot2::theme_bw() +
    ggplot2::geom_bar(stat = "identity", 
                      aes(fill = -log(qvalue))) +
    ggplot2::scale_fill_gradient(low = "#1E88E5", high = "#E53935")



rank_HCC_T_vs_Normal <- resLFC_HCC_T_vs_Normal %>%
    dplyr::select(external_gene_name, log2FoldChange) %>%
    dplyr::arrange(desc(log2FoldChange)) %>%
    tibble::deframe()

gsea_HCC_T_vs_Normal <- fgsea::fgsea(pathways = hallmarks,
                                     stats = rank_HCC_T_vs_Normal,
                                     scoreType = "std") %>%
    dplyr::arrange(padj)

fgsea::plotEnrichment(pathway = hallmarks[["HALLMARK_G2M_CHECKPOINT"]],
                      stats = rank_HCC_T_vs_Normal) + 
    labs(title = "HALLMARK_G2M_CHECKPOINT")



################################################################################
################################################################################

keep <- meta_data %>% dplyr::filter(Dose == "3x" | Treatment == "saline") %>% rownames()
new_count <- count_data[, keep, drop=FALSE]
new_meta <- meta_data[keep, , drop=FALSE]

dds <- DESeq2::DESeqDataSetFromMatrix(countData = new_count,
                                      colData = new_meta,
                                      design = ~Stage)

dds <- DESeq2::DESeq(dds)


res <- DESeq2::lfcShrink(dds,
                         contrast = c("Stage", "HCC_T", "Normal"),
                         type = "ashr") %>% as.data.frame()

res <- merge(res %>% tibble::rownames_to_column("ensembl_gene_id"), bm_res)


EnhancedVolcano::EnhancedVolcano(toptable = res,
                                 lab = res$external_gene_name,
                                 labSize = 3,
                                 labFace = "bold",
                                 drawConnectors = TRUE,
                                 x = "log2FoldChange",
                                 y = "padj",
                                 pCutoff = 10e-5,
                                 FCcutoff = 5)


resLFC_HCC_T_vs_Normal_3x <- res %>%
    dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)

resLFC_HCC_T_vs_Normal2_3x <- resLFC_HCC_T_vs_Normal_3x %>%
    dplyr::filter(abs(log2FoldChange) > 2)

up_reg_HCC_T_vs_Normal_3x <- resLFC_HCC_T_vs_Normal2_3x %>%
    dplyr::filter(log2FoldChange > 0)

down_reg_HCC_T_vs_Normal_3x <- resLFC_HCC_T_vs_Normal2_3x %>%
    dplyr::filter(log2FoldChange < 0)


up_kegg_HCC_T_vs_Normal_3x <- clusterProfiler::enrichKEGG(gene = up_reg_HCC_T_vs_Normal_3x$entrezgene_id,
                                                          organism = "mmu",
                                                          keyType = "kegg",
                                                          pvalueCutoff = 0.05,
                                                          qvalueCutoff = 0.05)

up_kegg_HCC_T_vs_Normal_sig_3x <- up_kegg_HCC_T_vs_Normal_3x@result %>%
    dplyr::filter(qvalue < 0.05)


ggplot2::ggplot(data = up_kegg_HCC_T_vs_Normal_sig_3x, 
                aes(x = Count, y = Description))+
    ggplot2::theme_bw() +
    ggplot2::geom_bar(stat = "identity", 
                      aes(fill = -log(qvalue))) +
    ggplot2::scale_fill_gradient(low = "#1E88E5", high = "#E53935")


up_go_HCC_T_vs_Normal_3x <- clusterProfiler::enrichGO(gene = up_reg_HCC_T_vs_Normal_3x$entrezgene_id,
                                                      OrgDb = org.Mm.eg.db,
                                                      keyType = "ENTREZID",
                                                      ont = "BP",
                                                      pvalueCutoff = 0.05,
                                                      qvalueCutoff = 0.05,
                                                      readable = TRUE)

up_go_HCC_T_vs_Normal_3x <- up_go_HCC_T_vs_Normal_3x@result %>%
    dplyr::filter(qvalue < 0.05)

ggplot2::ggplot(data = up_go_HCC_T_vs_Normal_3x %>% filter(qvalue < 10e-9), 
                aes(x = Count, y = Description))+
    ggplot2::theme_bw() +
    ggplot2::geom_bar(stat = "identity", 
                      aes(fill = -log(qvalue))) +
    ggplot2::scale_fill_gradient(low = "#1E88E5", high = "#E53935")


rank_HCC_T_vs_Normal_3x <- resLFC_HCC_T_vs_Normal_3x %>%
    dplyr::select(external_gene_name, log2FoldChange) %>%
    dplyr::arrange(desc(log2FoldChange)) %>%
    tibble::deframe()

gsea_HCC_T_vs_Normal_3x <- fgsea::fgsea(pathways = hallmarks,
                                        stats = rank_HCC_T_vs_Normal_3x,
                                        scoreType = "std") %>%
    dplyr::arrange(padj)

fgsea::plotEnrichment(pathway = hallmarks[["HALLMARK_G2M_CHECKPOINT"]],
                      stats = rank_HCC_T_vs_Normal_3x) + 
    labs(title = "HALLMARK_G2M_CHECKPOINT")

fgsea::plotEnrichment(pathway = hallmarks[["HALLMARK_BILE_ACID_METABOLISM"]],
                      stats = rank_HCC_T_vs_Normal_3x) + 
    labs(title = "HALLMARK_BILE_ACID_METABOLISM")


