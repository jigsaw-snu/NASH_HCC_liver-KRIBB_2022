# workflow.R contains plain analysis codes without functional encapsulations

RSRC <- "../../resources/"
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


# gene 
# remove low count 'genes' (count_value <= 15)
raw_count <- raw_count[apply(raw_count, MARGIN = 1, function(x) all(x > 10)), ]

# find and remove outlier 'genes' using MAD criteria


# feature selection (use highly variable genes)


# prepare meta data
# HCC : NASH diet + CCl4 treatment
treatment <- grepl("CCI4", colnames(raw_count)) | grepl("HCC", colnames(raw_count))  # logical vector for ccl4
dose <- grepl("3x", colnames(raw_count))  # logical vector for 3x dose
diet <- grepl("NASH", colnames(raw_count)) | grepl("HCC", colnames(raw_count))  # logical vector for NASH diet
tumor <- grepl("T", colnames(raw_count))  # logical vector for tumor existence

map_treatment <- c("TRUE" = "ccl4", "FALSE" = "saline")
map_dose <- c("TRUE" = "3x", "FALSE" = "1x")
map_diet <- c("TRUE" = "nash_diet", "FALSE" = "normal_chow")
map_tumor <- c("TRUE" = "tumor", "FALSE" = "non-tumor")

meta_data <- data.frame(Treatment = map_treatment[as.character(treatment)],
                        Dose = map_dose[as.character(dose)],
                        Diet = map_diet[as.character(diet)],
                        Tumor = map_tumor[as.character(tumor)])

rownames(meta_data) <- colnames(raw_count)

meta_data$Treatment <- factor(meta_data$Treatment, levels = c("saline", "ccl4"))
meta_data$Dose <- factor(meta_data$Dose, levels = c("1x", "3x"))
meta_data$Diet <- factor(meta_data$Diet, levels = c("normal_chow", "nash_diet"))
meta_data$Tumor <- factor(meta_data$Tumor, levels = c("non-tumor", "tumor"))




# 2. make DESeq object
library(DESeq2)

dds <- DESeq2::DESeqDataSetFromMatrix(countData =)




# 3. EDA