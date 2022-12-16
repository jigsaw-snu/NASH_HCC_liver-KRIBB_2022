library(readxl)
library(tidyverse)
library(DESeq2)
library(fgsea)

cnt <- readxl::read_xlsx("../../data/QC/cleaned_cnt.xlsx")  # outlier removed counts

# average values for duplicated genes
cnt <- cnt %>% 
    dplyr::group_by(external_gene_name) %>%
    dplyr::mutate_each(mean) %>%
    dplyr::distinct()

cnt <- cnt %>% tibble::column_to_rownames(var = "external_gene_name")

meta <- read.table("../../data/QC/meta_data.txt", header = TRUE)
meta <- meta %>% tibble::column_to_rownames(var = "Samples")


test <- DESeq2::estimateSizeFactorsForMatrix(cnt)