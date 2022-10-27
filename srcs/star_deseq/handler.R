# handler.R

library(tidyverse)
library(stringr)


data_path <- "../../data/star/"


BuildCountData <- function() {
    count_file_list <- list.files(data_path, pattern = "ReadsPerGene")
    
    rep_idx <- 1
    
    # merge all STAR outputs
    repeat {
        count_file_path = file.path(data_path, count_file_list[rep_idx])
        
        tmp_data <- read.table(count_file_path, skip = 4)
        tmp_data <- tmp_data %>% dplyr::select(c('V1', 'V2'))
        
        col_name <- stringr::str_split(count_file_list[rep_idx],
                                       pattern = "_ReadsPerGene")[[1]][1]
        
        colnames(tmp_data) <- c("ENSG", col_name)
        
        if (exists("count_data")) {
            count_data <- merge(count_data, tmp_data, by = "ENSG", all = TRUE)
        } else {
            count_data <- tmp_data
        }
        
        if (rep_idx == length(count_file_list)) break
        
        rep_idx = rep_idx + 1
    }
    
    # remove gene version
    count_data$ENSG <- sapply(
        count_data$ENSG,
        function(x) stringr::str_split(x, pattern = "\\.")[[1]][1]
    )
    
    count_data <- count_data %>% tibble::column_to_rownames(var = "ENSG")
    
    count_data[] <- sapply(count_data, as.numeric)
    
    count_data <- count_data[apply(count_data, 
                                   MARGIN = 1, 
                                   function(x) all(x != 0)), ]
    
    return(count_data)
}


BuildMetaData <- function() {
    meta_file_path <- file.path(data_path, "meta_data")
    
    meta_data <- read.table(meta_file_path, 
                            header = TRUE,
                            stringsAsFactors = TRUE)
    
    # ensure order of rownames are equal to colnames of count data
    if (!exists("count_data")) count_data <<- BuildCountData()
    
    meta_data <- meta_data[
        match(colnames(count_data), meta_data$Samples), 
    ]

    # relevel factors (columns)
    meta_data$Treatment <- factor(meta_data$Treatment,
                                  levels = c("oil", "ccl4_1x", "ccl4_3x"))
    meta_data$Diet <- factor(meta_data$Diet,
                             levels = c("ncd", "nashd"))
    meta_data$Tumor <- factor(meta_data$Tumor,
                              levels = c("non_tumor", "adj_tumor", "tumor"))
    meta_data$Group1 <- factor(meta_data$Group1,
                               levels = c("NORMAL", "NASH", "CCL4", "HCC"))
    meta_data$Group2 <- factor(meta_data$Group2,
                               levels = c("NCD_CORN_OIL", "NASH_CORN_OIL",
                                          "NCD_LOW_CCL4", "NCD_HIGH_CCL4",
                                          "NASH_LOW_CCL4", "NASH_HIGH_CCL4"))
    
    # ensure absence of rownames before column_to_rownames
    rownames(meta_data) <- NULL
    meta_data <- meta_data %>% tibble::column_to_rownames(var = "Samples")

    return(meta_data)
}


presets <- list(
    c(26, 27, 28, 23, 24, 25),  # Normal vs NASH
    c(26, 27, 28, 1, 2, 3, 4, 5, 6),  # Normal vs CCL4 (1x, 3x ; weak HCC)
    c(26, 27, 28, 7, 8, 9, 10, 11, 12, 13,
      14, 15, 16, 17, 18, 19, 20, 21, 22),  # Normal vs HCC (1x, 3x, AT, T)
    c(23, 24, 25, 1, 2, 3, 4, 5, 6),  # NASH vs CCL4 (1x, 3x ; weak HCC)
    c(23, 24, 25, 7, 8, 9, 10, 11, 12, 13,
      14, 15, 16, 17, 18, 19, 20, 21, 22),  # NASH vs HCC (1x, 3x, AT, T)
    c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
      13, 14, 15, 16, 17, 18, 19, 20, 21, 22),  # CCL4 vs HCC,
    c(7, 8, 9, 10, 15, 16, 17, 18,
      11, 12, 13, 14, 19, 20, 21, 22)  # HCC_AT vs HCC_T
)

names(presets) <- c(
    "NORMAL_vs_NASH", "NORMAL_vs_CCL4", "NORMAL_vs_HCC",
    "NASH_vs_CCL4", "NASH_vs_HCC", "CCL4_vs_HCC",
    "HCC_AT_vs_HCC_T"
)


ModifyData <- function(subset) {
    if (!exists("count_data")) count_data <<- BuildCountData()
    if (!exists("meta_data")) meta_data <<- BuildMetaData()
    
    res <- list(count_data[, subset], meta_data[subset, ])
    names(res) <- c("count_data", "meta_data")
    
    return(res)
}