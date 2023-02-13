# make various count data
CountExtractor <- function(data_path) {
    if (require("tidyverse")) {
    
        files <- list.files(data_path)
        
        # make templates using a single data
        template <- read.table(paste0(data_path, files[1]), header = TRUE)
        sample_name <- stringr::str_split(files[1], '\\.')[[1]][1]
        
        LENGTH <- template[, c(1, 3)]  # gene length
        colnames(LENGTH)[2] <- sample_name
        
        COUNT <- template[, c(1, 5)]  # gene count
        colnames(COUNT)[2] <- sample_name
        
        TPM <- template[, c(1, 6)]  # TPM
        colnames(TPM)[2] <- sample_name
        
        FPKM <- template[, c(1, 7)]  # FPKM
        colnames(FPKM)[2] <- sample_name
        
        # add data to template
        for (i in 2:length(files)) {
            data <- read.table(paste0(data_path, files[i]), header = TRUE)
            sample_name <- stringr::str_split(files[i], '\\.')[[1]][1]
            
            length <- data[, c(1, 3)]
            colnames(length)[2] <- sample_name
            
            count <- data[, c(1, 5)]
            colnames(count)[2] <- sample_name
            
            tpm <- data[, c(1, 6)]
            colnames(tpm)[2] <- sample_name
            
            fpkm <- data[, c(1, 7)]
            colnames(fpkm)[2] <- sample_name
            
            LENGTH <- merge(LENGTH, length, by = "gene_id")
            COUNT <- merge(COUNT, count, by = "gene_id")
            TPM <- merge(TPM, tpm, by = "gene_id")
            FPKM <- merge(FPKM, fpkm, by = "gene_id")
        }
        
        LENGTH <- na.omit(LENGTH)
        LENGTH$gene_id <- unlist(lapply(LENGTH$gene_id,
                                        function(x) stringr::str_split(x, '\\.')[[1]][1]),
                                 use.names = FALSE)
        LENGTH <- LENGTH %>% tibble::column_to_rownames("gene_id")
        
        COUNT <- na.omit(COUNT)
        COUNT$gene_id <- unlist(lapply(COUNT$gene_id,
                                       function(x) stringr::str_split(x, '\\.')[[1]][1]),
                                use.names = FALSE)
        COUNT <- COUNT %>% tibble::column_to_rownames("gene_id")
        
        TPM <- na.omit(TPM)
        TPM$gene_id <- unlist(lapply(TPM$gene_id,
                                     function(x) stringr::str_split(x, '\\.')[[1]][1]),
                              use.names = FALSE)
        TPM <- TPM %>% tibble::column_to_rownames("gene_id")
        
        FPKM <- na.omit(FPKM)
        FPKM$gene_id <- unlist(lapply(FPKM$gene_id,
                                      function(x) stringr::str_split(x, '\\.')[[1]][1]),
                               use.names = FALSE)
        FPKM <- FPKM %>% tibble::column_to_rownames("gene_id")
        
        res <- list(LENGTH, COUNT, TPM, FPKM)
        names(res) <- c('LENGTH', 'COUNT', 'TPM', 'FPKM')
        
        return(res)
        
    } else {
        print("required library is not loaded!\n")
    }
}




# query different gene ID
BuildGeneNameDB <- function(gene_list) {
    if (required("biomaRt")) {
        
        
        
    } else {
        print("required library is not loaded!\n")
    }
}




# replace gene ID with different IDs seamlessly
# (i.e. ENS <--> Entrez, ENS <--> Gene_Symbol)
