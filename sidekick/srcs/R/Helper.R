# make various count data
CountExtractor <- function(data_path, save_path) {
    # if saved file exists, just load it
    if ('CountData.Rds' %in% list.files(save_path)) {
        print("loading from Rds file...\n")
        
        return(readRDS(file.path(save_path, 'CountData.Rds')))
    }
    
    if (require("tidyverse")) {
    
        files <- list.files(data_path)
        
        # make templates using a single data
        template <- read.table(file.path(data_path, files[1]), header = TRUE)
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
            data <- read.table(file.path(data_path, files[i]), header = TRUE)
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
        
        # remove gene version and make it as rownames
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
        
        # remove zero-valued rows (should performed after 'column_to_rownames')
        # count data should not have zero-valued rows
        COUNT <- COUNT[apply(COUNT, MARGIN = 1, function(x) all(x != 0)), ]
        
        # LENGTH, TPM, FPKM data should have same rows as COUNT has
        LENGTH <- LENGTH[rownames(LENGTH) %in% rownames(COUNT), ]
        TPM <- TPM[rownames(TPM) %in% rownames(COUNT), ]
        FPKM <- FPKM[rownames(FPKM) %in% rownames(COUNT), ]
        
        res <- list(LENGTH, COUNT, TPM, FPKM)
        names(res) <- c('LENGTH', 'COUNT', 'TPM', 'FPKM')
        
        # save for next time
        saveRDS(res, file = file.path(save_path, 'CountData.Rds'))
        
        return(res)
        
    } else {
        stop("tidyverse library is not loaded!\n")
    }
}




# query different gene ID
BuildGeneNameDB <- function(gene_list, save_path) {
    # if saved file exists, jsut load it
    if ('GeneNameDB.Rds' %in% list.files(save_path)) {
        print("loading from Rds file...\n")
        
        return(readRDS(file.path(save_path, 'GeneNameDB.Rds')))
    }
    
    if (require("biomaRt")) {
        
        mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                 dataset = "mmusculus_gene_ensembl")
        
        # attributes : output features
        # filters : input features
        # values : input gene list
        # bm : dataframe returned by getBM function
        bm <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
                             filters = "ensembl_gene_id",
                             values = gene_list,
                             mart = mart)
        
        # Ensembl gene ID is referred as 'gene_id' in RSEM result
        colnames(bm) <- c("gene_id", "entrez", "symbol")
        
        # save for next time
        saveRDS(bm, file = file.path(save_path, 'GeneNameDB.Rds'))
        
        return(bm)
        
    } else {
        stop("biomaRt library is not loaded!\n")
    }
}




# replace gene ID with different IDs seamlessly
# (i.e. ENS <--> Entrez, ENS <--> Gene_Symbol)
ReplaceRowNames <- function(target_df, name_db, option) {
    # substitute_df should have columns for bot originanl row names and new row names
    
    if (require("tidyverse")) {
        
        if (option == 'ens2ent') {  # ensembl gene ID to entrez gene ID
            substitute_df <- name_db[, c(1, 2)]
            new_rowname <- "entrez"
        } else if (option == 'ens2sym') {  # ensembl gene ID to gene symbol
            substitute_df <- name_db[, c(1, 3)]
            new_rowname <- "symbol"
        } else {
            stop("Choose 'ens2ent' or 'ens2sym' as an option")
        }
        
        target_df <- target_df %>% tibble::rownames_to_column(var = "gene_id")
        target_df <- merge(target_df, substitute_df, by = "gene_id")
        target_df <- target_df[, !('gene_id' == colnames(target_df))]  # gene_id is no more needed
        
        # remove identically duplicated rows
        target_df <- target_df %>% dplyr::distinct()

        # sum multimapped rows (different Ensembl ID, same entrez/symbol)
        target_df <- target_df %>% 
            dplyr::group_by(!!rlang::sym(new_rowname)) %>%
            dplyr::mutate(dplyr::across(dplyr::everything(), sum)) %>%
            dplyr::distinct()
        
        test1 <<- target_df
        # remove unmapped rows (empty entrez/symbol) and unexpected empty/NA cells
        target_df <- subset(target_df, !is.na(new_rowname) & trimws(new_rowname) != '')
        test2 <<- target_df
        target_df <- target_df %>% tibble::column_to_rownames(var = new_rowname)
        
        return(target_df)
        
    } else {
        stop("tidyverse library is not loaded!")
    }
}





###
needed_packages <- c('a', 'b', 'c')

if (all(sapply(needed_packages, require, character.only = TRUE))) {
    ...
} else {
    ...
}
###




mutate_each(funs(mean), -(1:5)) %>% distinct