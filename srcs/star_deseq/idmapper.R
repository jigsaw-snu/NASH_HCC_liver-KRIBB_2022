# idmapper.R

library(biomaRt)


ens2sym <- function(dataset) {
    mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                             dataset = "mmusculus_gene_ensembl")
    
    # Exception handling for gene id column
    # in most case, there will be a gene id without gene version
    # column name referring gene id must be "ENSG"
    
    # 1. No Gene id
    if (!("ENSG" %in% colnames(dataset))) {
        print("No ENSG in columns")
        # check if rownames are ensemble gene id
        if (!is.null(rownames(dataset)) &
            +stringr::str_detect(rownames(dataset)[1], "ENSMUSG")) {
            dataset <- dataset %>% tibble::rownames_to_column(var = "ENSG")
        } else {
            stop("There seems no gene id column to convert!",
                 "\nCheck if the name of your gene id column is 'ENSG'",
                 "\nor Check if your dataset really contains gene id.")
        }
    }
    
    # 2. Gene id without version
    if (!stringr::str_detect(dataset$ENSG[1], '\\.')) {
        input_filter <- "ensembl_gene_id"
    } else {
    # 2. Gene id with version
        input_filter <- "ensembl_gene_id_version"
    }
    
    query <- biomaRt::getBM(
        # need input_filter as output to use it as a key to merge dataframes
        attributes = c(input_filter, "mgi_symbol"), # attributes : output filter
        filters = input_filter,  # filters : input filter
        values = dataset$ENSG,  # values : input value (vector)
        mart = mart  # mart : biomaRt object instantiated above
    )
    
    # re-set colnames for merging query result with input dataset
    colnames(query) <- c("ENSG", "gene_symbol")
    
    # if no matching gene symbol, drop that data
    res <- merge(dataset, query, by = "ENSG") 
    res <- res[, -1]
    res <- res[, c(colnames(res)[ncol(res)], colnames(res)[2:ncol(res) - 1])]
    
    return(res)
}


ens2entrez <- function() {
    mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                             dataset = "mmusculus_gene_ensembl")
    
    # Exception handling for gene id column
    # in most case, there will be a gene id without gene version
    # column name referring gene id must be "ENSG"
    
    # 1. No Gene id
    if (!("ENSG" %in% colnames(dataset))) {
        print("No ENSG in columns")
        # check if rownames are ensemble gene id
        if (!is.null(rownames(dataset)) &
            +stringr::str_detect(rownames(dataset)[1], "ENSMUSG")) {
            dataset <- dataset %>% tibble::rownames_to_column(var = "ENSG")
        } else {
            stop("There seems no gene id column to convert!",
                 "\nCheck if the name of your gene id column is 'ENSG'",
                 "\nor Check if your dataset really contains gene id.")
        }
    }
    
    # 2. Gene id without version
    if (!stringr::str_detect(dataset$ENSG[1], '\\.')) {
        input_filter <- "ensembl_gene_id"
    } else {
        # 2. Gene id with version
        input_filter <- "ensembl_gene_id_version"
    }
    
    query <- biomaRt::getBM(
        # need input_filter as output to use it as a key to merge dataframes
        attributes = c(input_filter, "entrezgene_id"), # attributes : output filter
        filters = input_filter,  # filters : input filter
        values = dataset$ENSG,  # values : input value (vector)
        mart = mart  # mart : biomaRt object instantiated above
    )
    
    # re-set colnames for merging query result with input dataset
    colnames(query) <- c("ENSG", "ENTREZ")
    
    # if no matching gene symbol, drop that data
    res <- merge(dataset, query, by = "ENSG") 
    res <- res[, -1]
    res <- res[, c(colnames(res)[ncol(res)], colnames(res)[2:ncol(res) - 1])]
    
    return(res)
}
