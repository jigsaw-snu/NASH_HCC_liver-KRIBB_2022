source("Helper.R")

RAWDATA_PATH <- "../../resources/rsem_count/"

DATA <- CountExtractor(RAWDATA_PATH)

GENE_NAME_DB <- BuildGeneNameDB(rownames(DATA$LENGTH))

