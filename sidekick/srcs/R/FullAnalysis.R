source("Helper.R")

DATA_PATH <- "../../resources"
RAWDATA_PATH <- file.path(DATA_PATH, "rsem_count")
FIGURE_PATH <- file.path(DATA_PATH, "figures")

DATA <- CountExtractor(RAWDATA_PATH, DATA_PATH)

GENE_NAME_DB <- BuildGeneNameDB(rownames(DATA$LENGTH), "../../resources")

COUNT_sym <- ReplaceRowNames(DATA$COUNT, GENE_NAME_DB, "ens2sym")
COUNT_ent <- ReplaceRowNames(DATA$COUNT, GENE_NAME_DB, "ens2ent")


