library(readxl)
library(biomaRt)
library(writexl)

up_kegg <- readxl::read_xlsx("../../resources/rdata/up_kegg_HCC_T_vs_Normal_sig.xlsx")
down_kegg <- readxl::read_xlsx("../../resources/rdata/down_kegg_HCC_T_vs_Normal_sig.xlsx")


# (i) pool all the entrez gene ids for query
up_entrez <- strsplit(up_kegg$geneID, split = '/')
down_entrez <- strsplit(down_kegg$geneID, split = '/')
entrez_genes <- c(unlist(up_entrez), unlist(down_entrez))

# (ii) query biomart and get gene symbols
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL")
dset <- biomaRt::useDataset(dataset = "mmusculus_gene_ensembl", mart = mart)
gene_map <- biomaRt::getBM(attributes = c("entrezgene_id", "external_gene_name"),
                           filters = "entrezgene_id",
                           values = entrez_genes,
                           mart = dset)

# (iii) replace entrez gene id with gene symbol
up_symbols <- lapply(up_entrez,
                     function(t) paste(gene_map[match(t, gene_map$entrezgene_id), 2], collapse = '/'))
down_symbols <- lapply(down_entrez,
                       function(t) paste(gene_map[match(t, gene_map$entrezgene_id), 2], collapse = '/'))


up_kegg$geneID <- as.character(up_symbols)
down_kegg$geneID <- as.character(down_symbols)

writexl::write_xlsx(up_kegg, path = "../../resources/rdata/up_kegg_HCC_T_vs_Normal_sig2.xlsx")
writexl::write_xlsx(down_kegg, path = "../../resources/rdata/down_kegg_HCC_T_vs_Normal_sig2.xlsx")
