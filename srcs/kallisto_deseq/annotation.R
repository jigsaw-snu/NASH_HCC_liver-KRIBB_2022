# Prepare Transcript-Gene Annotation Tables
library(org.Mm.eg.db)


BuildAnnotation <- function() {
    txdb <- org.Mm.eg.db
    
    # keytypes(txdb)  ## show list of available keytypes
    key <- AnnotationDbi::keys(txdb, keytype = "ENSEMBLTRANS")  # key : ENSMUST
    key2 <- AnnotationDbi::keys(txdb, keytype = "ENSEMBL")  # key : ENSMUSG
    
    # columns(txdb)  ## show list of retriveable data using 'select'
    tx2gene <- AnnotationDbi::select(txdb, 
                                     keys = key, 
                                     columns = "ENSEMBL", 
                                     keytype = "ENSEMBLTRANS")
    
    ens2sym <- AnnotationDbi::select(txdb,
                                     keys = key2,
                                     columns = "SYMBOL",
                                     keytype = "ENSEMBL")
    
    res = list(tx2gene, ens2sym)
    names(res) <- c("tx2gene", "ens2sym")
    
    return(res)
}
