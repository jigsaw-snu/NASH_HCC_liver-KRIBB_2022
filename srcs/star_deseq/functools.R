# functools.R

library(gprofiler2)
library(clusterProfiler)
library(org.Mm.eg.db)


ora <- function(sig_lfc, fc) {
    if (fc <= 0) {
        print("fc parameter should be a positive value.")
        stop()
    }
    
    if (!("ENSG" %in% colnames(sig_lfc))) {
        print("ENSEMBL gene id is needed.")
        stop()
    }
    
    up_list <- sig_lfc %>%
        dplyr::filter(log2FoldChange > fc) %>%
        dplyr::pull(ENSG)
    
    down_list <- sig_lfc %>%
        dplyr::filter(log2FoldChange < -fc) %>%
        dplyr::pull(ENSG)
    
    go_up <- clusterProfiler::enrichGO(
        gene = up_list,
        OrgDb = "org.Mm.eg.db",
        keyType = "ENSEMBL",
        ont = "BP",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 1000
    )
    
    go_down <- clusterProfiler::enrichGO(
        gene = down_list,
        OrgDb = "org.Mm.eg.db",
        keyType = "ENSEMBL",
        ont = "BP",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 1000
    )
    
    kegg_up <- clusterProfiler::enrichKEGG(
        gene = up_list,
        organism = "mmu",
        keyType = "kegg",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 1000
    )
    
    kegg_down <- clusterProfiler::enrichKEGG(
        gene = down_list,
        organism = "mmu",
        keyType = "kegg",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 1000
    )
    
    res <- list(
        go_up@result %>% dplyr::filter(p.adjust < 0.05 & qvalue < 0.05),
        go_down@result %>% dplyr::filter(p.adjust < 0.05 & qvalue < 0.05),
        kegg_up@result %>% dplyr::filter(p.adjust < 0.05 & qvalue < 0.05),
        kegg_down@result %>% dplyr::filter(p.adjust < 0.05 & qvalue < 0.05)
    )
    names(res) <- c("go_up", "go_down", "kegg_up", "kegg_down")
    
    return(res)
}


gsea <- function(sig_lfc, fc) {
    if (fc <= 0) {
        print("fc parameter should be a positive value.")
        stop()
    }
    
    if (!("ENSG" %in% colnames(sig_lfc))) {
        print("ENSEMBL gene id is needed.")
        stop()
    }
    
    up_list <- sig_lfc %>%
        dplyr::filter(log2FoldChange > fc) %>%
        dplyr::pull(log2FoldChange)
    
    names(up_list) <- sig_lfc %>%
        dplyr::filter(log2FoldChange > fc) %>%
        dplyr::pull(ENSG)
    
    up_list <- sort(up_list, decreasing = TRUE)
    
    down_list <- sig_lfc %>%
        dplyr::filter(log2FoldChange < -fc) %>%
        dplyr::pull(log2FoldChange)
    
    names(down_list) <- sig_lfc %>%
        dplyr::filter(log2FoldChange < -fc) %>%
        dplyr::pull(ENSG)
    
    down_list <- sort(down_list, decreasing = TRUE)
    
    go_up <- clusterProfiler::gseGO(
        geneList = up_list,
        OrgDb = "org.Mm.eg.db",
        keyType = "ENSEMBL",
        ont = "BP",
        minGSSize = 10,
        maxGSSize = 1000,
        scoreType = "pos"
    )
    
    go_down <- clusterProfiler::gseGO(
        geneList = down_list,
        OrgDb = "org.Mm.eg.db",
        keyType = "ENSEMBL",
        ont = "BP",
        minGSSize = 10,
        maxGSSize = 1000,
        scoreType = "pos"
    )
    
    kegg_up <- clusterProfiler::gseKEGG(
        geneList = up_list,
        organism = "mmu",
        keyType = "kegg",
        pvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 1000
    )
    
    kegg_down <- clusterProfiler::gseKEGG(
        geneList = down_list,
        organism = "mmu",
        keyType = "kegg",
        pvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 1000
    )
    
    res <-
}