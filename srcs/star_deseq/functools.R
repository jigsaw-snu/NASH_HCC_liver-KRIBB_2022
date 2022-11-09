# functools.R

library(gprofiler2)
library(clusterProfiler)
library(org.Mm.eg.db)


ora <- function(deseq_result, fc) {
    if (fc <= 0) {
        print("fc parameter should be a positive value.")
        stop()
    }
    
    # make sure no duplicates are in
    up_list_go <- deseq_result$sig_lfc %>%
        dplyr::filter(!duplicated(ENSG)) %>%
        dplyr::filter(log2FoldChange > fc) %>%
        dplyr::pull(ENSG)
    
    down_list_go <- deseq_result$sig_lfc %>%
        dplyr::filter(!duplicated(ENSG)) %>%
        dplyr::filter(log2FoldChange < -fc) %>%
        dplyr::pull(ENSG)
    
    up_list_kegg <- deseq_result$sig_lfc_with_ent %>%
        dplyr::filter(!duplicated(ENTREZ)) %>%
        dplyr::filter(log2FoldChange > fc) %>%
        dplyr::pull(ENTREZ)
    
    down_list_kegg <- deseq_result$sig_lfc_with_ent %>%
        dplyr::filter(!duplicated(ENTREZ)) %>%
        dplyr::filter(log2FoldChange < -fc) %>%
        dplyr::pull(ENTREZ)
    
    go_up <- clusterProfiler::enrichGO(
        gene = up_list_go,
        OrgDb = "org.Mm.eg.db",
        keyType = "ENSEMBL",
        ont = "BP",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 1000
    )
    
    go_down <- clusterProfiler::enrichGO(
        gene = down_list_go,
        OrgDb = "org.Mm.eg.db",
        keyType = "ENSEMBL",
        ont = "BP",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 1000
    )
    
    kegg_up <- clusterProfiler::enrichKEGG(
        gene = up_list_kegg,
        organism = "mmu",
        keyType = "ncbi-geneid",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 1000
    )
    
    kegg_down <- clusterProfiler::enrichKEGG(
        gene = down_list_kegg,
        organism = "mmu",
        keyType = "ncbi-geneid",
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


gsea <- function(deseq_result, fc) {
    if (fc <= 0) {
        print("fc parameter should be a positive value.")
        stop()
    }
    
    # make sure no duplicates are in 
    up_df_go <- deseq_result$sig_lfc %>%
        dplyr::filter(!duplicated(ENSG)) %>%
        dplyr::filter(log2FoldChange > fc) %>%
        dplyr::mutate(rank = rank(log2FoldChange, ties.method = "random"))
    
    up_list_go <- up_df_go %>% dplyr::pull(rank)
    
    names(up_list_go) <- up_df_go %>% dplyr::pull(ENSG)
    
    up_list_go <- sort(up_list_go, decreasing = TRUE)
    
    down_df_go <- deseq_result$sig_lfc %>%
        dplyr::filter(!duplicated(ENSG)) %>%
        dplyr::filter(log2FoldChange < -fc) %>%
        dplyr::mutate(rank = rank(log2FoldChange, ties.method = "random"))
    
    down_list_go <- down_df_go %>% dplyr::pull(rank)
    
    names(down_list_go) <- down_df_go %>% dplyr::pull(ENSG)
    
    down_list_go <- sort(down_list_go, decreasing = TRUE)
    
    up_df_kegg <- deseq_result$sig_lfc_with_ent %>%
        dplyr::filter(!duplicated(ENTREZ)) %>%
        dplyr::filter(log2FoldChange > fc) %>%
        dplyr::mutate(rank = rank(log2FoldChange, ties.method = "random"))
    
    up_list_kegg <- up_df_kegg %>% dplyr::pull(rank)
    
    names(up_list_kegg) <- up_df_kegg %>% dplyr::pull(ENTREZ)
    
    up_list_kegg <- sort(up_list_kegg, decreasing = TRUE)
    
    down_df_kegg <- deseq_result$sig_lfc_with_ent %>%
        dplyr::filter(!duplicated(ENTREZ)) %>%
        dplyr::filter(log2FoldChange < -fc) %>%
        dplyr::mutate(rank = rank(log2FoldChange, ties.method = "random"))
    
    down_list_kegg <- down_df_kegg %>% dplyr::pull(rank)
    
    names(down_list_kegg) <- down_df_kegg %>% dplyr::pull(ENTREZ)
    
    down_list_kegg <- sort(down_list_kegg, decreasing = TRUE)
    
    go_up <- clusterProfiler::gseGO(
        geneList = up_list_go,
        OrgDb = "org.Mm.eg.db",
        keyType = "ENSEMBL",
        ont = "BP",
        minGSSize = 10,
        maxGSSize = 1000,
        scoreType = "pos",
        pvalueCutoff = 0.05
    )
    
    go_down <- clusterProfiler::gseGO(
        geneList = down_list_go,
        OrgDb = "org.Mm.eg.db",
        keyType = "ENSEMBL",
        ont = "BP",
        minGSSize = 10,
        maxGSSize = 1000,
        scoreType = "pos",
        pvalueCutoff = 0.05
    )
    
    kegg_up <- clusterProfiler::gseKEGG(
        geneList = up_list_kegg,
        organism = "mmu",
        keyType = "ncbi-geneid",
        pvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 1000,
        scoreType = "pos"
    )
    
    kegg_down <- clusterProfiler::gseKEGG(
        geneList = down_list_kegg,
        organism = "mmu",
        keyType = "ncbi-geneid",
        pvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 1000,
        scoreType = "pos"
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
