# Box plot
showBoxplot <- function(x) {
    color <- c(rep("darkgreen", 6),
               rep("red", 16),
               rep("yellow", 3),
               rep("blue", 3))
    
    ggplot2::ggplot(data = tidyr::gather(data = x,
                                         key = Sample,
                                         value = Count),
                    aes(x = Sample, y = Count)) +
        ggplot2::geom_boxplot(fill = color, alpha = rep(0.3, 28)) + 
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = element_text(angle = 45, 
                                                  vjust = 1, 
                                                  hjust = 1))
}


# RLE plot
# 0. log transform the original matrix
# 1. For each gene, calculate its median expression across the samples
# 2. calculate deviations from this median
# 3. generate a boxplot of all the deviations for each sample
showRLEplot <- function(x) {  # x : gene by sample matrix
    y <- log(x+1)  # log(x + pseudo_count)
    med <- apply(y, 1, median)  # gene-wise median
    rle <- apply(y, 2, function(x) (x - med))
    
    color <- c(rep("darkgreen", 6),
               rep("red", 16),
               rep("yellow", 3),
               rep("blue", 3))
    
    # return data type of apply function is matrix
    # tidyr::gather needs dataframe
    ggplot2::ggplot(data = tidyr::gather(data = as.data.frame(rle), 
                                         key = Sample, 
                                         value = Count),
                    aes(x = Sample, y = Count)) + 
        ggplot2::geom_hline(yintercept = 0, color = "red", linewidth = 1) +
        ggplot2::geom_boxplot(fill = color, alpha = rep(0.3, 28)) +
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = element_text(angle = 45, 
                                                  vjust = 1, 
                                                  hjust = 1)) +
        xlab("Samples") + ylab("Relative Log Expression")
}


# PCA plot
showPCAplot <- function(trf_cnt, mdata, intgroup, ntop=500, names = FALSE, text_size = NULL) {
    # trf_cnt must be result of vst or rlog function with blind = TRUE
    rv <- matrixStats::rowVars(SummarizedExperiment::assay(trf_cnt))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- stats::prcomp(t(SummarizedExperiment::assay(trf_cnt)[select, ]), center = TRUE, scale. = TRUE)
    
    percent_var <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, digits = 2)
    pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3])
    pca_df <- merge(pca_df, mdata[, intgroup, drop = FALSE], by = "row.names") %>%
        tibble::column_to_rownames(var = "Row.names")
    
    base_color_levels <- levels(pca_df[[intgroup[1]]])
    color_compatibility <- unique(as.vector(pca_df[[intgroup[1]]]))
    color_labels <- base_color_levels[base_color_levels %in% color_compatibility]
    color_values <- scales::hue_pal()(length(color_labels))
    
    if (names) {
        if (is.null(text_size)) text_size <- 4  # default text size 4
        
        if (length(intgroup) != 1) {
            base_shape_levels <- levels(pca_df[[intgroup[2]]])
            shape_compatibility <- unique(as.vector(pca_df[[intgroup[2]]]))
            shape_labels <- base_shape_levels[base_shape_levels %in% shape_compatibility] 
            shape_values <- scales::shape_pal()(length(shape_labels))
            
            ggplot2::ggplot(pca_df, aes(x = PC1, y = PC2)) + 
                ggplot2::geom_point(size = 4, aes(color = pca_df[[intgroup[1]]], shape = pca_df[[intgroup[2]]])) +
                ggplot2::scale_color_manual(name = intgroup[1],
                                            labels = color_labels,
                                            values = color_values) +
                ggplot2::scale_shape_manual(name = intgroup[2],
                                            labels = shape_labels,
                                            values = shape_values) +
                ggplot2::xlab(paste0("PC1 : ", percent_var[1], '%')) +
                ggplot2::ylab(paste0("PC2 : ", percent_var[2], '%')) +
                ggplot2::theme_bw() +
                ggplot2::coord_fixed() +
                ggrepel::geom_text_repel(aes(label = rownames(pca_df)), size = text_size)
        } else {
            ggplot2::ggplot(pca_df, aes(x = PC1, y = PC2)) + 
                ggplot2::geom_point(size = 4, aes(color = pca_df[[intgroup[1]]], shape = pca_df[[intgroup[2]]])) +
                ggplot2::scale_color_manual(name = intgroup[1],
                                            labels = color_labels,
                                            values = color_values) +
                ggplot2::xlab(paste0("PC1 : ", percent_var[1], '%')) +
                ggplot2::ylab(paste0("PC2 : ", percent_var[2], '%')) +
                ggplot2::theme_bw() +
                ggplot2::coord_fixed() +
                ggrepel::geom_text_repel(aes(label = rownames(pca_df)), size = text_size)
        }
    } else {
        if (length(intgroup) != 1) {
            base_shape_levels <- levels(pca_df[[intgroup[2]]])
            shape_compatibility <- unique(as.vector(pca_df[[intgroup[2]]]))
            shape_labels <- base_shape_levels[base_shape_levels %in% shape_compatibility] 
            shape_values <- scales::shape_pal()(length(shape_labels))
            
            ggplot2::ggplot(pca_df, aes(x = PC1, y = PC2)) + 
                ggplot2::geom_point(size = 4, aes(color = pca_df[[intgroup[1]]], shape = pca_df[[intgroup[2]]])) +
                ggplot2::scale_color_manual(name = intgroup[1],
                                            labels = color_labels,
                                            values = color_values) +
                ggplot2::scale_shape_manual(name = intgroup[2],
                                            labels = shape_labels,
                                            values = shape_values) +
                ggplot2::xlab(paste0("PC1 : ", percent_var[1], '%')) +
                ggplot2::ylab(paste0("PC2 : ", percent_var[2], '%')) +
                ggplot2::theme_bw() +
                ggplot2::coord_fixed()
        } else {
            ggplot2::ggplot(pca_df, aes(x = PC1, y = PC2)) + 
                ggplot2::geom_point(size = 4, aes(color = pca_df[[intgroup[1]]], shape = pca_df[[intgroup[2]]])) +
                ggplot2::scale_color_manual(name = intgroup[1],
                                            labels = color_labels,
                                            values = color_values) +
                ggplot2::xlab(paste0("PC1 : ", percent_var[1], '%')) +
                ggplot2::ylab(paste0("PC2 : ", percent_var[2], '%')) +
                ggplot2::theme_bw() +
                ggplot2::coord_fixed()
        }
    }
}


# sample by sample heatmap
showSampleHeatmap <- function(trf_cnt, mdata, annot_top = NA, annot_left = NA) {
    corr <- stats::cor(SummarizedExperiment::assay(trf_cnt))
    
    if (is.nan(annot_top) & is.nan(annot_left)) {
        ComplexHeatmap::Heatmap(mat = corr,
                                top_annotation = annot_top,
                                left_annotation = annot_left,
                                rect_gp = grid::gpar(col = "white", lwd = 1),
                                heatmap_legend_param = list(title = "corr"))
    } else if (is.nan(annot_top)) {
        ComplexHeatmap::Heatmap(mat = corr,
                                top_annotation = annot_top,
                                rect_gp = grid::gpar(col = "white", lwd = 1),
                                heatmap_legend_param = list(title = "corr"))
    } else if (is.nan(annot_left)) {
        ComplexHeatmap::Heatmap(mat = corr,
                                left_annotation = annot_left,
                                rect_gp = grid::gpar(col = "white", lwd = 1),
                                heatmap_legend_param = list(title = "corr"))
    } else {
        ComplexHeatmap::Heatmap(mat = corr,
                                rect_gp = grid::gpar(col = "white", lwd = 1),
                                heatmap_legend_param = list(title = "corr"))
    }
}




# show sample correlation heatmap
annot_top <- ComplexHeatmap::HeatmapAnnotation(tumor = meta_data$Tumor,
                                               dose = meta_data$Dose,
                                               stage = meta_data$Stage,
                                               col = list(stage = c("Normal" = "#33FFB8",
                                                                    "NASH" = "#FFC300",
                                                                    "Fibrosis" = "#5E33FF",
                                                                    "HCC_AT" = "#33C1FF",
                                                                    "HCC_T" = "#900C3F"),
                                                          dose = c("1x" = "#CADE6E",
                                                                   "3x" = "#EFA751"),
                                                          tumor = c("tumor" = "#FF5733",
                                                                    "non-tumor" = "#51D5EF")))


annot_left <- ComplexHeatmap::rowAnnotation(tumor = meta_data$Tumor,
                                            dose = meta_data$Dose,
                                            stage = meta_data$Stage,
                                            col = list(stage = c("Normal" = "#33FFB8",
                                                                 "NASH" = "#FFC300",
                                                                 "Fibrosis" = "#5E33FF",
                                                                 "HCC_AT" = "#33C1FF",
                                                                 "HCC_T" = "#900C3F"),
                                                       dose = c("1x" = "#CADE6E",
                                                                "3x" = "#EFA751"),
                                                       tumor = c("tumor" = "#FF5733",
                                                                 "non-tumor" = "#51D5EF")))