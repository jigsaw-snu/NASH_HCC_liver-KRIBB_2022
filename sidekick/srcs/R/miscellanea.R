library(ggplot2)


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


# TPM calculation
