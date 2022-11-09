library(ggtext)
library(cowplot)

test_df <- data.frame(x = rep(0, 10), y = sample(-5:15, size = 10))

for (i in 1:10) {
    test_df_ <- data.frame(x = colnames(test_df[i, ]),
                           y = as.numeric(test_df[i, ]))
    
    if (test_df_$y[2] - test_df_$y[1] > 0) test_col = 'red'
    else test_col = 'blue'
    
    gg <- ggplot(test_df_, aes(x = x, y = y, group = 1)) +
        geom_line(color = test_col, size = 1, linetype = 'twodash') +
        geom_point(size = 5, aes(shape = x), show.legend = FALSE) +
        labs(title = paste0("test_p", 
                            i, 
                            ' : ',
                            "<span style='color: ", 
                            test_col, 
                            ";'>", 
                            test_df_$y[2] - test_df_$y[1])) +
        theme(axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              plot.title = element_markdown(hjust = 0.5, size = 10))
    
    assign(paste0('test_p', i), gg)
}

plot_grid(test_p1, test_p2, test_p3, test_p4, test_p5,
          #NULL, NULL, NULL, NULL, NULL,
          test_p6, test_p7, test_p8, test_p9, test_p10,
          nrow = 2, ncol = 5)
# rel_heights = c(rep(2, 5), rep(0.1, 5), rep(2, 5)))




modify_deseq_result <- function() {
    
}