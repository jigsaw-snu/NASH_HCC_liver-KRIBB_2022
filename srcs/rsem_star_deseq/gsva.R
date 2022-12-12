# Differential Analysis for Pathways using GSVA
# https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html#62_Differential_expression_at_pathway_level


################################################################################
# Example Codes
################################################################################
library(GSVA)
library(GSVAdata)  # for test dataset

data(leukemia)  # Expression data  ; load leukemia expression dataset from GSVAdata
data(c2BroadSets)  # Pathway data (from MSigDB)  ; load c2BroadSets pathway dataset from GSVAdata

leukemia_es <- GSVA::gsva(leukemia_eset, c2BroadSets, min.sz = 10, max.sz = 500)


library(limma)

mod <- stats::model.matrix(~ factor(leukemia_es$subtype))  # design matrix ; shape : (sample x factor)
colnames(mod) <- c("ALL", "MLL.vs.ALL")

fit <- limma::lmFit(leukemia_es, mod)  # fit gene-wise linear models ; in here, fit "pathway-wise" linear models
# above expression is equal to limma::lmFit(exprs(leukemia_es), mod)
# above expression is also equal to limma::lmFit(leukemia_es@assayData$exprs, mod)
# exprs(leukemia_es) (same as leukemia_es@assayData$exprs) has actual (pathway x samples) ES values

fit.e <- limma::eBayes(fit)  # compute (empirical Bayes) moderated t-statistics / F-statistics / log-odds

res <- limma::decideTests(fit.e, p.value = 0.01)

tt <- limma::topTable(fit.e, coef = 2, n = Inf)

DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.01]
DEpwys_es <- exprs(leukemia_es[DEpwys, ])
################################################################################




library(stringr)
library(tidyverse)
library(readxl)
library(GSVA)
library(limma)


cnt_data <- readxl:::read_xlsx("test.xlsx")
cnt_data <- cnt_data %>% tibble::column_to_rownames("...1")

groups <- sapply(colnames(cnt_data),
                 function(x) stringr::str_split(x, '_')[[1]][1])

test_fit <- limma::lmFit(cnt_data, stats::model.matrix(~factor(groups)))
test_fit.e <- limma::eBayes(test_fit)

test_res <- limma::decideTests(test_fit.e, p.value = 0.01)
test_table <- limma::topTable(test_fit.e, coef = 2, n = Inf)


