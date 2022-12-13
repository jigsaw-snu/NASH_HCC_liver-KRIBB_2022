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


PathDiff <- function(file_path) {
    cnt_data <- readxl:::read_xlsx(file_path)
    cnt_data <- cnt_data %>% tibble::column_to_rownames("...1")
    
    groups <- sapply(colnames(cnt_data),
                     function(x) stringr::str_split(x, '_')[[1]][1])
    
    test_model <- stats::model.matrix(~0 + factor(groups))
    colnames(test_model) <- groups[!duplicated(groups)]
    
    test_cont <- limma::makeContrasts(
        NCDvsALL = NCD - (CCI4 + HCC + NASH) / 3,
        CCI4vsALL = CCI4 - (NCD + NASH + HCC) / 3,
        NASHvsALL = NASH - (NCD + CCI4 + HCC) / 3,
        HCCvsALL = HCC - (NCD + CCI4 + NASH) / 3,
        levels = test_model
    )
    
    test_fit <- limma::lmFit(cnt_data, test_model)
    test_fit.c <- limma::contrasts.fit(test_fit, contrasts = test_cont)
    test_fit.e <- limma::eBayes(test_fit.c)
    
    test_res <- limma::decideTests(test_fit.e, p.value = 0.05)
    
    test_table1 <- limma::topTable(test_fit.e, coef = 1, n = Inf)  # NCD vs ALL
    test_table2 <- limma::topTable(test_fit.e, coef = 2, n = Inf)  # CCL4 vs ALL
    test_table3 <- limma::topTable(test_fit.e, coef = 3, n = Inf)  # NASH vs ALL
    test_table4 <- limma::topTable(test_fit.e, coef = 4, n = Inf)  # HCC vs ALL
    
    res <- list(test_fit, test_fit.c, test_fit.e, test_res,
                test_table1, test_table2, test_table3, test_table4)
    
    names(res) <- c("fit", "fit.c", "fit.e", "res",
                    "tbl1", "tbl2", "tbl3", "tbl4")
    
    return(res)
}

test1 <- PathDiff("./test.xlsx")
test2 <- PathDiff("./test2.xlsx")
test3 <- PathDiff("./test3.xlsx")


PathDiff2 <- function(file_path) {
    cnt_data <- readxl:::read_xlsx(file_path)
    cnt_data <- cnt_data %>% tibble::column_to_rownames("...1")
    
    groups <- sapply(colnames(cnt_data),
                     function(x) stringr::str_split(x, '_')[[1]][1])
    
    test_model <- stats::model.matrix(~0 + factor(groups))
    colnames(test_model) <- groups[!duplicated(groups)]
    
    test_cont <- limma::makeContrasts(
        NCDvsALL = NCD - (CCI4 + HCC + NASH) / 3,
        CCI4vsALL = CCI4 - (NCD + NASH + HCC) / 3,
        NASHvsALL = NASH - (NCD + CCI4 + HCC) / 3,
        HCCvsALL = HCC - (NCD + CCI4 + NASH) / 3,
        levels = test_model
    )
    
    test_fit <- limma::lmFit(cnt_data, test_model)
    test_fit.c <- limma::contrasts.fit(test_fit, contrasts = test_cont)
    test_fit.e <- limma::eBayes(test_fit.c)
    
    test_res <- limma::decideTests(test_fit.e, p.value = 0.05)
    
    test_table1 <- limma::topTable(test_fit.e, coef = 1, n = Inf)  # NCD vs ALL
    test_table2 <- limma::topTable(test_fit.e, coef = 2, n = Inf)  # CCL4 vs ALL
    test_table3 <- limma::topTable(test_fit.e, coef = 3, n = Inf)  # NASH vs ALL
    test_table4 <- limma::topTable(test_fit.e, coef = 4, n = Inf)  # HCC vs ALL
    
    res <- list(test_fit, test_fit.c, test_fit.e, test_res,
                test_table1, test_table2, test_table3, test_table4)
    
    names(res) <- c("fit", "fit.c", "fit.e", "res",
                    "tbl1", "tbl2", "tbl3", "tbl4")
    
    return(res)
}

test4 <- PathDiff(("./test4.xlsx"))
