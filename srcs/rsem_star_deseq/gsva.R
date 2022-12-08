# Differential Analysis for Pathways using GSVA
# https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html#62_Differential_expression_at_pathway_level


library(GSVA)

################################################################################
# Example Codes
################################################################################
library(GSVAdata)  # for test dataset

data(leukemia)  # Expression data  ; load leukemia expression dataset from GSVAdata
data(c2BroadSets)  # Pathway data (from MSigDB)  ; load c2BroadSets pathway dataset from GSVAdata

leukemia_es <- GSVA::gsva(leukemia_eset, c2BroadSets, min.sz=10, max.sz=500)


library(limma)

mod <- stats::model.matrix(~ factor(leukemia_es$subtype))  # design matrix ; shape : (sample x factor)
colnames(mod) <- c("ALL", "MLL.vs.ALL")

fit <- limma::lmFit(leukemia_es, mod)  # fit gene-wise linear models ; in here, fit "pathway-wise" linear models
fit <- limma::eBayes(fit)

res <- limma::decideTests(fit, p.value=0.01)

tt <- limma::topTable(fit, coef=2, n=Inf)
################################################################################