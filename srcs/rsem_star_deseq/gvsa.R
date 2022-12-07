# Differential Analysis for Pathways using GSVA
# https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html#62_Differential_expression_at_pathway_level


library(GSVA)

################################################################################
# Example Codes
################################################################################
library(GSVAdata)  # for test dataset

GSVAdata::data(leukemia)  # Expression data
GSVAdata::data(c2BroadSets)  # Pathway data (from MSigDB)

leukemia_es <- GSVA::gsva(leukemia_eset, c2BroadSets, min.sz=10, max.sz=500)

################################################################################