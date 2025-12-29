# ssGSEA analysis
#' @param exprData Gene_exprData: the gene expression profiles
#' @param geneSets GeneSets: list of gene sets

library(GSVA)
param <- ssgseaParam(exprData = Gene_exprData, geneSets = GeneSets)
scores <- gsva(param, verbose = TRUE)

# z-score scale
zs_scores = t(apply(scores,1,scale))

# Group patients with top 3 scores
sample_assign = apply(zs_scores,2,function(x){
  rownames(zs_scores)[order(x,decreasing = T)[1:3]]
  })

