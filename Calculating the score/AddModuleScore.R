# AddModuleScore
#' @param object seurat_obj: a Seurat object of malignant cell
#' @param features GeneSets: list of gene sets
#' @param name names(GeneSets)

library(Seurat)
score <- AddModuleScore(object = seurat_obj, features = GeneSets, name = names(GeneSets))


