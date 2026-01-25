# Identifying CO/ME patterns of consensus cancer cell states (CCSs)

# 1. Trajectory analysis
#' @param seurat_obj a Seurat object
#' @param order_gene a vector of features consisting of significantly upregulated genes identified in the differential analysis of CCSs

run_monocle <- function(seurat_obj, order_gene){
  
  library(monocle)
  library(Seurat)
  
  # Create a CellDateSet object
  expr_matrix <- as(as.matrix(seurat_obj@assays[["RNA"]]$counts), 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data = seurat_obj@meta.data)
  f_data <- data.frame(gene_short_name = row.names(seurat_obj),row.names = row.names(seurat_obj))
  fd <- new('AnnotatedDataFrame', data = f_data)
  cds <- newCellDataSet(expr_matrix, 
                        phenoData = pd, 
                        featureData = fd, 
                        lowerDetectionLimit = 0.5, 
                        expressionFamily = negbinomial.size())
  
  # Estimate size factor and dispersion
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  # Construct the trajectory
  cds <- setOrderingFilter(cds,order_gene)
  cds <- reduceDimension(cds,max_components = 2,method = 'DDRTree')
  cds <- orderCells(cds)
  
  return(cds)
}


# 2. Discretize pseudotime into equal-width bins
# seurat_obj: a Seurat object containing the "Monocle2_pse" column in its meta.data

bins <- round(2 * nrow(seurat_obj)^(1/3))
seurat_obj$Bin <- cut(seurat_obj$Monocle2_pse, 
                      breaks = bins, 
                      labels = paste0("bin", 1:bins), 
                      include.lowest = TRUE)


# 3. Binarize proportions using a Gaussian mixture model
#' @param obj a vector of CCS proportions across all bins
#' @param label a single CCS label

run_GMM <- function(obj, label){
  library(mclust)
  GMM = Mclust(obj, G = 2)
  
  # Plot
  pdf(paste0(label,"_GMM.pdf"),width = 5,height = 4)
  plot(GMM,"BIC")
  plot(GMM,"classification")
  plot(GMM,"density")
  plot(GMM,"uncertainty")
  dev.off()
  
  # Combine variables into a data frame
  result <- data.frame(classification = GMM$classification, score = obj)
  
  return(result)
}


# 4. Identify CO/ME patterns
# binary_matrix: a binary matrix with bins as rows and CCSs as columns

library(select)
sample_class <- rep("Bin", nrow(binary_matrix))
names(sample_class) <- rownames(binary_matrix)
alteration_class <- rep("Cell state", ncol(binary_matrix))
names(alteration_class) <- colnames(binary_matrix)
select_result <- select(M = binary_matrix, 
                        sample.class = sample_class,
                        alteration.class = alteration_class,
                        folder = 'select_result/', 
                        r.seed = 123, 
                        n.cores = 1,
                        min.feature.support = 1,
                        calculate_APC_threshold = T)

