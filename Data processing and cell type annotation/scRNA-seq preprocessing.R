# scRNA-seq data preprocessing

# Load R packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(rlang)
library(scater)


# Empty drops 
#' @param count_mat raw count matrix

run_emptyDrops <- function(count_mat,lower = 100) {
  
  library(DropletUtils)
  print(min(colSums(count_mat)))
  
  set.seed(100)
  out <- emptyDrops(count_mat, lower = lower)
  retained <- count_mat[,which(out$FDR < 0.05)]
  
  return(retained)
}


# Low complexity libraries
# count_mat: raw count matrix

seurat_obj <- CreateSeuratObject(counts = count_mat, project = "HNSCC", min.features = 200, names.field = 1, names.delim = "_")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
max_features <- sort(seurat_obj$nFeature_RNA, decreasing = T)[round(ncol(seurat_obj)*0.01)]
seurat_obj <- subset(seurat_obj, subset = percent.mt < 20 & nFeature_RNA < as.numeric(max_features))


# Doublet detection and removal
#' @param seurat_obj a Seurat object after filtering out low-quality cells

run_doubletFinder <- function(seurat_obj){
  library(DoubletFinder)
  library(Seurat)
  
  # clustering
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(object = seurat_obj))
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), verbose = F)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 1, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
    
  # pK identification 
  sweep.list <- paramSweep(seurat_obj, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn_keloid <- find.pK(sweep.stats) 
  mpK <- as.numeric(as.vector(bcmvn_keloid$pK[which.max(bcmvn_keloid$BCmetric)]))
  print(paste0("mpK:",mpK))
    
  # Homotypic doublet proportion
  homotypic.prop <- modelHomotypic(seurat_obj$seurat_clusters)  
  DoubletRate = ncol(seurat_obj)*8*1e-6 
  nExp_poi <- round(DoubletRate*ncol(seurat_obj))  
  # Heterotypic doublet proportion
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
  # Run doubletFinder
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
    
  # Integrate doublet prediction results
  seurat_obj@meta.data[,"DoubletFinder_result"] <- seurat_obj@meta.data[, paste("DF.classifications", 0.25, mpK, nExp_poi, sep = "_")]
  seurat_obj@meta.data$DoubletFinder_result[which(seurat_obj@meta.data$DoubletFinder_result == "Doublet" & seurat_obj@meta.data[, paste("DF.classifications", 0.25, mpK, nExp_poi.adj, sep = "_")] == "Singlet")] <- "Doublet-Low Confidience"
  seurat_obj@meta.data$DoubletFinder_result[which(seurat_obj@meta.data$DoubletFinder_result == "Doublet")] <- "Doublet-High Confidience"
  DoubletFinder_result <- seurat_obj$DoubletFinder_result
    
  return(DoubletFinder_result)
  }


# Doublet cluster definition by lineage markers
#' @param seurat_obj a Seurat object containing the "DoubletFinder_result" column in its meta.data
#' @param MajorCell a list of canonical lineage markers

scRNA_MajorCellMarkers <- function(seurat_obj,MajorCell){
  
  # Clustering
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(object = seurat_obj))
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), verbose = F)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 1, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
  
  # Calculate the proportion of potential doublets per cluster
    mat <- table(seurat_obj$seurat_clusters,seurat_obj$DoubletFinder_result)
    per_mat <- matrix(nrow = dim(mat)[1], ncol = dim(mat)[2])
    for(i in 1:dim(mat)[1]){
      for(j in 1:dim(mat)[2]){
        per_mat[i,j] <- round(mat[i,j]/sum(mat[i,]), 2) 
      }
    }
    colnames(per_mat) <- colnames(mat)
    rownames(per_mat) <- rownames(mat)
    write.table(per_mat, "doublet.percentage.txt", quote = F, sep = "\t", col.names = T,row.names = T)
  
  # Visualize clusters, potential doublets on UMAP, and check canonical lineage marker expression via a dot plot.
  pdf("preprocess_cluster.pdf",width = 22,height = 20)
  p1 <- UMAPPlot(seurat_obj, group.by ="seurat_clusters", raster=FALSE, label=T, label.size = 6)
  p2 <- UMAPPlot(seurat_obj, group.by ="DoubletFinder_result", cols =c("red","black","grey"), raster=FALSE) 
  p3 <- DotPlot(seurat_obj, features = MajorCell)+
    RotatedAxis()+
    scale_x_discrete("")+
    scale_y_discrete("")+
   
  print((p1|p2)/p3 + plot_layout(heights = c(1, 2)))
  dev.off()
  
  return(seurat_obj)
}


# Remove genes expressed in fewer than 5 cells and mitochondrial genes
# seurat_obj: a Seurat object

seurat_obj <- seurat_obj[which(rowSums(seurat_obj[["RNA"]]$counts != 0) >= 5),]
seurat_obj <- seurat_obj[-which(substr(rownames(seurat_obj), 1, 3) == "MT-"),]


# Normalization
# seurat_obj: a Seurat object
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
















