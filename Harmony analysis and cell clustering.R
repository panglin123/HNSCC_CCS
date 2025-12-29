# seurat_obj: a Seurat object of malignant cells containing the "Sample" column in its meta.data

# Run harmony

library(Seurat)
library(harmony)
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "Sample", plot_convergence = T)

# Cell clustering

seurat_obj <- FindNeighbors(seurat_obj, dims =1:43, reduction = "harmony")
seurat_obj <- FindClusters(seurat_obj, resolution = 2)
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:43, verbose = FALSE)


