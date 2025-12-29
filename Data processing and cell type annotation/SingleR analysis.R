# SingleR analysis

# seurat_obj: a Seurat object

library(Seurat)
library(celldex)
library(SingleR)
hpca.se <- celldex::HumanPrimaryCellAtlasData()
seurat_obj_SingleR <- GetAssayData(seurat_obj, slot="data") 
hesc <- SingleR(test = seurat_obj_SingleR, ref = hpca.se, labels = hpca.se$label.main)
seurat_obj$celltype <- hesc$labels


