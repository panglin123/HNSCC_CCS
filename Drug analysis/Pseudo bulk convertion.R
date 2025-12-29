# load R packages
library(xlsx)
library(Seurat)
library(tidyverse)

# Pseudo bulk convertion
# seurat_obj: a Seurat object of all patients' malignant cells containing the "Sample" and "Cell_state" columns in its meta.data 

# Convert pseudo-bulk according to Sample and Cell_state

pseudo_seu <- AggregateExpression(object = seurat_obj, slot = "counts", assays = "RNA", group.by = c("Sample","Cell_state"))
pseudo_seu <- CreateSeuratObject(counts = as.matrix(pseudo_seu$RNA), min.cells = 0, min.features = 0)
pseudo_seu <- NormalizeData(pseudo_seu, normalization.method = "LogNormalize", scale.factor = 10000)
pseudo_seu <- as.matrix(pseudo_seu@assays$RNA$data)


# Sampling pseudo sample
# Set the random seed in advance

sampled_sample <- matrix(nrow = 10)
for (state in CCS_state) {
  temp_sample <- grep(paste0(state,"$"), colnames(pseudo_seu), value = T)
  if(length(temp_sample) < 10){
    state_sample <- c(temp_sample, rep("", times = (10-length(temp_sample))))
  }else{
    state_sample <- sample(temp_sample, 10, replace = F) 
  }
  sampled_sample <- cbind(sampled_sample, state_sample)
}
sampled_sample <- sampled_sample[,-1]
sampled_sample <- as.data.frame(sampled_sample)
colnames(sampled_sample) <- pseudo_state
pseudo_bulk <- pseudo_seu[, setdiff(unlist(sampled_sample), "")]
