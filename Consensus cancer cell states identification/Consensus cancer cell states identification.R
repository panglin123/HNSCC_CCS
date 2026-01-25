# Identify consensus cancer cell states
# seurat_obj: a Seurat object of malignant cells after running Harmony

library(Seurat)

# Differential analysis

levels(seurat_obj)
Markers <- FindAllMarkers(seurat_obj, group.by = "seurat_clusters", only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
Markers <- Markers[Markers$p_val_adj < 0.05,] 


# Merge clusters from Harmony into Meta-clusters (MCs)
# harmony_cluster_name: a character vector of Harmony cluster names (C0 to C47)
# ncol: the number of Harmony clusters, ncol = 48

library(dplyr)
library(tibble)
Markers_sorted <- Markers %>% 
  group_by(cluster) %>% 
  arrange(desc(abs(avg_log2FC)), p_val_adj, .by_group = TRUE)
Processed_list <- Markers_sorted %>% group_by(cluster) %>%
  group_map(~ {
    if(nrow(.x) >= 50) {
      head(.x, 50)
    } else {
      .x 
    }
  }) %>% 
  set_names(unique(Markers_sorted$cluster))
DEG_df <- as.data.frame(matrix(NA, nrow = 50, ncol = ncol))
for(i in 1:length(Processed_list)){
  DEG_df[,i] <- Processed_list[[i]]$gene
}
colnames(DEG_df) <- harmony_cluster_name
DEG_intersect <- apply(DEG_df , 2, function(x) apply(DEG_df , 2, function(y) length(intersect(x,y)))) 

Min_intersect_initial <- 10  
Min_intersect_cluster <- 10
Min_group_size <- 0   
Sorted_intersection <- sort(apply(DEG_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1) ) , decreasing = TRUE)

Cluster_list <- list()  
MC_list <- list()
k <- 1
Curr_cluster <- c()
DEG_intersect_original  <- DEG_intersect

while (Sorted_intersection[1]>Min_group_size) {  
  
  Curr_cluster <- c(Curr_cluster, names(Sorted_intersection[1]))

  Genes_MC <- DEG_df[, names(Sorted_intersection[1])]

  DEG_df <- DEG_df[, -match(names(Sorted_intersection[1]), colnames(DEG_df))]  

  Intersection_with_Genes_MC <- sort(apply(DEG_df, 2, function(x) length(intersect(Genes_MC,x))) , decreasing = TRUE) 

  History <- Genes_MC 
  
  while (Intersection_with_Genes_MC[1] >= Min_intersect_cluster) {  
    
    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MC)[1])
    Genes_MC_temp   <- sort(table(c(History , DEG_df[,names(Intersection_with_Genes_MC)[1]])), decreasing = TRUE)  
    Genes_at_border <- Genes_MC_temp[which(Genes_MC_temp == Genes_MC_temp[50])]  
    
    library(dplyr)
    library(tibble)
    if (length(Genes_at_border)>1){
      Genes_curr_cluster_score <- c()
      for (i in Curr_cluster) {
        Q <- Processed_list[[i]] %>% 
          filter(gene %in% names(Genes_at_border)) %>% 
          select(gene, avg_log2FC) %>% 
          deframe()
        Genes_curr_cluster_score <- c(Genes_curr_cluster_score, Q)
      }
      Genes_curr_cluster_score_sort <- sort(Genes_curr_cluster_score, decreasing = TRUE)
      Genes_curr_cluster_score_sort <- Genes_curr_cluster_score_sort[unique(names(Genes_curr_cluster_score_sort))]   
      
      Genes_MC_temp <- c(names(Genes_MC_temp[which(Genes_MC_temp > Genes_MC_temp[50])]), names(Genes_curr_cluster_score_sort))
      
    } else {
      Genes_MC_temp <- names(Genes_MC_temp)[1:50] 
    }
    
    History <- c(History, DEG_df[,names(Intersection_with_Genes_MC)[1]]) 
    Genes_MC <- Genes_MC_temp[1:50]
    DEG_df <- DEG_df[,-match(names(Intersection_with_Genes_MC)[1], colnames(DEG_df))]  
    Intersection_with_Genes_MC <- sort(apply(DEG_df, 2, function(x) length(intersect(Genes_MC,x))) , decreasing = TRUE) 
  }
  
  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MC_list[[paste0("MC_",k)]] <- Genes_MC
  k <- k+1
  DEG_intersect <- DEG_intersect[-match(Curr_cluster,rownames(DEG_intersect)), -match(Curr_cluster,colnames(DEG_intersect))]  
  Sorted_intersection <- sort(apply(DEG_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)), decreasing = TRUE)   
  
  Curr_cluster <- c()
  print(dim(DEG_intersect)[2])
}













