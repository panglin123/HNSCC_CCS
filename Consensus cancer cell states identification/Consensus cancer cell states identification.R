# Identify consensus cancer cell states
# seurat_obj: a Seurat object of malignant cells after running Harmony

library(Seurat)

# Differential analysis

levels(seurat_obj)
markers <- FindAllMarkers(seurat_obj, group.by = "seurat_clusters", only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
markers <- markers[markers$p_val_adj < 0.05,] 


# Merge clusters from Harmony into Meta-clusters (MCs)
# harmony_cluster_name: a character vector of Harmony cluster names (C0 to C47)
# ncol: the number of Harmony clusters, ncol = 48

library(dplyr)
library(tibble)
markers_sorted <- markers %>% 
  group_by(cluster) %>% 
  arrange(desc(abs(avg_log2FC)), p_val_adj, .by_group = TRUE)
processed_list <- markers_sorted %>% group_by(cluster) %>%
  group_map(~ {
    if(nrow(.x) >= 50) {
      head(.x, 50)
    } else {
      .x 
    }
  }) %>% 
  set_names(unique(markers_sorted$cluster))
deg_df <- as.data.frame(matrix(NA, nrow = 50, ncol = ncol))
for(i in 1:length(processed_list)){
  deg_df[,i] <- processed_list[[i]]$gene
}
colnames(deg_df) <- harmony_cluster_name
deg_intersect <- apply(deg_df , 2, function(x) apply(deg_df , 2, function(y) length(intersect(x,y)))) 

Min_intersect_initial <- 10  
Min_intersect_cluster <- 10
Min_group_size <- 0   
Sorted_intersection <- sort(apply(deg_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1) ) , decreasing = TRUE)

Cluster_list <- list()  
MC_list <- list()
k <- 1
Curr_cluster <- c()
deg_intersect_original  <- deg_intersect

while (Sorted_intersection[1]>Min_group_size) {  
  
  Curr_cluster <- c(Curr_cluster, names(Sorted_intersection[1]))

  Genes_MC <- deg_df[, names(Sorted_intersection[1])]

  deg_df <- deg_df[, -match(names(Sorted_intersection[1]), colnames(deg_df))]  

  Intersection_with_Genes_MC <- sort(apply(deg_df, 2, function(x) length(intersect(Genes_MC,x))) , decreasing = TRUE) 

  NMF_history <- Genes_MC 
  
  while (Intersection_with_Genes_MC[1] >= Min_intersect_cluster) {  
    
    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MC)[1])
    Genes_MC_temp   <- sort(table(c(NMF_history , deg_df[,names(Intersection_with_Genes_MC)[1]])), decreasing = TRUE)  
    Genes_at_border <- Genes_MC_temp[which(Genes_MC_temp == Genes_MC_temp[50])]  
    
    library(dplyr)
    library(tibble)
    if (length(Genes_at_border)>1){
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        Q <- processed_list[[i]] %>% 
          filter(gene %in% names(Genes_at_border)) %>% 
          select(gene, avg_log2FC) %>% 
          deframe()
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score, Q)
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score, decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   
      
      Genes_MC_temp <- c(names(Genes_MC_temp[which(Genes_MC_temp > Genes_MC_temp[50])]), names(Genes_curr_NMF_score_sort))
      
    } else {
      Genes_MC_temp <- names(Genes_MC_temp)[1:50] 
    }
    
    NMF_history <- c(NMF_history, deg_df[,names(Intersection_with_Genes_MC)[1]]) 
    Genes_MC <- Genes_MC_temp[1:50]
    deg_df <- deg_df[,-match(names(Intersection_with_Genes_MC)[1], colnames(deg_df))]  
    Intersection_with_Genes_MC <- sort(apply(deg_df, 2, function(x) length(intersect(Genes_MC,x))) , decreasing = TRUE) 
  }
  
  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MC_list[[paste0("MC_",k)]] <- Genes_MC
  k <- k+1
  deg_intersect <- deg_intersect[-match(Curr_cluster,rownames(deg_intersect)), -match(Curr_cluster,colnames(deg_intersect))]  
  Sorted_intersection <- sort(apply(deg_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)), decreasing = TRUE)   
  
  Curr_cluster <- c()
  print(dim(deg_intersect)[2])
}












