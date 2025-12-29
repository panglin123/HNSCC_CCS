# Hierarchical clustering analysis
# Gene_exprData: a gene expression matrix, rows are genes, columns are samples

# Calculate the Euclidean matrix between the samples
dist_samples <- dist(t(Gene_exprData), method = "euclidean") 
dist_genes <- dist(Gene_exprData, method = "euclidean") 

# Perform hierarchical clustering
hc_samples <- hclust(dist_samples, method = "ward.D2") 
hc_genes <- hclust(dist_genes, method = "average")

# Cut the clustering tree, define the sample clusters, and prepare the heatmap annotation columns
sample_clusters <- cutree(hc_samples, k = 4)  
annotation_col <- data.frame(Cluster = factor(sample_clusters))
annotation_col$Cluster <- paste("Cluster",annotation_col$Cluster, sep = "")
annotation_col$Cluster <- factor(annotation_col$Cluster, levels = paste("Cluster",1:4,sep = ""))
library(pheatmap)
p1 <- pheatmap(Gene_exprData,
               cluster_rows = TRUE,   
               cluster_cols = hc_samples, 
               clustering_distance_rows = "euclidean",
               clustering_method = "ward.D2",
               annotation_col = annotation_col,    
               show_rownames = FALSE, 
               show_colnames = FALSE,  
               color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
               scale = "row") 
