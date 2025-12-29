# inferCNV analysis

# Run inferCNV
#' @param seurat_obj a Seurat object containing the "celltype" column in its meta.data 
#' @param gene_order gene annotation file path
#' @param ref_group_name a character vector specifying which cell types in the "celltype" column are designated as the reference
#' @param out_dir output path

scRNA_inferCNV <-function(seurat_obj, 
                          gene_order, 
                          ref_group_name,
                          out_dir){
  library(infercnv)
  library(Seurat)
  matrix_counts <- as.matrix(seurat_obj[["RNA"]]$counts)
  
  # create infercnv object
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = matrix_counts,
                                       annotations_file = data.frame(seurat_obj$celltype), 
                                       gene_order_file = gene_order, 
                                       ref_group_names = ref_group_name)
  
  # run inferCNV
  infercnv_obj <- infercnv::run(infercnv_obj, 
                                cutoff = 0.1, 
                                out_dir = out_dir,
                                denoise = TRUE,
                                analysis_mode = "samples",
                                cluster_by_groups = FALSE, 
                                k_obs_groups = 1,
                                HMM = FALSE,
                                write_phylo = FALSE,
                                write_expr_matrix = FALSE)
}


# Generate cluster metric plots for epithelial cluster
# infercnv_obj: a infercnv object after running scRNA_inferCNV
# geneFile: a dataframe containing the gene_name (gene symbol) and chr (chromosome) columns

# Load R packages
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Import CNV profiles and construct annotations
expr <- infercnv_obj@expr.data
normal_loc <- unlist(infercnv_obj@reference_grouped_cell_indices)
test_loc <- unlist(infercnv_obj@observation_grouped_cell_indices)
anno.df <- data.frame(
  CB = c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class = c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
)
head(anno.df)
gn <- rownames(expr)
rownames(geneFile) <-  geneFile$gene_name
sub_geneFile <-  geneFile[intersect(gn,geneFile$gene_name),]
expr <- expr[intersect(gn,geneFile$gene_name),]
sub_geneFile$chr <- paste("chr",sub_geneFile$chr,sep = "")

# Cluster cells
set.seed(123)
kmeans.result <- kmeans(t(expr), 6)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB <- rownames(kmeans_df)
kmeans_df <- kmeans_df%>%inner_join(anno.df,by="CB") 
kmeans_df_s <- arrange(kmeans_df,kmeans_class) 
rownames(kmeans_df_s) <- kmeans_df_s$CB
kmeans_df_s$CB <- NULL
kmeans_df_s$kmeans_class <- as.factor(kmeans_df_s$kmeans_class) 

# Define heatmap annotations and color matching
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v <- RColorBrewer::brewer.pal(8, "Dark2")[1:6] 
names(color_v) <- as.character(1:6)
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("test"="red","normal" = "blue"),kmeans_class=color_v))

# Plot
pdf("infercnv.pdf",width = 15,height = 10)
Heatmap(t(expr)[rownames(kmeans_df_s),], 
        col = colorRamp2(c(0.85,1,1.15), c("#377EB8","#F0F0F0","#E41A1C")), 
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        column_split = factor(sub_geneFile$chr, paste("chr",1:22,sep = "")), 
        column_gap = unit(2, "mm"),
        heatmap_legend_param = list(title = "Modified expression",
                                    direction = "vertical",
                                    title_position = "leftcenter-rot",
                                    at = c(0.4,1,1.6),
                                    legend_height = unit(3, "cm")),
        top_annotation = top_anno,
        left_annotation = left_anno,
        row_title = NULL,
        column_title = NULL)
dev.off()


