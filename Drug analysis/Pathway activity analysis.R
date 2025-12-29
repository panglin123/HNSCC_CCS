# Load R packages
library(clusterProfiler)
library(org.Hs.eg.db)

# Obtain KEGG pathways related to target genes
# state_target: Significantly up-regulated CCS target genes

kegg_data <- download_KEGG("hsa")
for (gene in state_target) {
  Etz_id <- mapIds(x = org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "SYMBOL")
  
  hsa_path_id <- kegg_data$KEGGPATHID2EXTID[which(kegg_data$KEGGPATHID2EXTID$to == Etz_id),"from"]
  
  pathway_gene <- lapply(hsa_path_id, function(x){
    filter(kegg_data$KEGGPATHID2EXTID, from %in% x)
  })
  
  target_pathway <- sapply(hsa_path_id, function(x){
    kegg_data$KEGGPATHID2NAME$to[kegg_data$KEGGPATHID2NAME$from == x]
    })
}


# Differential Testing of Target Pathway Activity
# pathway_GSEAscore: activity of each patient in various target pathways
# test_samples/control_samples: group patients between test and control according to ssGSEA top 3 score

temp_wilcox <-  apply(pathway_GSEAscore, 1, function(x){
  wilcox.test(x[test_samples], x[control_samples], conf.int = T, alternative = "greater")[c("p.value", "estimate")]
})
