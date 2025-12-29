# GO/KEGG/HALLMARK enrichment analysis
# gene: each signature gene set of CCS

library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

data(geneList, package = "DOSE") 
gene_df <- bitr(gene,
                fromType = "SYMBOL", 
                toType = c("ENSEMBL", "ENTREZID"), 
                OrgDb = org.Hs.eg.db)

# GO
enrich_GO <- enrichGO(gene = gene_df$ENTREZID, 
                      universe = names(geneList), 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID",
                      ont = "BP",
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)


# KEGG
enrich_KEGG <- enrichKEGG(gene = gene_df$ENTREZID,   
                          organism = "hsa", 
                          keyType = "kegg", 
                          minGSSize = 1, 
                          maxGSSize = 500,
                          pvalueCutoff = 0.05,  
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)


# HALLMARK
library(dplyr)
library(msigdbr)
Hsa_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
Hsa_hallmark_TERM2GENE <- Hsa_hallmark[, c("gs_name", "gene_symbol")]

universe_symbols <- mapIds(org.Hs.eg.db, 
                           keys = names(geneList), 
                           column = "SYMBOL", 
                           keytype = "ENTREZID")
enrich_Hallmark <- enricher(gene = gene, 
                            TERM2GENE = Hsa_hallmark_TERM2GENE,
                            universe = unname(na.omit(universe_symbols)), 
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05
)

