# Target gene comparison

# Obtain drug target genes
# result_drugs: differential sensitivity drug
# Dtar_Dban: drug target genes in the DrugBank database

tar_list <- list()
for (drug in result_drugs) {
  tar_list <- unique(Dtar_Dban$target_gene_name[grep(paste0("^",drug,"$"), x = Dtar_Dban$drug_name, ignore.case = T)])
}
names(tar_list) <- names(result_drugs)


# Difference examination
# validation_data: 8 processed HNSCC transcriptome datasets
# target: drug target genes in tar_list
# test_samples/control_samples: group patients between test and control according to ssGSEA top 3 score

wilcox_res <- wilcox.test(x = unlist(validation_data[target, test_samples]),
                          y = unlist(validation_data[target, control_samples]),
                          conf.int = T,
                          alternative = "greater")
wilcox_p <- wilcox_res$p.value
wilcox_HL <- wilcox_res$estimate
