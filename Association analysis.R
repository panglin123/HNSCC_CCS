# Association analysis: pearson
# ssGSEA_Res_A: a ssGSEA score matrix, rows are samples, columns are CCSs and Cluster
# ssGSEA_Res_B: a ssGSEA score matrix, rows are samples, columns are CCSs and Cluster

library(stats)
# calculate the average expression of CSSs  
ssGSEA_Res_A_mean <- aggregate(. ~ Cluster, data = ssGSEA_Res_A, FUN = mean)
ssGSEA_Res_B_mean <- aggregate(. ~ Cluster, data = ssGSEA_Res_B, FUN = mean)

# Pearson
corr <- cor(t(ssGSEA_Res_A_mean), t(ssGSEA_Res_B_mean), method = "pearson")




