# ConsensusCluster analysis
# ssGSEA_Res: a ssGSEA score matrix, rows are CCSs, columns are samples

ssGSEA_Res = sweep(ssGSEA_Res, 1, apply(ssGSEA_Res, 1, median, na.rm = T))

library(ConsensusClusterPlus)
title=tempdir()
results = ConsensusClusterPlus(d = ssGSEA_Res, 
                               maxK = 20, 
                               reps = 500, 
                               pItem = 0.8, 
                               pFeature = 1, 
                               title = "Consensus_cluster", 
                               clusterAlg = "hc", 
                               distance = "pearson", 
                               seed = 123, 
                               tmyPal = NULL, 
                               writeTable = FALSE, 
                               plot = "pdf")