# Survival analysis
# survival_dat: a dataframe with time, status and geneset_mean columns
# col_cluster: a dataframe with ConsensusCluster analysis result

library(survival)
library(survminer)
# The average expression for each CCS signatures was calculated for survival analysis. 
survival_dat$group <- ifelse(survival_dat$geneset_mean > median(survival_dat$geneset_mean, na.rm = TRUE),"high","low")
f1 <- survfit(Surv(time, status)~ group, data = survival_dat)
pdf("CCS.survival.pdf",width = 8,height = 8)
ggsurvplot(f1, 
           data = survival_dat,
           surv.median.line = "hv",
           legend.title = "Group",
           legend.labs = c("high","low"), 
           pval = TRUE,            
           conf.int = F,
           risk.table = T)
dev.off()

# Survival analysis was conducted using CSSs.
survival_dat$group <- col_cluster
f1 <- survfit(Surv(time, status)~ group, data = survival_dat)
pdf("CSS.survival.pdf",width = 8,height = 8)
ggsurvplot(f1, 
           data = survival_dat,
           surv.median.line = "hv",
           legend.title = "Group",
           legend.labs = c("CSS1", "CSS2", "CSS3", "CSS4", "CSS5", "CSS6"), 
           pval = TRUE,            
           conf.int = F,
           risk.table = T)
dev.off()


