library(rcompanion)
library(tidyverse)

# Oncopredict drug response prediction
#' @param testExpr gene expression profile of the sample to be predicted
#' @param Dir path to the processed drug response data folder
#' @param background_DB the background database used for prediction: CTPR, PRISM, GDSC1, GDSC2
#' @param response_type the type of drug response to predict: IC50, EC50, AUC
#' @param cancer_type the TCGA abbreviation for the type of cancer

run_oncoPREDICT <- function(testExpr, Dir, background_DB, response_type, cancer_type){
  # 1.Load processed data
  load(paste0(paste(Dir,background_DB,"processed",response_type,sep = "/"),".RData"))
  expr <- expr[,intersect(colnames(expr),unlist(unname(ccl_type[cancer_type])))]
  resp <- resp[intersect(rownames(resp),unlist(unname(ccl_type[cancer_type]))),]
  # 2.Filter drugs
  # 2.1 Remove drugs whose auc become Inf after powerTransform 
  if (any(is.infinite(exp(resp)))) {
    resp <- resp[,-as.numeric(names(table((which(is.infinite(exp(resp)),arr.ind = T))[,2])))]
  }
  # 2.2 Remove drugs that can not be processed by powerTransform 
  err_loc <- c()
  for (i in 1:ncol(resp)) {
    tryCatch(
      car::powerTransform(resp[,i]),
      error = function(x){
        err_loc<<-c(err_loc,i)
        TRUE
      })}
  if(length(err_loc)!=0){
    resp <- resp[,-err_loc]
  }
  # 3.Oncopredict drug response prediction
  onc <- oncoPredict::calcPhenotype(trainingExprData = expr,
                                  trainingPtype = resp,
                                  testExprData = testExpr,
                                  batchCorrect = 'eb', 
                                  powerTransformPhenotype = T,
                                  removeLowVaryingGenes = 0.2,
                                  minNumSamples = 10, 
                                  printOutput = TRUE, 
                                  folder = F,
                                  removeLowVaringGenesFrom = 'rawData' )
  return(onc)
}


# Differential sensitivity examination
#' @param Response_profiles object after running run_oncoPREDICT
#' @param groups patient groups, such as CSS or CCS
#' @param label tested label, such as CSS1 or EpiSen_I

resp_wilcox_pseu <- function(Response_profiles, groups, label){
  resp_dif <- apply(Response_profiles, 2, function(x){
    wilcox.test(x[groups$label == label],
                x[groups$label != label], alternative = "less")[c("statistic","p.value")]
  })
  resp_dif <- as.data.frame(t(sapply(resp_dif, "[", i = 1:max(sapply(resp_dif, length)))))
  resp_dif$p.adjust <- p.adjust(resp_dif$p.value, method = "fdr")
  return(resp_dif)
}

# Pay attention to the calculation order of Cliff's Delta 
temp_Cdelta <- apply(Response_profiles, 2, function(x){
  cliffDelta(x = x[groups$label != label], 
             y = x[groups$label == label])
})

# Filter drugs
resp_dif$Cliffd <- unlist(temp_Cdelta[[label]])
resp_dif <- resp_dif[resp_dif$p.adjust < 0.05 & resp_dif$Cliffd > 0.33,]

# label_support_table：statistics of supporting drug database for each label
result_drugs <- label_support_table[, colSums(label_support_table) >= 2, drop = F]

