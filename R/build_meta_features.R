#' Feature reduction
#'
#' This function computes meta features from the best selected features by building a weighted mean of a group of pairwise (spearman) correlated features.
#' @param attrib_training Matrix of atrribute values with attribute names as row names and training sample names as column names.
#' @param attrib_test Matrix of atrribute values with attribute names as row names and test sample names as column names.
#' @param corrT Threshold to decide whether two features are correlated. Defaults to 0.6.
#' @param top Vector of selected attribute names.
#' @param scores Either score1 or score_combo. Vector of survival distance scores or combo with attribute names as names.
#' @keywords feature reduction
#' @keywords correlation
#' @return Matrix of meta feature values with sample names as row names.
#' @export
#' @examples
#' build_meta_feature()
build_meta_features <- function(attrib_training, attrib_test, corrT = 0.6, top, scores){
  corrMatrix <- cor(attrib_training, method = "spearman")
  diag(corrMatrix) <- 0
  lenC <- dim(corrMatrix)[1]
  attrib_training_corr <- attrib_training

  delete <- c()
  newFeatures <- matrix(NA, ncol = top, nrow = dim(attrib_training_corr)[1])
  rownames(newFeatures) <- rownames(attrib_training_corr)
  newFeatures_test <- matrix(NA, ncol = top, nrow = dim(attrib_test)[1])
  rownames(newFeatures_test) <- rownames(attrib_test)
  n <- 1
  for(f in 1:lenC){

    name <- rownames(corrMatrix)[f]
    filteredCorrMatrixNames <- rownames(corrMatrix)
    if(length(delete)>0){
      filteredCorrMatrixNames <- rownames(corrMatrix)[- which(rownames(corrMatrix) %in% delete)]
    }
    if(is.element(name, filteredCorrMatrixNames)){
      filteredCorrMatrix <- corrMatrix[filteredCorrMatrixNames,filteredCorrMatrixNames]
      geneTrue <- filteredCorrMatrix[name,]>corrT
      newFeatures[,n] <- attrib_training_corr[,name]
      newFeatures_test[,n] <- attrib_test[,name]

      if(sum(geneTrue) > 0){
        groupNames <- c(name, names(geneTrue[geneTrue]))
        delete <- c(delete, groupNames)

        scoresWeights <- scores[groupNames]/sum(scores[groupNames])
        newFeatures[,n] <- rowSums(t(t(attrib_training_corr[,groupNames]) * scoresWeights))
        newFeatures_test[,n] <- rowSums(t(t(attrib_test[,groupNames]) * scoresWeights))
      }
      n <- n + 1
    }
  }

  newFeatures <- newFeatures[, !apply(is.na(newFeatures), 2, all)]
  newFeatures_test <- newFeatures_test[, !apply(is.na(newFeatures_test), 2, all)]

  features <- list(training = newFeatures, test = newFeatures_test)
  return(features)
}
