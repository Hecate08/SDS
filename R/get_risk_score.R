#' Risk Score
#'
#' This function computes the risk score for a ridge regression model with the selected features.
#' @param features_training Matrix of selected (meta-)feature values with training sample names as row names.
#' @param features_test Matrix of selected (meta-)feature values with test sample names as row names.
#' @param survObject_training The survival data of the training sample in form of a survival object
#' @param alpha Alpha value which declares ridge regression (alpha = 0). For lasso set alpha = 1. Defaults to 0.
#' @keywords risk prediction
#' @return List of one risk score per sample for training samples and test samples
#' @export
#' @examples
#' get_risk_score()
get_risk_score <- function(features_training,features_test,survObject_training,alpha=0){
  riskScoreTest <- matrix(NA, ncol = 10, nrow = dim(features_test)[1])
  riskScoreTraining <- matrix(NA, ncol = 10, nrow = dim(features_training)[1])
  for(f in 1:10){

    cvglmfit <- cv.glmnet(features_training,survObject_training,family = "cox", nfolds = 5, alpha = alpha)

    glmPredTest <- predict(object = cvglmfit, newx = features_test, s = "lambda.min", type = "link")
    riskScoreTest[,f] <- exp(glmPredTest)

    glmPredTraining <- predict(object = cvglmfit, newx = features_training, s = "lambda.min", type = "link")
    riskScoreTraining[,f] <- exp(glmPredTraining)

  }
  meanRiskScore <- list(training = rowMeans(riskScoreTraining,na.rm=TRUE), test = rowMeans(riskScoreTest,na.rm=TRUE))

  return(meanRiskScore)
}
