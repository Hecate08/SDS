# functions

Score1Function <- function(filteredExp, survv){
  
  RowVar <- function(x, ...) {
    rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
  }
  sigma2 = RowVar(filteredExp)
  # sorted timepoints of events
  timepoints <- sort(unique(survv[survv$overall.survival.indicator==1,]$OVERALL.SURVIVAL))
  # vector of weights for each gene
  sds <- rep(0, each=dim(filteredExp)[1])
  
  for (t in timepoints){
    m <- rowMeans(subset(filteredExp, select = rownames(survv[survv$OVERALL.SURVIVAL >= t,])))
    m[is.nan(m)] = 0 # this transforms all NaN into 0
    n = dim(subset(filteredExp, select = rownames(survv[survv$OVERALL.SURVIVAL >= t,])))[2]
    
    
    mh <- rowMeans(subset(filteredExp, select = rownames(survv[survv$OVERALL.SURVIVAL <= t & survv$overall.survival.indicator==1,])))
    mh[is.nan(mh)] = 0 # this transforms all NaN into 0
    nh = dim(subset(filteredExp, select = rownames(survv[survv$OVERALL.SURVIVAL <= t & survv$overall.survival.indicator==1,])))[2]
    
    # function
    # look only at patients with event at time t
    patients <- rownames(survv[survv$OVERALL.SURVIVAL == t & survv$overall.survival.indicator == 1,])
    
    w1 = rep(0, each=dim(filteredExp)[1])
    for (p in patients){
      w1 = w1 + n*(filteredExp[,p] - m)^2 + nh*(filteredExp[,p] - mh)^2
    }
    w1 = w1/length(patients)
    sds <- sds + w1
    
  }
  sds <- sds/sigma2
  return(sds)
}

Score2Function <- function(trainingExp,trainingClin){
  lenE = dim(trainingExp)[1]
  R2adj = rep(NA,lenE)
  for(f in 1:lenE){
    fit = lm(trainingExp[f,] ~ .,data = trainingClin)
    R2adj[f] = summary(fit)$adj.r.squared
  }
  R2adj_ok = as.numeric(gsub(NaN, 1, R2adj))
  names(R2adj_ok) = rownames(trainingExp)
  
  score2 = 1 - R2adj_ok
  
  return(score2)
}

scoreCombination = function(score1,score2,weight){
  #print(combo)
  # normalizing
  mu = mean(score1)
  std = sqrt(var(score1))
  score1_normal = (score1 - mu)/std
  names(score1_normal) = names(score1)
  #plot(density(score1_normal))
  
  mu = mean(score2[score2 != 0])
  std = sqrt(var(score2[score2 != 0]))
  score2_normal = (score2 - mu)/std
  names(score2_normal) = names(score2)
  
  #top 50 genes filtered by correlation
  sum_normal = score1_normal + score2_normal*weight
  names(sum_normal) = names(score1)
  
  return(sum_normal)
}

corrGroupFilter = function(expression, expression_test, corrT, top, score1){
  corrMatrix = cor(expression, method = "spearman")
  diag(corrMatrix) = 0
  lenC = dim(corrMatrix)[1]
  expression_corr = expression
  
  delete = c()
  newFeatures = matrix(NA, ncol = top, nrow = dim(expression_corr)[1])
  rownames(newFeatures) = rownames(expression_corr)
  newFeatures_test = matrix(NA, ncol = top, nrow = dim(expression_test)[1])
  rownames(newFeatures_test) = rownames(expression_test)
  n = 1
  for(f in 1:lenC){
    
    name = rownames(corrMatrix)[f]
    filteredCorrMatrixNames = rownames(corrMatrix)
    if(length(delete)>0){
      filteredCorrMatrixNames = rownames(corrMatrix)[- which(rownames(corrMatrix) %in% delete)]
    }
    if(is.element(name, filteredCorrMatrixNames)){
      filteredCorrMatrix = corrMatrix[filteredCorrMatrixNames,filteredCorrMatrixNames]
      geneTrue = filteredCorrMatrix[name,]>corrT
      newFeatures[,n] = expression_corr[,name]
      newFeatures_test[,n] = expression_test[,name]
      
      if(sum(geneTrue) > 0){
        groupNames = c(name, names(geneTrue[geneTrue]))
        delete = c(delete, groupNames)
        
        scoresWeights = score1[groupNames]/sum(score1[groupNames])
        newFeatures[,n] = rowSums(t(t(expression_corr[,groupNames]) * scoresWeights))
        newFeatures_test[,n] = rowSums(t(t(expression_test[,groupNames]) * scoresWeights))
      }
      n = n + 1
    }
  }
  
  newFeatures = newFeatures[, !apply(is.na(newFeatures), 2, all)]
  newFeatures_test = newFeatures_test[, !apply(is.na(newFeatures_test), 2, all)]
  
  features = list(training = newFeatures, test = newFeatures_test)
  return(features)
}

meanRiskScoreFunction = function(dataTrai,dataTe,survObject,alp){
  riskScoreTest = matrix(NA, ncol = 10, nrow = dim(dataTe)[1])
  riskScoreTraining = matrix(NA, ncol = 10, nrow = dim(dataTrai)[1])
  for(f in 1:10){
    #print(f)
    cvglmfit = cv.glmnet(dataTrai,survObject,family = "cox", nfolds = 5, alpha = alp)
    
    glmPredTest = predict(object = cvglmfit, newx = dataTe, s = "lambda.min", type = "link")
    riskScoreTest[,f] = exp(glmPredTest)
    #riskScoreTest[,f] = glmPredTest
    
    glmPredTraining = predict(object = cvglmfit, newx = dataTrai, s = "lambda.min", type = "link")
    riskScoreTraining[,f] = exp(glmPredTraining)
    #riskScoreTraining[,f] = glmPredTraining
    
  }
  meanRiskScore = list(training = rowMeans(riskScoreTraining,na.rm=TRUE), test = rowMeans(riskScoreTest,na.rm=TRUE))
  
  return(meanRiskScore)
}

meanCIndexFunction = function(dataTrai,dataTe,te,trai,ty,surv,survObject,al){
  cIndexeS1 = rep(NA,10)
  for(f in 1:10){
    #print(f)
    cvglmfit = cv.glmnet(dataTrai,survObject,family = "cox", nfolds = 5, alpha = al)
    
    if(ty == "test"){
      glmPred = predict(object = cvglmfit, newx = dataTe, s = "lambda.min", type = "link")
      cindex_validation = concordance.index(glmPred, surv.time = surv[te,]$OVERALL.SURVIVAL,
                                            surv.event = surv[te,]$overall.survival.indicator, 
                                            method = "noether")
      cIndexeS1[f] = cindex_validation$c.index
    } else {
      glmPred = predict(object = cvglmfit, newx = dataTrai, s = "lambda.min", type = "link")
      cindex_validation = concordance.index(glmPred, surv.time = surv[trai,]$OVERALL.SURVIVAL,
                                            surv.event = surv[trai,]$overall.survival.indicator, 
                                            method = "noether")
      cIndexeS1[f] = cindex_validation$c.index
      
    } 
  }
  meanCIndex = mean(cIndexeS1,na.rm=TRUE)
  return(meanCIndex)
}
