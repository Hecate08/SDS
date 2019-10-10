# libraries and arguments
args <- commandArgs(TRUE)
i = args[1]
disease = args[2]

path = paste0("~/projects/survivalOmics/data/",disease,"/")


expressio = readRDS(file = paste0(path,"filteredExp_",disease))
filteredClinDummy = readRDS(file = paste0(path,"filteredClin_dummy_",disease))
survv = readRDS(file = paste0(path,"survival_",disease))
trainingSize2 = floor((dim(survv)[1]/3)*2)



library(survival)
#install.packages("glmnet")
library(glmnet)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("survcomp", version = "3.8")
library(survcomp)
library(matrixStats)


# functions

survivalFisherScore1 <- function(filteredExp, survv){
  
  survObjectTraining <- Surv(survv$OVERALL.SURVIVAL,
                             survv$overall.survival.indicator) 
  lengthE = dim(filteredExp)[1]
  sim = rep(NA, lengthE)
  for(i in 1:lengthE){
    gene = filteredExp[i,]
    genefit = summary(coxph(survObjectTraining ~ gene))
    sim[i] = abs(genefit$coefficients[,"z"])
  }
  return(sim)
}

survivalReliefScore2 <- function(trainingExp,trainingClin){
  lenE = dim(trainingExp)[1]
  R2adj = rep(NA,lenE)
  for(f in 1:lenE){
    fit = lm(trainingExp[f,] ~ .,data = trainingClin)
    R2adj[f] = summary(fit)$adj.r.squared
  }
  R2adj_ok = as.numeric(gsub(NaN, 1, R2adj))
  names(R2adj_ok) = rownames(trainingExp)
  
  score2 = 1 - R2adj_ok
  #score2 = R2adj_ok
  return(score2)
}
corrGroupFilter = function(expres_top50, expres_top50_test, corrT, top, score1){
  corrMatrix = cor(expres_top50, method = "spearman")
  diag(corrMatrix) = 0
  lenC = dim(corrMatrix)[1]
  expres_top50_corr = expres_top50
  
  delete = c()
  newFeatures = matrix(NA, ncol = top, nrow = dim(expres_top50_corr)[1])
  rownames(newFeatures) = rownames(expres_top50_corr)
  newFeatures_test = matrix(NA, ncol = top, nrow = dim(expres_top50_test)[1])
  rownames(newFeatures_test) = rownames(expres_top50_test)
  n = 1
  for(f in 1:lenC){
    #print(f)
    name = rownames(corrMatrix)[f]
    filteredCorrMatrixNames = rownames(corrMatrix)
    if(length(delete)>0){
      filteredCorrMatrixNames = rownames(corrMatrix)[- which(rownames(corrMatrix) %in% delete)]
    }
    if(is.element(name, filteredCorrMatrixNames)){
      filteredCorrMatrix = corrMatrix[filteredCorrMatrixNames,filteredCorrMatrixNames]
      geneTrue = filteredCorrMatrix[name,]>corrT
      newFeatures[,n] = expres_top50_corr[,name]
      newFeatures_test[,n] = expres_top50_test[,name]
      
      if(sum(geneTrue) > 0){
        groupNames = c(name, names(geneTrue[geneTrue]))
        delete = c(delete, groupNames)
        
        scoresWeights = score1[groupNames]/sum(score1[groupNames])
        newFeatures[,n] = rowSums(t(t(expres_top50_corr[,groupNames]) * scoresWeights))
        newFeatures_test[,n] = rowSums(t(t(expres_top50_test[,groupNames]) * scoresWeights))
      }
      n = n + 1
    }
  }
  
  newFeatures = newFeatures[, !apply(is.na(newFeatures), 2, all)]
  newFeatures_test = newFeatures_test[, !apply(is.na(newFeatures_test), 2, all)]
  
  #  topnew = top
  #  if(top > length(filteredCorrMatrixNames)){
  #    print(paste("topnew:", length(filteredCorrMatrixNames)))
  #    topnew = length(filteredCorrMatrixNames)
  #  }
  #  score1_new2 = sort(score1[filteredCorrMatrixNames],decreasing = TRUE)
  
  #  score1_expres_top50 = expres_top50[,names(score1_new2[1:topnew])]
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



# variables
expFil = "cor"
corrT = 0.6
top = 75
cIndexUnP = rep(NA, 2)
print(paste("seed",i))
set.seed(i)
training <- sample(rownames(survv), size = trainingSize2)
test <- rownames(survv)[!(rownames(survv) %in% training)]
survObjectTraining <- Surv(survv[training,]$OVERALL.SURVIVAL,
                           survv[training,]$overall.survival.indicator)  

# scores
print("scores")
score1 = survivalFisherScore1(filteredExp = expressio[,training], survv = survv[training,])
names(score1) = rownames(expressio)
score2 = survivalReliefScore2(trainingExp = expressio[,training], trainingClin = filteredClinDummy[training,])
names(score2) = rownames(expressio)


################################
### score 1 alone 
################################
print("sds")
# score1 top genes
score1_top50_temp = sort(score1, decreasing = TRUE)[1:top]
topGenes = names(score1_top50_temp)
expres_top50 = t(expressio[topGenes,training])
expres_top50_test = t(expressio[topGenes,test])
score1_expres_top50 = list(training = expres_top50[,1:top], expres_top50_test[,colnames(expres_top50[,1:top])])
if(expFil == "cor"){
  score1_expres_top50 = corrGroupFilter(expres_top50 = expres_top50, expres_top50_test = expres_top50_test, corrT = corrT, top = top, score1 = score1)
}



#######################
### score1 c.Index     ###
dataTrainingGLM <- score1_expres_top50$training
dataTestGLM <- score1_expres_top50$test

cIndexUnP[1] <- meanCIndexFunction(dataTrai = dataTrainingGLM,
                                   dataTe = dataTestGLM,
                                   te = test,
                                   trai = training,
                                   ty = "test",
                                   surv = survv,
                                   survObject = survObjectTraining, al = 0)








################################
### score 1 + clin
################################
print("sds + clin")
#######################
### score1 risk     ###
dataTrainingGLM <- score1_expres_top50$training
dataTestGLM <- score1_expres_top50$test
riskScoreList <- meanRiskScoreFunction(dataTrai = dataTrainingGLM,
                                       dataTe = dataTestGLM, 
                                       survObject = survObjectTraining, 
                                       alp = 0)


#######################
### score1 + clin   ###
dataTraining <- data.frame(filteredClinDummy[training,], riskScore = riskScoreList$training)
dataTest <- data.frame(filteredClinDummy[test,], riskScore = riskScoreList$test)
cIndexeS1 = rep(NA,10)
for(f in 1:10){
  if(disease == "bladder"){
    cfit1 <- coxph(survObjectTraining ~ age + gender_male + tumor_stage_low + riskScore, dataTraining)
  } else {
    cfit1 <- coxph(survObjectTraining ~ age + gender_male + tumor_stage_stage_iii + tumor_stage_stage_ii + tumor_stage_stage_iv + riskScore, dataTraining)
  }
  
  pred1Test <- predict(object = cfit1, newdata = dataTest, type = "lp")
  cindex_validation = concordance.index(pred1Test, surv.time = survv[test,]$OVERALL.SURVIVAL,
                                        surv.event = survv[test,]$overall.survival.indicator, 
                                        method = "noether")
  cIndexeS1[f] = cindex_validation$c.index
  
}
cIndexUnP[2] <- mean(cIndexeS1,na.rm=TRUE)


saveRDS(cIndexUnP, file = paste0(disease,"_SDSunp_simple_spearman_seed",i,".rds"))

