# libraries and arguments
library(survival)
library(glmnet)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("survcomp", version = "3.8")
library(survcomp)
library(matrixStats)
source("SDS_functions.R")



# variables
args <- commandArgs(TRUE)
i = args[1]
disease = args[2]
expFil = "cor"
corrT = 0.6
top = 75
cIndexUnP = rep(NA, 2)



# data sets
path = paste0("~/projects/survivalOmics/data/",disease,"/")

expressio = readRDS(file = paste0(path,"filteredExp_",disease))
filteredClinDummy = readRDS(file = paste0(path,"filteredClin_dummy_",disease))
survv = readRDS(file = paste0(path,"survival_",disease))
trainingSize2 = floor((dim(survv)[1]/3)*2)



# cross validation data sets
print(paste("seed",i))
set.seed(i)
training <- sample(rownames(survv), size = trainingSize2)
test <- rownames(survv)[!(rownames(survv) %in% training)]
survObjectTraining <- Surv(survv[training,]$OVERALL.SURVIVAL,
                           survv[training,]$overall.survival.indicator)  


# scores
print("scores")
score1 = Score1Function(filteredExp = expressio[,training], survv = survv[training,])
names(score1) = rownames(expressio)
score2 = Score2Function(trainingExp = expressio[,training], trainingClin = filteredClinDummy[training,])
names(score2) = rownames(expressio)



################################
### score 1 + score 2 + clin
################################

fro = 0
tro = 5
weights = seq(from=fro, to=tro, by=0.1)
cIndexe = matrix(NA, nrow = length(weights), ncol = 2)
for(weightUnp in weights){
  print(weightUnp)
  
  # scoreCombo top genes 
  topGenes = names(score1)  
  score_combo = scoreCombination(score1 = score1, score2 = score2, weightUnp)
  topGenes = names(sort(score_combo, decreasing = TRUE))[1:top]

  expression = t(expressio[topGenes,training])
  expression_test = t(expressio[topGenes,test])
  score1_expression = list(training = expression[,1:top], test = expression_test[,colnames(expression[,1:top])])
  if(expFil == "cor"){
    score1_expression = corrGroupFilter(expression = expression, expression_test = expression_test, corrT = corrT, top = top, score1 = score_combo)
  }
  
  
  
  #######################
  ### scoreCombo risk     ###
  dataTrainingGLM <- score1_expression$training
  dataTestGLM <- score1_expression$test
  riskScoreList <- meanRiskScoreFunction(dataTrai = dataTrainingGLM,
                                         dataTe = dataTestGLM, 
                                         survObject = survObjectTraining, 
                                         alp = 0)
  
  #######################
  ### scoreCombo + clin ###
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
  
  cIndexe[((weightUnp*10)-(10*fro))+1,1] = weightUnp
  cIndexe[((weightUnp*10)-(10*fro))+1,2] = mean(cIndexeS1,na.rm=TRUE)
}

saveRDS(cIndexe, file = paste0(disease,"_SDSC_computingCluster_seed",i,".rds"))
