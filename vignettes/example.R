library(pROC)

Rcpp::sourceCpp('src/dmvnrm_arma.cpp')
Rcpp::sourceCpp('src/dmvnorm_arma.cpp')
source('SAMGEP.R')
source('simulate.R')

## Test just E step (for C++ development)

dat <- readRDS('Estep_testDat_100pats.rds')
fitted <- readRDS('Estep_fittedModel.rds')

system.time(Estep_result <- Estep_partial(dat,fitted,nX=10))
print(paste('Sanity Check: AUC =',auc(dat$Y,Estep_result)))

## Test entire procedure 

nTrain <- 1000
nTest <- 100
embedDim <- 10

simulated <- simulate(trainData[[1]],trainData[[2]],6:160,nTrain,nTest,embedDim)
samgep_result <- samgep(dat_train=simulated$trainSet,dat_test=simulated$testSet,Cindices=6:160,
                        V=simulated$embeddings,observed=1:100,nCores=10)

idx <- which(simulated$testSet$ID %in% 101:1000)
samgep_auc <- auc(simulated$testSet$Y[idx], samgep_result$margMix[idx])
mgpsup_auc <- auc(simulated$testSet$Y[idx], samgep_result$margSup[idx])
mgpsemisup_auc <- auc(simulated$testSet$Y[idx], samgep_result$margSemisup[idx])
print(paste('AUCs:',samgep_auc,mgpsup_auc,mgpsemisup_auc))

  
