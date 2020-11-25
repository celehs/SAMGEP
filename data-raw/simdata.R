library(bindata)

simulate <- function(trainDat,totDat,Xinds,nTrain,nTest,embedDim,rho=0.4,
                     yt_gen='correct',xy_gen='lognormal',nInformative=NA){
  
  # Fit P(Y0), Yt | Yt-1,T,HID
  keep <- which(trainDat$T != 0)
  tDiff1 <- as.integer(trainDat$T[keep]-trainDat$T[keep-1]==1)
  if (yt_gen == 'correct'){
    pY0 <- mean(trainDat$Y[trainDat$T==0])
    transCoefs <- glm(trainDat$Y[keep] ~ trainDat$Y[keep-1] + trainDat$logHID[keep] + trainDat$T[keep] + log(trainDat$T[keep]+1) +
                        tDiff1 + trainDat$Y[keep-1]*tDiff1, family='quasibinomial')$coefficients
  }
  else if (yt_gen == 'complex'){
    pY0 <- mean(trainDat$Y[trainDat$T==0])
    transCoefs <- glm(trainDat$Y[keep] ~ trainDat$Y[keep-1] + trainDat$logHID[keep] + trainDat$T[keep] + log(trainDat$T[keep]+1) + tDiff1 +
                        trainDat$Y[keep-1]*trainDat$T[keep] + trainDat$Y[keep-1]*log(trainDat$T[keep]+1) + trainDat$Y[keep-1]*tDiff1, family='quasibinomial')$coefficients
  }
  else if (yt_gen == 'independent'){
    pY0 <- mean(trainDat$Y)
    transCoefs <- glm(trainDat$Y[keep] ~ trainDat$Y[keep-1] + trainDat$logHID[keep] + tDiff1 + trainDat$Y[keep-1]*tDiff1, family='quasibinomial')$coefficients
  }
  else{
    stop('yt_gen must be correct, complex, or independent')
  }
  
  
  # Select informative variables
  if (!is.na(nInformative)){
    cors <- sapply(Xinds,function(i){cor(trainDat[,i],trainDat$Y,method='spearman')})
    informative <- Xinds[order(cors,decreasing=TRUE)[1:nInformative]]
  }
  else{
    informative <- Xinds
  }
  
  
  # Fit I(X>0) | Y
  nonZeroModels <- sapply(Xinds,function(i){
    if (i %in% informative){
      mod <- glm(I(trainDat[,i]>0)~trainDat$Y+I(log(trainDat$T+1))+trainDat$logHID, family='quasibinomial')
      c(mod$coefficients,rstudent(mod))
    }
    else{
      mod <- glm(I(trainDat[,i]>0)~I(log(trainDat$T+1))+trainDat$logHID, family='quasibinomial')
      c(mod$coefficients[1],0,mod$coefficients[2:3],rstudent(mod))
    }
  })
  nonZeroCoefs <- t(nonZeroModels[1:4,])
  nonZeroResids <- nonZeroModels[-c(1:4),]
  nonZeroCovs <- 1/(nrow(nonZeroResids)-1) * t(nonZeroResids) %*% nonZeroResids + 1e-6*diag(ncol(nonZeroResids))
  sigs <- sqrt(diag(nonZeroCovs))
  nonZeroCors <- nonZeroCovs / (sigs %*% t(sigs))
  
  
  # Fit logX | X>0,Y
  logCountModels <- sapply(Xinds,function(i){
    keep <- which(trainDat[,i] > 0)
    if (i %in% informative){
      mod <- lm(trainDat[keep,i]~trainDat$Y[keep]+I(log(trainDat$T[keep]+1))+trainDat$logHID[keep])
      c(mod$coefficients[1]+1,mod$coefficients[2:4],summary(mod)$sigma)
    }
    else{
      mod <- lm(trainDat[keep,i]~I(log(trainDat$T[keep]+1))+trainDat$logHID[keep])
      c(mod$coefficients[1]+1,0,mod$coefficients[2:3],summary(mod)$sigma)
    }
  })
  logCountResids <- sapply(Xinds,function(i){
    mod <- lm(trainDat[,i]~trainDat$Y+I(log(trainDat$T+1))+trainDat$logHID)
    mod$residuals
  })
  logCountCoefs <- t(logCountModels[1:4,])
  logCountCovs <- 1/(nrow(logCountResids)-1) * t(logCountResids) %*% logCountResids + 1e-6*diag(ncol(logCountResids))
  sigs <- sqrt(diag(logCountCovs))
  logCountCovs <- logCountCovs / (sigs %*% t(sigs)) * (logCountModels[5,] %*% t(logCountModels[5,]))
  
  
  # Sample logHIDs and Ts
  IDs <- sample(unique(totDat$ID),nTrain+nTest,replace=TRUE)
  idx <- unlist(sapply(IDs,function(id){which(totDat$ID==id)}))
  
  simulated <- data.frame('logHID'=totDat$logHID[idx],'T'=totDat$T[idx])
  simulated$H <- exp(simulated$logHID) - 1
  simulated$ID <- unlist(sapply(1:length(IDs),function(i){rep(i,sum(totDat$ID==IDs[i]))}))
  
  
  # Simulate Y
  simulated$Y <- unlist(sapply(unique(simulated$ID),function(id){
    matches <- which(simulated$ID==id)
    ys = c(rbinom(1,1,pY0),rep(NA,length(matches)-1))
    if (length(matches) >= 2){
      for (i in 2:length(matches)){
        tDiff <- as.integer(simulated$T[matches[i]]-simulated$T[matches[i-1]])
        if (yt_gen == 'correct'){
          pY <- expit(transCoefs %*% c(1,ys[i-1],simulated$logHID[matches[i]],simulated$T[matches[i]],
                                       log(simulated$T[matches[i]]+1),tDiff,tDiff*ys[i-1]))
        }
        else if (yt_gen == 'complex'){
          pY <- expit(transCoefs %*% c(1,ys[i-1],simulated$logHID[matches[i]],simulated$T[matches[i]],
                                       log(simulated$T[matches[i]]+1),tDiff,simulated$T[matches[i]]*ys[i-1],
                                       log(simulated$T[matches[i]]+1)*ys[i-1],tDiff*ys[i-1]))
        }
        else if (yt_gen == 'independent'){
          pY <- expit(transCoefs %*% c(1,ys[i-1],simulated$logHID[matches[i]],tDiff,tDiff*ys[i-1]))
        }
        ys[i] <- rbinom(1,1,pY)
      }
    }
    ys
  }))
  
  
  # Simulate X
  covs <- rbind(1,simulated$Y,log(simulated$T+1),simulated$logHID)
  
  X <- t(matrix(unlist(sapply(unique(simulated$ID),function(id){
    if (id %% 1000 == 0){
      print(paste('Generating X',id))
    }
    keep <- which(simulated$ID==id)
    muNonZero <- expit(nonZeroCoefs %*% covs[,keep])
    muLogCount <- logCountCoefs %*% covs[,keep]
    xs <- array(dim=dim(muNonZero))
    ts <- simulated$T[keep]
    
    if (xy_gen == 'lognormal'){
      xs[,1] <- rmvbin(1,muNonZero[,1],sigma=nonZeroCors) *
        mixtools::rmvnorm(1,muLogCount[,1],sigma=logCountCovs)
    }
    else if (xy_gen == 'logt'){
      xs[,1] <- rmvbin(1,muNonZero[,1],sigma=nonZeroCors) *
        rmvt(1,sigma=logCountCovs,df=5,delta=muLogCount[,1])
    }
    else{
      stop('xy_gen must be either lognormal or logt')
    }
    xs[is.na(xs[,1]),1] <- 0
    
    if (ncol(muNonZero) > 1){
      for (i in 2:ncol(muNonZero)){
        rho_i <- rho * 0.5^as.integer(ts[i]-ts[i-1]>1)
        
        if (xy_gen == 'lognormal'){
          xs[,i] <- rmvbin(1,muNonZero[,i]+rho_i*((xs[,i-1]!=0)-muNonZero[,i-1]),sigma=nonZeroCors) *
            mixtools::rmvnorm(1,muLogCount[,i]+rho_i*(xs[,i-1]-muLogCount[,i-1]),sigma=logCountCovs)
        }
        else if (xy_gen == 'logt'){
          xs[,i] <- rmvbin(1,muNonZero[,i]+rho_i*((xs[,i-1]!=0)-muNonZero[,i-1]),sigma=nonZeroCors) *
            rmvt(1,sigma=logCountCovs,df=5,delta=muLogCount[,i]+rho_i*(xs[,i-1]-muLogCount[,i-1]))
        }
        xs[is.na(xs[,i]),i] <- 0
      }     
    }
    xs
  })),nrow=dim(nonZeroCoefs)[1]))
  
  
  # Produce train and test sets
  simulated <- cbind(simulated,X)
  trainIDs <- sample(unique(simulated$ID),nTrain)
  testIDs <- setdiff(unique(simulated$ID),trainIDs)
  train <- simulated[simulated$ID %in% trainIDs,]
  test <- simulated[simulated$ID %in% testIDs,]
  remove(simulated)
  
  
  # Simulate embeddings
  X_nonZero <- X>0
  remove(X)
  X_cooccur <- t(X_nonZero) %*% X_nonZero
  X_cooccur_colSums <- colSums(X_cooccur)
  X_PMI <- log(X_cooccur*sum(X_cooccur)/(X_cooccur_colSums %*% t(X_cooccur_colSums)))
  X_PMI[is.nan(X_PMI)] <- 0
  X_PMI[X_PMI < -10] <- -10
  embeddings <- eigen(X_PMI)$vectors[,1:embedDim]
  
  
  return(list('trainSet'=train,'testSet'=test,'embeddings'=embeddings))
}

set.seed(123)

nTrain <- 500
nTest <- 100
embedDim <- 10
load('~/inst/extdata/trainData.RData')

simdata <- simulate(trainData[[1]], trainData[[2]], 6:160, nTrain, nTest, embedDim)
str(simdata)

usethis::use_data(simdata, overwrite = TRUE)

