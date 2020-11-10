#setwd('~/Documents/HMS/Biostatistics PhD/Research/Cai Lab Summer 2019/EHR')

library(glmnet)
library(stats)
library(flexmix)
library(PRROC)
library(base)
library(e1071)
library(mixtools)
library(mvtnorm)
library(lbfgs)
library(nlme)
library(mixtools)
library(expm)
library(pROC)
library(bindata)
library(foreach)

# source("../Scripts/library_v2.R")  ###
# source("../Scripts/library_v3.r")  ###
# source("../Scripts/Library.R")       ###AUC\\\


logit <- function(x){log(x/(1-x))}
expit <- function(x){1/(1+exp(-x))}


fitGLS <- function(dat,nX=11,r=1){
  repInds <- 1 + which(dat$ID[-1]==dat$ID[-nrow(dat)])
  
  coefs <- lapply(1:nX,function(i){
    model <- gls(as.formula(paste0('X',i,' ~ Y+H+T+I(log(T+1))+Y:H+Y:T+Y:I(log(T+1))')), data=dat,
                 weights=varComb(varPower(form=~H), varFixed(~pInv)))
    list('coefs'=model$coefficients,'std.errors'=model$varBeta,'sigma'=model$sigma,'resids'=model$residuals,
         'varHID'=unlist(model$modelStruct$varStruct),'varVar'=model$apVar)
  })
  
  coefs$beta_hat <- t(sapply(1:nX,function(i){coefs[[i]]$coefs}))
  coefs$sigma <- sapply(1:nX,function(i){coefs[[i]]$sigma}) * sqrt((nrow(dat)-6)/(sum(dat$pY)-6))
  coefs$varHID <- sapply(1:nX,function(i){unlist(coefs[[i]]$varHID)})
  
  residuals_normalized <- sapply(1:nX,function(i){
    coefs[[i]]$resids / (dat$HID^coefs$varHID[i])
  })
  covMat <- 1/(sum(dat$pY)-6) * (t(residuals_normalized*dat$pY) %*% residuals_normalized)
  sigma <- sqrt(diag(covMat))
  coefs$corrMat <- covMat / (sigma %*% t(sigma))
  
  # Vector auto-regression (1 lag)
  coefs$autoCoefs <- r * sapply(1:nX,function(i){
    rep(lm(residuals_normalized[repInds,i]~residuals_normalized[repInds-1,i]:as.factor(dat$T[repInds]-dat$T[repInds-1]==1),
              weights=(dat$pY[repInds]*dat$pY[repInds-1]))$coef[-1],2)
  })
  
  return(coefs)
}


trainTransitionMatrix <- function(train){
  counts <- matrix(rep(0,4),2,2)
  rownames(counts) <- c('Start0','Start1'); colnames(counts) <- c('End0','End1')
  
  for (i in 2:nrow(train)){
    if (train$ID[i] == train$ID[i-1]){
      counts <- counts + train$pY[i] * outer(c(1-train$Y[i-1],train$Y[i-1]), c(1-train$Y[i],train$Y[i]))
    }
  }
  
  transition <- counts / rowSums(counts)
  return(list('counts'=counts,'tmat'=transition))
}


Estep_partial <- function(dat,trained,nX=11){
  priorModel <- trained$priorModel; likeModel <- trained$likeModel; transCoefs <- trained$transCoefs
  #  a <- likeModel$autoCoefs
  #  amin <- apply(a,2,min)
  a1 <- likeModel$autoCoefs[3:4,]; a2 <- likeModel$autoCoefs[1:2,]
  amin1 <- apply(a1,2,min); amin2 <- apply(a2,2,min)
  prior_fitted <- predict(priorModel,dat,type='response')
  
  post <- unlist(sapply(unique(dat$ID), function(patient){
    encounters <- which(dat$ID == patient)
    Ti <- dat$T[encounters]; DDi <- dat$DD[encounters]; Hi <- mean(dat$H[encounters])
    prior_pat <- prior_fitted[encounters]
    
    #    print(likeModel$beta_hat)
    #    print(rbind(1,0,Hi,DDi,0,0))
    
    mu0 <- likeModel$beta_hat %*% rbind(1,0,Hi,Ti,DDi,0,0,0)
    mu1 <- likeModel$beta_hat %*% rbind(1,1,Hi,Ti,DDi,Hi,Ti,DDi)
    mus <- list(mu0,mu1)
    
    sigsq_base <- (likeModel$sigma * unique(dat$HID[encounters])^likeModel$varHID)^2
    sigsq_prior <- sapply(1:nX,function(i){c(1,1,Hi,mean(Ti),mean(DDi),Hi,mean(Ti),mean(DDi)) %*% likeModel[[i]]$std.errors %*% c(1,1,Hi,mean(Ti),mean(DDi),Hi,mean(Ti),mean(DDi))})
    sig <- sqrt(sigsq_base + sigsq_prior)
    Sigma <- likeModel$corrMat * (sig %*% t(sig))
    A1 <- sqrt(1-amin1^2) %*% t(sqrt(1-amin1^2)); A2 <- sqrt(1-amin2^2) %*% t(sqrt(1-amin2^2))
    Sigma_cond1 <- A1 * Sigma; Sigma_cond2 <- A2 * Sigma
    #    Sigma_cond1 <- Sigma; Sigma_cond2 <- Sigma
    #    A <- sqrt(1-amin^2) %*% t(sqrt(1-amin^2))
    #    Sigma_cond <- A * Sigma
    
    Xi <- as.matrix(dat[encounters,paste0('X',1:nX)])
    
    sapply(1:length(Ti),function(i){
      if (length(Ti)==1){
        logprior <- log(c(1-prior_pat[i],prior_pat[i]))
        logprobs <- logprior + c(logdmvnorm(Xi[i,],mus[[1]][,i],Sigma), logdmvnorm(Xi[i,],mus[[2]][,i],Sigma))
        probs <- exp(logprobs - max(logprobs))
        c(probs[2]/sum(probs),rep(0,4))
      }
      else if (i==1){
        X_comp <- c(Xi[i,],Xi[i+1,])
        prior_i <- c(1-prior_pat[i],prior_pat[i])
        if (Ti[i+1]-Ti[i]==1){
          a <- a1; Sigma_cond <- Sigma_cond1
        }
        else{
          a <- a2; Sigma_cond <- Sigma_cond2
        }
        
        logprobs <- sapply(1:2,function(yn){sapply(1:2,function(yc){
#          p <- expit(transCoefs %*% c(1,yc-1,Hi,DDi[i+1],as.integer(Ti[i+1]-Ti[i]==1),as.integer(Ti[i+1]-Ti[i]==1)*(yc-1)))
          p <- expit(transCoefs %*% c(1,yc-1,Hi,Ti[i+1],log(Ti[i+1]+1),as.integer(Ti[i+1]-Ti[i]==1),as.integer(Ti[i+1]-Ti[i]==1)*(yc-1)))
          logprior <- log(prior_i[yc] * c(1-p,p)[yn])
          muC <- mus[[yc]][,i]; muN <- mus[[yn]][,i+1] + a[1+(yc==yn),]*(Xi[i,]-muC)
          loglike <- logdmvnorm(Xi[i,],muC,Sigma) + logdmvnorm(Xi[i+1,],muN,Sigma_cond)
          logprior + loglike
        })})
        probs <- exp(logprobs - max(logprobs))
        c(sum(probs[2,])/sum(probs),rep(0,4))
      }
      else if (i==length(Ti)){
        X_comp <- c(Xi[i-1,],Xi[i,])
        prior_i <- c(1-prior_pat[i-1],prior_pat[i-1])
        if (Ti[i]-Ti[i-1]==1){
          a <- a1; Sigma_cond <- Sigma_cond1
        }
        else{
          a <- a2; Sigma_cond <- Sigma_cond2
        }
        
        logprobs <- sapply(1:2,function(yn){sapply(1:2,function(yc){
#          p <- expit(transCoefs %*% c(1,yc-1,Hi,DDi[i],as.integer(Ti[i]-Ti[i-1]==1),as.integer(Ti[i]-Ti[i-1]==1)*(yc-1)))
          p <- expit(transCoefs %*% c(1,yc-1,Hi,Ti[i],log(Ti[i]+1),as.integer(Ti[i]-Ti[i-1]==1),as.integer(Ti[i]-Ti[i-1]==1)*(yc-1)))
          logprior <- log(prior_i[yc] * c(1-p,p)[yn])
          muC <- mus[[yc]][,i-1]; muN <- mus[[yn]][,i] + a[1+(yc==yn),]*(Xi[i-1,]-muC)
          loglike <- logdmvnorm(Xi[i-1,],muC,Sigma) + logdmvnorm(Xi[i,],muN,Sigma_cond)
          logprior + loglike
        })})
        probs <- exp(logprobs - max(logprobs))
        c(sum(probs[,2])/sum(probs), probs/sum(probs))
      }
      else{
        X_comp <- c(Xi[i-1,],Xi[i,],Xi[i+1,])
        prior_i <- c(1-prior_pat[i-1],prior_pat[i-1])
        if (Ti[i]-Ti[i-1]==1){
          aL <- a1; Sigma_condL <- Sigma_cond1
        }
        else{
          aL <- a2; Sigma_condL <- Sigma_cond2
        }
        if (Ti[i+1]-Ti[i]==1){
          aN <- a1; Sigma_condN <- Sigma_cond1
        }
        else{
          aN <- a2; Sigma_condN <- Sigma_cond2
        }
        
        logprobs <- array(sapply(1:2,function(yn){sapply(1:2,function(yc){sapply(1:2,function(yl){
#          p1 <- expit(transCoefs %*% c(1,yl-1,Hi,DDi[i],as.integer(Ti[i]-Ti[i-1]==1),as.integer(Ti[i]-Ti[i-1]==1)*(yl-1)))
#          p2 <- expit(transCoefs %*% c(1,yc-1,Hi,DDi[i+1],as.integer(Ti[i+1]-Ti[i]==1),as.integer(Ti[i+1]-Ti[i]==1)*(yc-1)))
          p1 <- expit(transCoefs %*% c(1,yl-1,Hi,Ti[i],log(Ti[i]+1),as.integer(Ti[i]-Ti[i-1]==1),as.integer(Ti[i]-Ti[i-1]==1)*(yl-1)))
          p2 <- expit(transCoefs %*% c(1,yc-1,Hi,Ti[i+1],log(Ti[i+1]+1),as.integer(Ti[i+1]-Ti[i]==1),as.integer(Ti[i+1]-Ti[i]==1)*(yc-1)))
          logprior <- log(prior_i[yl] * c(1-p1,p1)[yc] * c(1-p2,p2)[yn])
          muL <- mus[[yl]][,i-1]; muC <- mus[[yc]][,i] + aL[1+(yl==yc),]*(Xi[i-1,]-muL); muN <- mus[[yn]][,i+1] + aN[1+(yc==yn),]*(Xi[i,]-muC)
          loglike <- logdmvnorm(Xi[i-1,],muL,Sigma) + logdmvnorm(Xi[i,],muC,Sigma_condL) + logdmvnorm(Xi[i+1,],muN,Sigma_condN)
          logprior + loglike
        })})}), dim=c(2,2,2))
        probs <- exp(logprobs - max(logprobs))
        c(sum(probs[,2,])/sum(probs), probs[,,1]/sum(probs[,,1]))
      }
    })
  }))
  output <- post[seq(1,length(post),5)]
  attr(output,'post2') <- array(post[-seq(1,length(post),5)],dim=c(2,2,length(output)))
  output
}


Mstep <- function(train,train_interTemp=NULL,nX=11,r=1){
  priorModel <- glm(Y~H+T+I(log(T+1)),weights=pY,data=train,family='quasibinomial')
  if (is.null(train_interTemp)){
    repInds <- 1 + which(train$ID[-1]==train$ID[-nrow(train)])
    transCoefs <- glm(train$Y[repInds]~train$Y[repInds-1]+train$H[repInds]+train$T[repInds]+log(train$T[repInds]+1)+
                        as.factor(train$T[repInds]-train$T[repInds-1]==1)+as.factor(train$T[repInds]-train$T[repInds-1]==1):train$Y[repInds-1],family='quasibinomial')$coefficients
  }
  else{
    transCoefs <- glm(Ycurr~Yprev+H+T+I(log(T+1))+as.factor(T-Tprev==1)+Yprev:as.factor(T-Tprev==1),
                      weights=pY,data=train_interTemp,family='quasibinomial')$coefficients
  }
  likeModel <- fitGLS(train,nX=nX,r=r)
  return(list('priorModel'=priorModel, 'transCoefs'=transCoefs, 'likeModel'=likeModel))
}


# Full Markov implementation
Estep_full <- function(dat,trained,nX=11){
  priorModel <- trained$priorModel; likeModel <- trained$likeModel; transCoefs <- trained$transCoefs
  #  a <- likeModel$autoCoefs
  #  amin <- apply(a,2,min)
  a1 <- likeModel$autoCoefs[3:4,]; a2 <- likeModel$autoCoefs[1:2,]
  amin1 <- apply(a1,2,min); amin2 <- apply(a2,2,min)
  prior_fitted <- predict(priorModel,dat,type='response')
  
  post1 <- unlist(sapply(unique(dat$ID), function(patient){
    encounters <- which(dat$ID == patient)
    Ti <- dat$T[encounters]
    DDi <- dat$DD[encounters]
    Hi <- unique(dat$H[encounters])
    prior_pat <- prior_fitted[encounters]
    
    mu0 <- likeModel$beta_hat %*% rbind(1,0,Hi,Ti,DDi,0,0,0)
    mu1 <- likeModel$beta_hat %*% rbind(1,1,Hi,Ti,DDi,Hi,Ti,DDi)
    mus <- list(mu0,mu1)
    
    sigsq_base <- (likeModel$sigma * unique(dat$HID[encounters])^likeModel$varHID)^2
    sigsq_prior <- sapply(1:nX,function(i){c(1,1,Hi,mean(Ti),mean(DDi),Hi,mean(Ti),mean(DDi)) %*% likeModel[[i]]$std.errors %*% c(1,1,Hi,mean(Ti),mean(DDi),Hi,mean(Ti),mean(DDi))})
    sig <- sqrt(sigsq_base + sigsq_prior)
    Sigma <- likeModel$corrMat * (sig %*% t(sig))
    A1 <- sqrt(1-amin1^2) %*% t(sqrt(1-amin1^2)); A2 <- sqrt(1-amin2^2) %*% t(sqrt(1-amin2^2))
    Sigma_cond1 <- A1 * Sigma; Sigma_cond2 <- A2 * Sigma
    #    A <- sqrt(1-amin^2) %*% t(sqrt(1-amin^2))
    #    Sigma_cond <- A * Sigma
    
    Xi <- as.matrix(dat[encounters,paste0('X',1:nX)])
    
    fwd <- matrix(0,2,length(encounters))
    f_prev <- c(1-prior_pat[1],prior_pat[1])
    for (i in 1:length(encounters)){
      if (i==1){
        fwd[,i] <- f_prev * c(dmvnorm(Xi[i,],mu0[,i],Sigma),
                              dmvnorm(Xi[i,],mu1[,i],Sigma))
      }
      else{
        if (Ti[i]-Ti[i-1]==1){
          a <- a1; Sigma_cond <- Sigma_cond1
        }
        else{
          a <- a2; Sigma_cond <- Sigma_cond2
        }
#        p <- expit(cbind(1,c(0,1),Hi,DDi[i],as.integer(Ti[i]-Ti[i-1]==1),as.integer(Ti[i]-Ti[i-1]==1)*c(0,1)) %*% transCoefs)
        p <- expit(cbind(1,c(0,1),Hi,Ti[i],log(Ti[i]+1),as.integer(Ti[i]-Ti[i-1]==1),as.integer(Ti[i]-Ti[i-1]==1)*c(0,1)) %*% transCoefs)
        transition <- cbind(1-p,p)
        probs <- (t(transition) * rbind(f_prev,f_prev))
        loglike <- matrix(c(logdmvnorm(Xi[i,],mu0[,i]+a[2,]*(Xi[i-1,]-mu0[,i-1]),Sigma_cond),
                            logdmvnorm(Xi[i,],mu1[,i]+a[1,]*(Xi[i-1,]-mu0[,i-1]),Sigma_cond),
                            logdmvnorm(Xi[i,],mu0[,i]+a[1,]*(Xi[i-1,]-mu1[,i-1]),Sigma_cond),
                            logdmvnorm(Xi[i,],mu1[,i]+a[2,]*(Xi[i-1,]-mu1[,i-1]),Sigma_cond)),2,2)
        probs <- probs * exp(loglike - max(loglike))
        fwd[,i] <- probs %*% c(1,1)
      }
      
      fwd[,i] <- fwd[,i]/sum(fwd[,i])
      f_prev <- fwd[,i]
    }
    
    bkw <- matrix(0,2,length(encounters))
    b_next <- bkw[,length(encounters)] <- c(1,1)
    if (length(encounters) > 1){
      for (i in (length(encounters)-1):1){
        if (Ti[i+1]-Ti[i]==1){
          a <- a1; Sigma_cond <- Sigma_cond1
        }
        else{
          a <- a2; Sigma_cond <- Sigma_cond2
        }
#        p <- expit(cbind(1,c(0,1),Hi,DDi[i+1],as.integer(Ti[i+1]-Ti[i]==1),as.integer(Ti[i+1]-Ti[i]==1)*c(0,1)) %*% transCoefs)
        p <- expit(cbind(1,c(0,1),Hi,Ti[i+1],log(Ti[i+1]+1),as.integer(Ti[i+1]-Ti[i]==1),as.integer(Ti[i+1]-Ti[i]==1)*c(0,1)) %*% transCoefs)
        transition <- cbind(1-p,p)
        
        probs <- (transition * rbind(b_next,b_next))
        loglike <- matrix(c(logdmvnorm(Xi[i+1,],mu0[,i+1]+a[2,]*(Xi[i,]-mu0[,i]),Sigma_cond),
                            logdmvnorm(Xi[i+1,],mu1[,i+1]+a[1,]*(Xi[i,]-mu0[,i]),Sigma_cond),
                            logdmvnorm(Xi[i+1,],mu0[,i+1]+a[1,]*(Xi[i,]-mu1[,i]),Sigma_cond),
                            logdmvnorm(Xi[i+1,],mu1[,i+1]+a[2,]*(Xi[i,]-mu1[,i]),Sigma_cond)),2,2)
        probs <- probs * exp(loglike - max(loglike))
        
        bkw[,i] <- probs %*% c(1,1)
        bkw[,i] <- bkw[,i]/sum(bkw[,i])
        b_next <- bkw[,i]
      }
    }
    
    post <- fwd[2,]*bkw[2,]/colSums(fwd*bkw)
    post[colSums(fwd*bkw)==0] <- prior_pat[colSums(fwd*bkw)==0]
    post
  }))
  attr(post1,'post2') <- array(cbind(0,sapply(2:length(post1),function(i){
    c(1-post1[i-1],post1[i-1]) %*% t(c(1-post1[i],post1[i]))
  })),dim=c(2,2,length(post1)))
  post1
}


cumulative_prob <- function(dat,trained,nX=11){
  priorModel <- trained$priorModel; likeModel <- trained$likeModel; transCoefs <- trained$transCoefs
  #  a <- likeModel$autoCoefs
  #  amin <- apply(a,2,min)
  a1 <- likeModel$autoCoefs[3:4,]; a2 <- likeModel$autoCoefs[1:2,]
  amin1 <- apply(a1,2,min); amin2 <- apply(a2,2,min)
  prior_fitted <- predict(priorModel,dat,type='response')
  
  cumprob <- unlist(sapply(unique(dat$ID), function(patient){
    encounters <- which(dat$ID == patient)
    Ti <- dat$T[encounters]
    DDi <- dat$DD[encounters]
    Hi <- unique(dat$H[encounters])
    prior_pat <- prior_fitted[encounters]
    
    mu0 <- likeModel$beta_hat %*% rbind(1,0,Hi,Ti,DDi,0,0,0)
    mu1 <- likeModel$beta_hat %*% rbind(1,1,Hi,Ti,DDi,Hi,Ti,DDi)
    mus <- list(mu0,mu1)
    
    sigsq_base <- (likeModel$sigma * unique(dat$HID[encounters])^likeModel$varHID)^2
    sigsq_prior <- sapply(1:nX,function(i){c(1,1,Hi,mean(Ti),mean(DDi),Hi,mean(Ti),mean(DDi)) %*% likeModel[[i]]$std.errors %*% c(1,1,Hi,mean(Ti),mean(DDi),Hi,mean(Ti),mean(DDi))})
    sig <- sqrt(sigsq_base + sigsq_prior)
    Sigma <- likeModel$corrMat * (sig %*% t(sig))
    A1 <- sqrt(1-amin1^2) %*% t(sqrt(1-amin1^2)); A2 <- sqrt(1-amin2^2) %*% t(sqrt(1-amin2^2))
    Sigma_cond1 <- A1 * Sigma; Sigma_cond2 <- A2 * Sigma
    #    A <- sqrt(1-amin^2) %*% t(sqrt(1-amin^2))
    #    Sigma_cond <- A * Sigma
    
    Xi <- as.matrix(dat[encounters,paste0('X',1:nX)])
    
    fwd <- matrix(0,2,length(encounters))
    f_prev <- c(1-prior_pat[1],prior_pat[1])
    survival <- rep(0,length(encounters))
    for (i in 1:length(encounters)){
      if (i==1){
        fwd[,i] <- f_prev * c(dmvnorm(Xi[i,],mu0[,i],Sigma),
                              dmvnorm(Xi[i,],mu1[,i],Sigma))
        survival[i] <- fwd[1,i]
      }
      else{
        if (Ti[i]-Ti[i-1]==1){
          a <- a1; Sigma_cond <- Sigma_cond1
        }
        else{
          a <- a2; Sigma_cond <- Sigma_cond2
        }
#        p <- expit(cbind(1,c(0,1),Hi,DDi[i],as.integer(Ti[i]-Ti[i-1]==1),as.integer(Ti[i]-Ti[i-1]==1)*c(0,1)) %*% transCoefs)
        p <- expit(cbind(1,c(0,1),Hi,Ti[i],log(Ti[i]+1),as.integer(Ti[i]-Ti[i-1]==1),as.integer(Ti[i]-Ti[i-1]==1)*c(0,1)) %*% transCoefs)
        transition <- cbind(1-p,p)
        probs <- (t(transition) * rbind(f_prev,f_prev))
        loglike <- matrix(c(logdmvnorm(Xi[i,],mu0[,i]+a[2,]*(Xi[i-1,]-mu0[,i-1]),Sigma_cond),
                            logdmvnorm(Xi[i,],mu1[,i]+a[1,]*(Xi[i-1,]-mu0[,i-1]),Sigma_cond),
                            logdmvnorm(Xi[i,],mu0[,i]+a[1,]*(Xi[i-1,]-mu1[,i-1]),Sigma_cond),
                            logdmvnorm(Xi[i,],mu1[,i]+a[2,]*(Xi[i-1,]-mu1[,i-1]),Sigma_cond)),2,2)
        like <- exp(loglike - max(loglike))
        probs <- probs * like
        fwd[,i] <- probs %*% c(1,1)
        survival[i] <- survival[i-1] * transition[1,1] * like[1,1]
      }
      
      survival[i] <- survival[i]/sum(fwd[,i])
      fwd[,i] <- fwd[,i]/sum(fwd[,i])
      f_prev <- fwd[,i]
    }
    
    bkw <- matrix(0,2,length(encounters))
    b_next <- bkw[,length(encounters)] <- c(1,1)
    if (length(encounters) > 1){
      for (i in (length(encounters)-1):1){
        if (Ti[i+1]-Ti[i]==1){
          a <- a1; Sigma_cond <- Sigma_cond1
        }
        else{
          a <- a2; Sigma_cond <- Sigma_cond2
        }
#        p <- expit(cbind(1,c(0,1),Hi,DDi[i+1],as.integer(Ti[i+1]-Ti[i]==1),as.integer(Ti[i+1]-Ti[i]==1)*c(0,1)) %*% transCoefs)
        p <- expit(cbind(1,c(0,1),Hi,Ti[i],log(Ti[i]+1),as.integer(Ti[i+1]-Ti[i]==1),as.integer(Ti[i+1]-Ti[i]==1)*c(0,1)) %*% transCoefs)
        transition <- cbind(1-p,p)
        
        probs <- (transition * rbind(b_next,b_next))
        loglike <- matrix(c(logdmvnorm(Xi[i+1,],mu0[,i+1]+a[2,]*(Xi[i,]-mu0[,i]),Sigma_cond),
                            logdmvnorm(Xi[i+1,],mu1[,i+1]+a[1,]*(Xi[i,]-mu0[,i]),Sigma_cond),
                            logdmvnorm(Xi[i+1,],mu0[,i+1]+a[1,]*(Xi[i,]-mu1[,i]),Sigma_cond),
                            logdmvnorm(Xi[i+1,],mu1[,i+1]+a[2,]*(Xi[i,]-mu1[,i]),Sigma_cond)),2,2)
        probs <- probs * exp(loglike - max(loglike))
        
        bkw[,i] <- probs %*% c(1,1)
        bkw[,i] <- bkw[,i]/sum(bkw[,i])
        b_next <- bkw[,i]
      }
    }
    
    surv <- survival*bkw[1,]/colSums(fwd*bkw)
    1-surv
  }))
  
  cumprob
}


EM <- function(train,observedPats,test=NULL,maxIt=1,r=0.8,tol=0.01,Estep=Estep_partial){
  print('Starting EM')
  
  observedIndices <- which(train$ID %in% observedPats)
  unobservedIndices <- setdiff(seq(nrow(train)),observedIndices)
  train$pY <- train$pWeights <- 1
  prediction <- aucs <- NULL
  lastY <- train$Y[unobservedIndices]
  
  trained_sup <- trained_semisup <- Mstep(train[observedIndices,],r=r)
  if (!is.null(test)){
    prediction <- Estep(test,trained_sup)
    print(paste('No E AUC =',auc(test$Y,prediction))) 
    aucs <- c(auc(test$Y,prediction),rep(0,maxIt))
  }
  
  for (it in 1:maxIt){
    
    # E-step
    prediction <- Estep(train[unobservedIndices,],trained_semisup)
    transProbs <- attr(prediction,'post2')
    
    train_augmented <- rbind(train[observedIndices,],train[unobservedIndices,],train[unobservedIndices,])
    train_augmented$pY <- c(rep(1,length(observedIndices)),1-prediction,prediction)
    train_augmented$pWeights <- 1/train_augmented$pY
    train_augmented$Y <- c(train$Y[observedIndices],rep(0,length(unobservedIndices)),rep(1,length(unobservedIndices)))
    train_augmented$ID <- c(train$ID[observedIndices],train$ID[unobservedIndices],paste0(train$ID[unobservedIndices],'.2'))
    train_augmented <- train_augmented[train_augmented$pWeights<1e8,]
    
    train_interTemp <- rbind(train[observedIndices[-1],],train[unobservedIndices[-1],],train[unobservedIndices[-1],],
                             train[unobservedIndices[-1],],train[unobservedIndices[-1],])
    train_interTemp$Ycurr <- c(train$Y[observedIndices[-1]],rep(0,length(unobservedIndices)-1),rep(0,length(unobservedIndices)-1),
                               rep(1,length(unobservedIndices)-1),rep(1,length(unobservedIndices)-1))
    train_interTemp$Yprev <- c(train$Y[observedIndices[-length(observedIndices)]],rep(0,length(unobservedIndices)-1),
                               rep(1,length(unobservedIndices)-1),rep(0,length(unobservedIndices)-1),rep(1,length(unobservedIndices)-1))
    train_interTemp$pY <- c(as.integer(train$ID[observedIndices[-1]]==train$ID[observedIndices[-length(observedIndices)]]),
                            transProbs[1,1,-1],transProbs[2,1,-1],transProbs[1,2,-1],transProbs[2,2,-1])
    train_interTemp$Tprev <- c(train$T[observedIndices[-length(observedIndices)]],
                               rep(train$T[unobservedIndices[-length(unobservedIndices)]],4))
    
    # M-step    
    trained_semisup <- Mstep(train_augmented,train_interTemp,r=r)
    
    if (!is.null(test)){
      testpred <- Estep(test,trained_semisup)
      print(paste('AUC Iteration',it,'=',auc(test$Y,testpred)))
      aucs[it+1] <- auc(test$Y,testpred)      
    }
    
    if (all(abs(prediction-lastY) < tol)){
      break
    }
    lastY <- prediction
  }
  
  return(list('fitted_semisup'=trained_semisup, 'fitted_sup'=trained_sup, 'trainComplete'=train,
              'trainUnobserved'=train[-observedIndices,], 'testPrediction'=prediction, 'aucs'=aucs))
}


Simulate <- function(C,V,trueDat,Estep=Estep_partial){
  d <- ncol(C)
  ratios <- sapply(1:ncol(C),function(i){
    (mean(C[trueDat$Y==1,i])-mean(C[trueDat$Y==0,i]))/sd(C[,i])
  })
  W <- sign(ratios) * (abs(ratios)>0.25) * rnorm(d,0,.25)
  
  X <- as.matrix(C) %*% (W*as.matrix(V))
  
  dat <- as.data.frame(cbind(trueDat[,c('ID','Y','T','DD','HID','H','X11')],X))
  colnames(dat) <- c('ID','Y','T','DD','HID','H','X11',paste0('X',1:10))
  dat$pY <- dat$pWeights <- 1
  
  margprobs  <- Estep(dat,Mstep(dat,r=0.5))
  Ysim <- unlist(sapply(unique(dat$ID),function(id){
    inds <- which(dat$ID==id)
    margprob <- margprobs[inds]
    sigma <- sapply(1:length(inds),function(i){
      sapply(1:length(inds),function(j){
        ifelse(i==j,1,0.8^abs(i-j))
      })})
    rmvbin(1,margprob,sigma=sigma)
  }))
  dat$Y <- Ysim
  
  return(list('data'=dat,'w'=W))
}


Simulate_noW <- function(N,Tn,trained,nX=11,dist='normal'){
  priorModel <- trained$priorModel; likeModel <- trained$likeModel; transCoefs <- trained$transCoefs
  Tns <- rpois(N,Tn)
  
  result <- foreach(i=1:N, .combine=rbind) %do% {
    if (i%%100==0){print(i)}
    HID <- 1 + rchisq(1,2.2)
    H <- log(HID+1)
    
    DD <- c(0,rep(NA,Tns[i]-1))
    T <- c(1,rep(NA,Tns[i]-1))
    for (j in 2:Tns[i]){
      delta <- rexp(1,3)
      DD[j] <- DD[j-1] + delta
      T[j] <- T[j-1] + ceiling(2*delta)
    }
    
    Yn <- c(rbinom(1,1,expit(c(1,0,H)%*%priorModel$coefficients)),rep(NA,Tns[i]-1))
    for (j in 2:Tns[i]){
      Yn[j] <- rbinom(1,1,expit(transCoefs%*%c(1,Yn[j-1],H,DD[j],T[j]-T[j-1],Yn[j-1]*(T[j]-T[j-1]))))
    }
    
    Hn <- as.matrix(cbind(1,Yn,H,DD,Yn*H,Yn*DD))
    mu <- c(likeModel$beta_hat %*% t(Hn))
    
    sig <- likeModel$sigma * HID^likeModel$varHID
    Sig0 <- likeModel$corrMat * (sig %*% t(sig))
    Amean <- colMeans(likeModel$autoCoefs[3:4,])
    A <- sqrt(Amean) %*% t(sqrt(Amean))
    Sigma <- matrix(0,Tns[i]*nX,Tns[i]*nX)
    for (j in 1:Tns[i]){
      for (k in 1:Tns[i]){
        Sigma[(1+(j-1)*nX):(j*nX), (1+(k-1)*nX):(k*nX)] <- Sig0 * A^abs(T[j]-T[k])
      }
    }
    
    if (dist == 'normal'){
      Xn <- t(matrix(rmvnorm(1,mu,Sigma),nrow=nX,ncol=Tns[i]))
    }
    else if (dist == 't'){
      Xn <- t(matrix(rmvt(1,delta=mu,sigma=Sigma,df=5),nrow=nX,ncol=Tns[i]))
    }
    else{
      stop(error='dist must be either normal or t')
    }
    
    cbind(i,HID,H,T,log(DD+1),Yn,Xn,1,1)
  }
  
  colnames(result) <- c('ID','HID','H','T','DD','Y',paste0('X',seq(nX)),'pWeights','pY')
  return(as.data.frame(result))
}


lineSearch <- function(train,observedPats,test=NULL,nCrosses=5,alphas=seq(0,1,.1),r=0.8,Estep=Estep_partial){
  print('Starting line search')
  
  if (length(alphas) == 1){
    alpha <- alphas
  }
  else{
    n <- length(observedPats)
    observedPats <- sample(observedPats)
    valPats_overall <- lapply(1:nCrosses, function(i){observedPats[round(n*(i-1)/nCrosses+1):round(n*i/nCrosses)]})
    
    alpha_results <- as.matrix(sapply(1:nCrosses,function(i){
      print(paste('Starting cross',i))
      
      validatePats <- valPats_overall[[i]]
      trainPats <- setdiff(observedPats,validatePats)
      validateIndices <- which(train$ID %in% validatePats)
      
      tryCatch({
        em <- EM(train,trainPats,r=r,Estep=Estep)
        supervised <- Estep(train[validateIndices,],em$fitted_sup)
        semisupervised <- Estep(train[validateIndices,],em$fitted_semisup)
        
        sapply(alphas,function(alpha){
          mixture <- alpha*semisupervised + (1-alpha)*supervised
          auc(train$Y[validateIndices],mixture)
        })
      },
      error=function(e){
        rep(NA,length(alphas))
      })
    }))
    alpha <- alphas[which.max(rowMeans(alpha_results,na.rm=T))]
  }
  
  if (!is.null(test)){
    em <- EM(train,observedPats,test,r=r,Estep=Estep)
    supervised <- Estep(test,em$fitted_sup)
    semisupervised <- Estep(test,em$fitted_semisup)
    mixture <- alpha*semisupervised + (1-alpha)*supervised
    resultSup <- evaluate(test$Y,supervised)
    resultSemisup <- evaluate(test$Y,semisupervised)
    resultMix <- evaluate(test$Y,mixture)
    
    cumSup <- cumulative_prob(test,em$fitted_sup)
    cumSemisup <- cumulative_prob(test,em$fitted_semisup)
    cumMixture <- alpha*cumSemisup + (1-alpha)*cumSup
  }
  else{
    supervised <- semisupervised <- mixture <- resultSup <- resultSemisup <- resultMix <- cumSup <- cumSemisup <- cumMixture <- NULL
  }
  
  return(list('alpha'=alpha,'prediction'=mixture,
              'margSup'=supervised,'margSemisup'=semisupervised,'margMix'=mixture,
              'resultSup'=resultSup,'resultSemisup'=resultSemisup,'resultMix'=resultMix,
              'cumSup'=cumSup,'cumSemisup'=cumSemisup,'cumMix'=cumMixture))
}


cv.r <- function(train,observedPats,nCrosses=5,rs=seq(0,1,.1),Estep=Estep_partial){
  print('Starting cross validation of r')
  n <- length(observedPats)
  observedPats <- sample(observedPats)
  valPats_overall <- lapply(1:nCrosses, function(i){observedPats[round(n*(i-1)/nCrosses+1):round(n*i/nCrosses)]})
  
  grid <- sapply(1:nCrosses,function(i){
    print(paste('Starting cross',i))
    
    validatePats <- valPats_overall[[i]]
    trainPats <- setdiff(observedPats,validatePats)
    
    sapply(rs,function(r){
      fitted_M <- Mstep(train[train$ID %in% trainPats,],r=r)
      supervised <- Estep(train[train$ID %in% validatePats,],fitted_M)
      print(paste(r,auc(train$Y[train$ID %in% validatePats],supervised)))
      auc(train$Y[train$ID %in% validatePats],supervised)
    })
  })
  means <- rowMeans(grid)
  
  return(list('results'=grid,'r_opt'=rs[which.max(means)]))
}


objective_w <- function(w,args,lambda=0){
  C <- args$C; y <- args$y; V <- args$V
  X <- as.matrix(C) %*% (w*V)
  X0 <- X[y==0,]
  X1 <- X[y==1,]
  mu0 <- colMeans(X0)
  mu1 <- colMeans(X1)
  epsilon0 <- t(X0) - mu0
  epsilon1 <- t(X1) - mu1
  Sigma <- 1/(nrow(X)-1) * (epsilon0%*%t(epsilon0) + epsilon1%*%t(epsilon1) + 1e-6*diag(length(mu0)))
  
  return(c(-1 * (mu1 - mu0) %*% solve(Sigma) %*% (mu1 - mu0) + lambda * sum(abs(w))))
}


numericGradientDescent <- function(x0, f, args=NULL, constIndex=1, alphas=c(1e-4,2e-4,5e-4,1e-3,2e-3,5e-3,.01,.02,.05,.1,.2,.5,1),
                                   lambda=0, maxIt=100, tol=1e-4){
  for (it in 1:maxIt){
    grad <- nl.grad(x0,f,heps = .Machine$double.eps^(1/3),args,lambda)
    if (any(x0==0)){grad[x0==0] <- sign(grad[x0==0]) * sapply(abs(grad[x0==0])-lambda, function(xi){max(xi,0)})}
    grad[constIndex] <- 0
    grad <- grad/c(sqrt(grad%*%grad))
    
    alphaDeltas <- sapply(alphas,function(alpha){
      deltas <- alpha*grad*sqrt(length(x0))
      x1 <- sapply(1:length(x0),function(i){
        if (x0[i]>0){max(x0[i]-deltas[i],0)}
        else if (x0[i]<0){min(x0[i]-deltas[i],0)}
        else{x0[i]-deltas[i]}
      })
      f(x1, args, lambda)
    })
    
    if (all(alphaDeltas >= c(f(x0,args,lambda)-tol))){
      break
    }
    else{
      alpha <- alphas[which.min(alphaDeltas)]
      deltas <- alpha*grad*sqrt(length(x0))
      x0 <- sapply(1:length(x0),function(i){
        if (x0[i]>0){max(x0[i]-deltas[i],0)}
        else if (x0[i]<0){min(x0[i]-deltas[i],0)}
        else{x0[i]-deltas[i]}
      })
#      print(paste('Descent iteration',it,'objective =',f(x0,args,lambda)))
    }
  }
  
  return(x0)
}


cv.lambda <- function(C,y,V,w0=NULL,nCrosses=5,lambdas=NULL,surrIndex=122){
  print('Starting lambda cross val')
  
  if (is.null(w0)){
    w0 <- glm(y~C,family='quasibinomial')$coefficients[-1]
    w0[is.na(w0)] <- 0
  }
  if (is.null(lambdas)){
    lambdas <- 10^seq(-3,-0.5,0.5)
  }
  
  n <- nrow(C)
  observedPats <- sample(n)
  valPats_overall <- lapply(1:nCrosses, function(i){observedPats[round(n*(i-1)/nCrosses+1):round(n*i/nCrosses)]})
  
  #grid <- sapply(1:nCrosses,function(i){
  grid <- foreach(i=1:nCrosses, .combine=cbind, .packages=c('lbfgs','glmnet','nlme','mixtools','expm','mvtnorm','lme4','geepack','foreach','doParallel','pROC','flexmix','PRROC','pROC','nloptr','bindata','randomForest','MASS','mixtools')) %do% {
    print(paste('Starting cross',i))
   
    validatePats <- valPats_overall[[i]]
    trainPats <- setdiff(observedPats,validatePats)
    
    sapply(lambdas,function(lambda){
      w_opt <- numericGradientDescent(w0,objective_w,args=list('C'=C[trainPats,],'y'=y[trainPats],'V'=V),
                                      constIndex=surrIndex,lambda=lambda,maxIt=50,tol=1e-4)
      objective_w(w_opt,args=list('C'=C[validatePats,],'y'=y[validatePats],'V'=V),lambda=0)
    })
  }
  #  })
  means <- rowMeans(grid)
  
  lambda_opt <- lambdas[which.min(means)]
  print(paste('lambda_opt =',lambda_opt))
  
  return(list('results'=grid,'lambda_opt'=lambda_opt))
}


hmmgp <- function(dat_train=NULL,dat_test=NULL,Cindices=NULL,w=NULL,w0=NULL,V=NULL,observed=NULL,
                  Estep=Estep_partial,Xtrain=NULL,Xtest=NULL,alpha=NULL,r=NULL,lambda=NULL,surrIndex=NULL){
  if (is.null(observed)){
    observed <- unique(dat_train$ID)
  }
  
  if (is.null(Xtrain)){
    Ctrain <- as.matrix(dat_train[,Cindices])
    Ctest <- as.matrix(dat_test[,Cindices])
    observedIndices <- which(dat_train$ID %in% observed)
    
    # Optimize w
    if (is.null(w)){
      if (is.null(w0)){
        w0 <- glm(dat_train$Y[observedIndices]~Ctrain[observedIndices,],family='quasibinomial')$coefficients[-1]
        w0[is.na(w0)] <- 0
        if (is.null(surrIndex)){
          surrIndex <- which.max(w0)
        }
        w0 <- w0/abs(w0[surrIndex])
      }
      w0 <- numericGradientDescent(w0, objective_w, args=list('C'=Ctrain[observedIndices,],'y'=dat_train$Y[observedIndices],'V'=V),
                                   constIndex=surrIndex, lambda=0, maxIt=100, tol=1e-3)
      
      # Optimize lambda
      if (is.null(lambda)){
        lambda <- cv.lambda(Ctrain[observedIndices,],dat_train$Y[observedIndices],V,w0,surrIndex=surrIndex)$lambda_opt
      }
      
      w <- numericGradientDescent(w0, objective_w, args=list('C'=Ctrain[observedIndices,],'y'=dat_train$Y[observedIndices],'V'=V),
                                  constIndex=surrIndex, lambda=lambda, maxIt=100, tol=1e-4)
    }
    
    # Define Xtrain, Xtest
    CWVtrain <- Ctrain %*% (w * V)
    CWVtest <- Ctest %*% (w * V)
    Xtrain <- data.frame(ID=dat_train$ID,Y=dat_train$Y,DD=log(dat_train$T+1),HID=exp(dat_train$H)-1,H=dat_train$H,T=dat_train$T)
    Xtrain <- cbind(Xtrain,CWVtrain); colnames(Xtrain)[-c(1:6)] <- paste0('X',seq(ncol(V)))
    Xtest <- data.frame(ID=dat_test$ID,Y=dat_test$Y,DD=log(dat_test$T+1),HID=exp(dat_test$H)-1,H=dat_test$H,T=dat_test$T)
    Xtest <- cbind(Xtest,CWVtest); colnames(Xtest)[-c(1:6)] <- paste0('X',seq(ncol(V)))
    Xtrain$pY <- Xtrain$pWeights <- 1
  }
  
  # Optimize r
  if (is.null(r)){
    r <- cv.r(Xtrain,observed,Estep=Estep)$r_opt
  }
  
  # Optimize alpha and predict
  if (is.null(alpha)){
    result <- lineSearch(Xtrain,observed,Xtest,r=r,Estep=Estep)
  }
  else{
    result <- lineSearch(Xtrain,observed,Xtest,alphas=alpha,r=r,Estep=Estep)
  }
  alpha <- result$alpha
  
  # Result
  return(list('w_opt'=w,'r_opt'=r,'alpha_opt'=alpha,'lambda_opt'=lambda,
              'margSup'=result$margSup,'margSemisup'=result$margSemisup,'margMix'=result$margMix,
              'resultSup'=result$resultSup,'resultSemisup'=result$resultSemisup,'resultMix'=result$resultMix,
              'cumSup'=result$cumSup,'cumSemisup'=result$cumSemisup,'cumMix'=result$cumMix))
}


# EVALUATION

ABC <- function(f1,f2,x,ID){
  keep <- which(!is.na(f1) & !is.na(f2))
  f1 <- f1[keep]; f2 <- f2[keep]; x <- x[keep]; ID <- ID[keep]
  mean(sapply(unique(ID),function(id){
    keep <- which(ID==id)
    f1 <- f1[keep]; f2 <- f2[keep]; x <- x[keep]
    if (length(x) < 2){NA}
    else{
      sum(sapply(2:length(x),function(i){
        if((f2[i-1]-f1[i-1])*(f2[i]-f1[i]) >= 0){
          (x[i]-x[i-1]) * mean(abs(f2[i-1]-f1[i-1]), abs(f2[i]-f1[i]))
        }
        else{
          0.5 * (x[i]-x[i-1]) * mean(abs(f2[i-1]-f1[i-1]), abs(f2[i]-f1[i]))
        }
      }))
    }
  }), na.rm=TRUE)
}


binToCount <- function(y,ID){
  unlist(sapply(unique(ID),function(id){
    y <- y[ID==id]
    count <- y
    if (length(count) > 1){
      for (i in 2:length(count)){
        count[i] <- count[i-1] + as.integer(y[i]*(1-y[i-1]))
      }
    }
    count
  }))
}

margToCum <- function(marg_hat,ID){
  out <- c(1-marg_hat[1], rep(NA,length(marg_hat)-1))
  for (i in 2:length(marg_hat)){
    if (ID[i] == ID[i-1]){
      out[i] <- out[i-1] * (1-marg_hat[i])
    }
    else{
      out[i] <- 1-marg_hat[i]
    }
  }
  return(1-out)
}

#evaluateFinal <- function(dat,marg_hat,cum_hat,prevalence,prevID){
#  ID <- dat$ID
#  cumID <- sapply(unique(ID),function(id){max(cum_hat[ID==id])})
#  maxMargProbs <- sapply(unique(ID),function(id){max(marg_hat[ID==id])})
#  yID <- as.integer(sapply(unique(ID),function(id){sum(dat$Y[ID==id])>0}))
#  yCount <- binToCount(dat$Y,dat$ID)
  
  # AUCs and Fscores @ FPR = 0.05
#  result_allT <- ROC.Est.FUN(dat$Y, marg_hat, yy0=0.5, fpr0=0.05)
#  result_allT$Fscore <- 2/(1/result_allT$TPR + 1/result_allT$PPV)
  
  # Area between curves
#  yCum <- c(dat$Y[1], rep(NA,nrow(dat)-1))
#  for (i in 2:nrow(dat)){yCum[i] <- ifelse(dat$ID[i]==dat$ID[i-1], max(dat$Y[i],yCum[i-1]), dat$Y[i])}
#  result_ID <- ROC.Est.FUN(yID, maxMargProbs, yy0=0.5, fpr0=0.05)
#  result_ID$Fscore <- 2/(1/result_ID$TPR + 1/result_ID$PPV)
#  result_ID$cumIDAUC <- auc(yID,cumID)
#  result_ID$cumABC <- ABC(yCum,cum_hat,dat$T,ID)
  
  # Thresholding at different values
#  prediction50 <- round(marg_hat)
#  prediction50ID <- as.integer(sapply(unique(ID),function(id){sum(prediction50[ID==id])>0}))
#  predictionPrev <- as.integer(marg_hat > quantile(marg_hat,1-prevalence))
#  predictionPrevID <- as.integer(marg_hat > quantile(maxMargProbs,1-prevID))
  
  # Count ABCs
#  result_ID$countABC50 <- ABC(yCount,binToCount(prediction50,dat$ID),dat$T,dat$ID)
#  result_ID$countABCprev <- ABC(yCount,binToCount(predictionPrev,dat$ID),dat$T,dat$ID)
#  result_ID$countABCprevID <- ABC(yCount,binToCount(predictionPrevID,dat$ID),dat$T,dat$ID)
  
  # TPR/FPR etc of thresholded values
#  result_allT$FPRPrev <- sum(predictionPrev*(1-dat$Y)) / sum(1-dat$Y)
#  result_allT$TPRPrev <- sum(predictionPrev*dat$Y) / sum(dat$Y)
#  result_allT$PPVPrev <- sum(predictionPrev*dat$Y) / sum(predictionPrev)
#  result_allT$NPVPrev <- sum((1-predictionPrev)*(1-dat$Y)) / sum(1-predictionPrev)
#  result_allT$FscorePrev <- 2/(1/result_allT$TPRPrev + 1/result_allT$PPVPrev)
#  predictionPrev_byID <- as.integer(sapply(unique(ID),function(id){sum(predictionPrev[ID==id])>0}))
#  result_ID$FPRPrev <- sum(predictionPrev_byID*(1-yID)) / sum(1-yID)
#  result_ID$TPRPrev <- sum(predictionPrev_byID*yID) / sum(yID)
#  result_ID$PPVPrev <- sum(predictionPrev_byID*yID) / sum(predictionPrevID)
#  result_ID$NPVPrev <- sum((1-predictionPrev_byID)*(1-yID)) / sum(1-predictionPrev_byID)
#  result_ID$FscorePrev <- 2/(1/result_ID$TPRPrev + 1/result_ID$PPVPrev)
  
#  result_allT$FPRPrevID <- sum(predictionPrevID*(1-dat$Y)) / sum(1-dat$Y)
#  result_allT$TPRPrevID <- sum(predictionPrevID*dat$Y) / sum(dat$Y)
#  result_allT$PPVPrevID <- sum(predictionPrevID*dat$Y) / sum(predictionPrevID)
#  result_allT$NPVPrevID <- sum((1-predictionPrevID)*(1-dat$Y)) / sum(1-predictionPrevID)
#  result_allT$FscorePrevID <- 2/(1/result_allT$TPRPrevID + 1/result_allT$PPVPrevID)
#  predictionPrevID_byID <- as.integer(sapply(unique(ID),function(id){sum(predictionPrevID[ID==id])>0}))
#  result_ID$FPRPrevID <- sum(predictionPrevID_byID*(1-yID)) / sum(1-yID)
#  result_ID$TPRPrevID <- sum(predictionPrevID_byID*yID) / sum(yID)
#  result_ID$PPVPrevID <- sum(predictionPrevID_byID*yID) / sum(predictionPrevID_byID)
#  result_ID$NPVPrevID <- sum((1-predictionPrevID_byID)*(1-yID)) / sum(1-predictionPrevID_byID)
#  result_ID$FscorePrevID <- 2/(1/result_ID$TPRPrevID + 1/result_ID$PPVPrevID)
  
  # How far off are the T estimates
#  actual_bySubject <- sapply(unique(ID),function(id){
#    matches <- ID==id
#    ifelse (sum(dat$Y[matches]) > 0, min(dat$T[matches & dat$Y==1]), NA)
#  })
#  prediction50_bySubject <- sapply(unique(ID),function(id){
#    matches <- ID==id
#    ifelse (sum(prediction50[matches]) > 0, min(dat$T[matches & prediction50==1]), NA)
#  })
#  predictionPrev_bySubject <- sapply(unique(ID),function(id){
#    matches <- ID==id
#    ifelse (sum(predictionPrev[matches]) > 0, min(dat$T[matches & predictionPrev==1]), NA)
#  })
#  predictionPrevID_bySubject <- sapply(unique(ID),function(id){
#    matches <- ID==id
#    ifelse (sum(predictionPrevID[matches]) > 0, min(dat$T[matches & predictionPrevID==1]), NA)
#  })
  
#  result_ID$Tmed50 <- median(prediction50_bySubject-actual_bySubject,na.rm=TRUE)
#  result_ID$Tmed50Abs <- median(abs(prediction50_bySubject-actual_bySubject),na.rm=TRUE)
#  result_ID$TmedPrev <- median(predictionPrev_bySubject-actual_bySubject,na.rm=TRUE)
#  result_ID$TmedPrevAbs <- median(abs(predictionPrev_bySubject-actual_bySubject),na.rm=TRUE)
#  result_ID$TmedPrevID <- median(predictionPrevID_bySubject-actual_bySubject,na.rm=TRUE)
#  result_ID$TmedPrevIDAbs <- median(abs(predictionPrevID_bySubject-actual_bySubject),na.rm=TRUE)
  
#  result_ID$rangeAccuracies50 <- sapply(0:12,function(delta){mean(prediction50_bySubject-actual_bySubject <= delta, na.rm=TRUE)})
#  result_ID$rangeAccuraciesPrev <- sapply(0:12,function(delta){mean(predictionPrev_bySubject-actual_bySubject <= delta, na.rm=TRUE)})
#  result_ID$rangeAccuraciesPrevID <- sapply(0:12,function(delta){mean(predictionPrevID_bySubject-actual_bySubject <= delta, na.rm=TRUE)})
  
#  return(list("result_allT"=result_allT,"result_ID"=result_ID))
#}


evaluateFinal <- function(dat,marg_hat,cum_hat,prevalence){
  ID <- dat$ID
  cumID <- sapply(unique(ID),function(id){max(cum_hat[ID==id])})
  maxMargProbs <- sapply(unique(ID),function(id){max(marg_hat[ID==id])})
  yID <- as.integer(sapply(unique(ID),function(id){sum(dat$Y[ID==id])>0}))
  yCount <- binToCount(dat$Y,dat$ID)
  
  # AUCs and Fscores @ FPR = 0.05
  result_allT <- ROC.Est.FUN(dat$Y, marg_hat, yy0=0.5, fpr0=0.05)
  result_allT$Fscore <- 2/(1/result_allT$TPR + 1/result_allT$PPV)
  
  # Area between curves
  yCum <- c(dat$Y[1], rep(NA,nrow(dat)-1))
  for (i in 2:nrow(dat)){yCum[i] <- ifelse(dat$ID[i]==dat$ID[i-1], max(dat$Y[i],yCum[i-1]), dat$Y[i])}
  result_ID <- ROC.Est.FUN(yID, maxMargProbs, yy0=0.5, fpr0=0.05)
  result_ID$Fscore <- 2/(1/result_ID$TPR + 1/result_ID$PPV)
  result_ID$cumIDAUC <- auc(yID,cumID)
  result_ID$cumABC <- ABC(yCum,cum_hat,dat$T,ID)
  
  # Thresholding at different values
  predictionPrev <- as.integer(marg_hat > quantile(marg_hat,1-prevalence))
 
  # Count ABCs
  result_ID$countABCprev <- ABC(yCount,binToCount(predictionPrev,dat$ID),dat$T,dat$ID)

  # TPR/FPR etc of thresholded values
  result_allT$FPRPrev <- sum(predictionPrev*(1-dat$Y)) / sum(1-dat$Y)
  result_allT$TPRPrev <- sum(predictionPrev*dat$Y) / sum(dat$Y)
  result_allT$PPVPrev <- sum(predictionPrev*dat$Y) / sum(predictionPrev)
  result_allT$NPVPrev <- sum((1-predictionPrev)*(1-dat$Y)) / sum(1-predictionPrev)
  result_allT$FscorePrev <- 2/(1/result_allT$TPRPrev + 1/result_allT$PPVPrev)
  predictionPrev_byID <- as.integer(sapply(unique(ID),function(id){sum(predictionPrev[ID==id])>0}))
  result_ID$FPRPrev <- sum(predictionPrev_byID*(1-yID)) / sum(1-yID)
  result_ID$TPRPrev <- sum(predictionPrev_byID*yID) / sum(yID)
  result_ID$PPVPrev <- sum(predictionPrev_byID*yID) / sum(predictionPrev_byID)
  result_ID$NPVPrev <- sum((1-predictionPrev_byID)*(1-yID)) / sum(1-predictionPrev_byID)
  result_ID$FscorePrev <- 2/(1/result_ID$TPRPrev + 1/result_ID$PPVPrev)
  
  return(list("result_allT"=result_allT,"result_ID"=result_ID))
}


## TEST USE ##

#cumProbs <- sapply(unique(ID),function(id){max(cumulative[ID==id])})
#maxMargProbs <- sapply(unique(ID),function(id){max(marginal_partial[ID==id])})
#yID <- as.integer(sapply(unique(ID),function(id){sum(labeled$Y[ID==id])>0}))

#observed <- sample(labeled$ID,100)
#trained <- Mstep(labeled[labeled$ID %in% observed,],r=0.5)
#marginal_partial <- Estep_partial(labeled[!(labeled$ID %in% observed),],trained)
#cumulative <- cumulative_prob(labeled[!(labeled$ID %in% observed),],trained)

#partial_results <- evaluate(labeled[!(labeled$ID %in% observed),],marginal_partial,cumulative)

#phecode_cumhat <- as.numeric(sapply(1:nrow(CC_comb_cohort),function(i){sum(CC_comb_cohort$PheCode.335_[which(CC_comb_cohort$PatientNum[1:i]==CC_comb_cohort$PatientNum[i])])>0}))
#cuiMS_cumhat <- as.numeric(sapply(1:nrow(CC_comb_cohort),function(i){sum(CC_comb_cohort$CUI.C0026769[which(CC_comb_cohort$PatientNum[1:i]==CC_comb_cohort$PatientNum[i])])>0}))
#cuiRecurr_cumhat <- as.numeric(sapply(1:nrow(CC_comb_cohort),function(i){sum(CC_comb_cohort$CUI.C0277556[which(CC_comb_cohort$PatientNum[1:i]==CC_comb_cohort$PatientNum[i])])>0}))

#phecode_results <- evaluate(labeled,log(log(CC_comb_cohort$PheCode.335_+1)+1),phecode_cumhat)
#cuiMS_results <- evaluate(labeled,log(log(CC_comb_cohort$CUI.C0026769+1)+1),cuiMS_cumhat)
#cuiRelapse_results <- evaluate(labeled,log(log(CC_comb_cohort$CUI.C0277556+1)+1),cuiRecurr_cumhat)


