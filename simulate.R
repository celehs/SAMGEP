simulate_samgep <- function(N,Tn,trained,nX=11,dist='normal'){
  priorModel <- trained$priorModel; likeModel <- trained$likeModel; transCoefs <- trained$transCoefs
  Tns <- rpois(N,Tn)
  
  result <- foreach(i=1:N, .combine=rbind) %do% {
    H <- 1 + rchisq(1,2.2)
    Hlog <- log(HID+1)
    
    Tlog <- c(0,rep(NA,Tns[i]-1))
    T <- c(1,rep(NA,Tns[i]-1))
    for (j in 2:Tns[i]){
      delta <- rexp(1,3)
      Tlog[j] <- Tlog[j-1] + delta
      T[j] <- T[j-1] + ceiling(2*delta)
    }
    
    Yn <- c(rbinom(1,1,expit(c(1,Hlog,T[1],Tlog[1])%*%priorModel$coefficients)),rep(NA,Tns[i]-1))
    for (j in 2:Tns[i]){
      Yn[j] <- rbinom(1,1,expit(transCoefs%*%c(1,Yn[j-1],Hlog,Tlog[j],T[j],T[j]-T[j-1]==1,Yn[j-1]*(T[j]-T[j-1]==1))))
    }
    
    Hn <- as.matrix(cbind(1,Yn,Hlog,T,Tlog,Yn*H,Yn*T,Yn*Tlog))
    mu <- c(likeModel$beta %*% t(Hn))
    
    sig <- likeModel$sigma * HID^likeModel$alpha
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
    
    cbind(i,HID,H,T,log(Tlog+1),Yn,Xn,1,1)
  }
  
  colnames(result) <- c('ID','H','Hlog','T','Tlog','Y',paste0('X',seq(nX)),'pInv','pY')
  return(as.data.frame(result))
}


simulate_part2 <- function(C,V,trueDat,Estep=Estep_partial){
  d <- ncol(C)
  ratios <- sapply(1:ncol(C),function(i){
    (mean(C[trueDat$Y==1,i])-mean(C[trueDat$Y==0,i]))/sd(C[,i])
  })
  W <- sign(ratios) * (abs(ratios)>0.25) * rnorm(d,0,.25)
  
  X <- as.matrix(C) %*% (W*as.matrix(V))
  
  dat <- as.data.frame(cbind(trueDat[,c('ID','Y','T','Tlog','H','Hlog','X11')],X))
  colnames(dat) <- c('ID','Y','T','Tlog','H','Hlog','X11',paste0('X',1:10))
  dat$pY <- dat$pInv <- 1
  
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

