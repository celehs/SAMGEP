source('SAMGEP.R')
source('simulate.R')
load('trained_model.RData')

simulated <- simulate_samgep(1000,25,trained)
samgep_prediction <- samgep(dat_train=simulated,Cindices=7:17,
                            observed=1:100,V=simulated$)


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
