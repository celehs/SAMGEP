`%do%` <- foreach::`%do%`
`%dopar%` <- foreach::`%dopar%`
utils::globalVariables(c("it","lambda","r","pY","i"))

# SAMGEP.R: Contains samgep function. See Ahuja et al. (2020), BioArxiv for details.
# Author: Yuri Ahuja
# Last Updated: 11/25/2020

# INPUT:
# dat_train = Raw training data set, including patient IDs (ID), healthcare utilization feature (H) and censoring time (C) (optional if Xtrain is supplied)
# dat_test = Raw testing data set, including patient IDs (ID), a healthcare utilization feature (H) and censoring time (C) (optional)
# Cindices = Column indices of EHR feature counts in dat_train/dat_test (optional if Xtrain is supplied)
# w = Pre-optimized EHR feature weights (optional if Xtrain is supplied)
# w0 = Initial (i.e. partially optimized) EHR feature weights (optional if Xtrain is supplied)
# V = nFeatures x nEmbeddings embeddings matrix (optional if Xtrain is supplied)
# observed = IDs of patients with observed outcome labels (optional if Xtrain is supplied)
# nX = Number of embedding features (defaults to 10)
# Estep = E-step function to use (Estep_partial or Estep_full; defaults to Estep_partial)
# Xtrain = Embedded training data set, including patient IDs (ID), healthcare utilization feature (H) and censoring time (C) (optional)
# Xtest = Embedded testing data set, including patient IDs (ID), healthcare utilization feature (H) and censoring time (C) (optional)
# alpha = Relative weight of semi-supervised to supervised MGP predictors in SAMGEP ensemble (optional)
# r = Scaling factor of inter-temporal correlation (optional)
# lambda = L1 regularization hyperparameter for feature weight (w) optimization (optional)
# surrIndex = Index (within Cindices) of primary surrogate index for outcome event (optional)
# nCores = Number of cores to use for parallelization (defaults to 1)


fitGLS <- function(dat, nX = 10, r = 1, covs=NULL) {
  repInds <- which(duplicated(dat$ID))

  coefs <- lapply(1:nX, function(i) {
    tryCatch({
      model <- nlme::gls(as.formula(paste0("X", i, " ~ Y+Hlog+T+Tlog+Y:Hlog+Y:T+Y:Tlog")),
                         data = dat,
                         weights = varComb(varPower(form = ~H), varFixed(~pInv))
      )
      list(
        "beta" = model$coefficients, "alpha" = unlist(model$modelStruct$varStruct), "sigma" = model$sigma,
        "resids" = model$residuals, "std.errors" = model$varBeta
      )
    }, error=function(e){
      model <- nlme::gls(as.formula(paste0("X", i, " ~ Y+Hlog+T+Tlog+Y:Hlog+Y:T+Y:Tlog")),
                         data = dat,
                         weights = varFixed(~pInv)
      )
      list(
        "beta" = model$coefficients, "alpha" = 0, "sigma" = model$sigma,
        "resids" = model$residuals, "std.errors" = model$varBeta
      )
    })

  })

  coefs$beta <- t(sapply(1:nX, function(i) {
    coefs[[i]]$beta
  }))
  coefs$alpha <- sapply(1:nX, function(i) {
    coefs[[i]]$alpha
  })
  coefs$sigma <- sapply(1:nX, function(i) {
    coefs[[i]]$sigma
  }) * sqrt((nrow(dat) - 8) / (sum(dat$pY) - 8))
  coefs$std.errors <- abind(lapply(1:nX, function(i) {
    coefs[[i]]$std.errors
  }), along = 3)

  resids_norm <- sapply(1:nX, function(i) {
    coefs[[i]]$resids / (dat$H^coefs[[i]]$alpha)
  })
  covMat <- 1 / (sum(dat$pY) - 8) * (t(resids_norm * dat$pY) %*% resids_norm)
  sigma <- sqrt(diag(covMat))
  coefs$corrMat <- covMat / (sigma %*% t(sigma))

  # Vector auto-regression (1 lag)
  coefs$autoCoefs <- r * sapply(1:nX, function(i) {
    rep(lm(resids_norm[repInds, i] ~ resids_norm[repInds - 1, i]:as.factor(dat$T[repInds] - dat$T[repInds - 1] == 1),
      weights = (dat$pY[repInds] * dat$pY[repInds - 1])
    )$coef[-1], 2)
  })

  return(coefs)
}


trainTransitionMatrix <- function(train) {
  counts <- matrix(rep(0, 4), 2, 2)
  rownames(counts) <- c("Start0", "Start1")
  colnames(counts) <- c("End0", "End1")

  for (i in 2:nrow(train)) {
    if (train$ID[i] == train$ID[i - 1]) {
      counts <- counts + train$pY[i] * outer(c(1 - train$Y[i - 1], train$Y[i - 1]), c(1 - train$Y[i], train$Y[i]))
    }
  }

  transition <- counts / rowSums(counts)
  return(list("counts" = counts, "tmat" = transition))
}


Estep_partial <- function(dat, trained, nX = 10) {
  expit <- function(x) {
    1 / (1 + exp(-x))
  }

  priorModel <- trained$priorModel
  likeModel <- trained$likeModel
  transCoefs <- trained$transCoefs
  prior_fitted <- predict(priorModel, dat, type = "response")

  a1 <- likeModel$autoCoefs[3:4, ]
  a2 <- likeModel$autoCoefs[1:2, ]
  amin1 <- apply(a1, 2, min)
  amin2 <- apply(a2, 2, min)
  
  post <- c()
  for (pat in unique(dat$ID)) {
    # Find patient-specific parameters
    patIdx <- which(dat$ID == pat)
    Ti <- dat$T[patIdx]
    Tlogi <- dat$Tlog[patIdx]
    Hi <- unique(dat$H[patIdx])
    Hlogi <- unique(dat$Hlog[patIdx])
    prior_i <- prior_fitted[patIdx]
    
    # Compute mu | Y=0; mu | Y=1
    mu0 <- likeModel$beta %*% rbind(1, 0, Hlogi, Ti, Tlogi, 0, 0, 0)
    mu1 <- likeModel$beta %*% rbind(1, 1, Hlogi, Ti, Tlogi, Hlogi, Ti, Tlogi)
    mus <- list(mu0, mu1)
    
    # Compute marginal, conditional sigma_squareds
    sigsq_like <- (likeModel$sigma * Hi^likeModel$alpha)^2
    sigsq_prior <- sapply(1:nX, function(i) {
      c(1, 1, Hlogi, mean(Ti), mean(Tlogi), Hlogi, mean(Ti), mean(Tlogi)) %*%
        likeModel$std.errors[, , i] %*% c(1, 1, Hlogi, mean(Ti), mean(Tlogi), Hlogi, mean(Ti), mean(Tlogi))
    })
    sig <- sqrt(sigsq_like + sigsq_prior)
    Sigma <- likeModel$corrMat * (sig %*% t(sig))
    A1 <- sqrt(1 - amin1^2) %*% t(sqrt(1 - amin1^2))
    A2 <- sqrt(1 - amin2^2) %*% t(sqrt(1 - amin2^2))
    Sigma_cond1 <- A1 * Sigma
    Sigma_cond2 <- A2 * Sigma
    
    Xi <- as.matrix(dat[patIdx, paste0("X", 1:nX)])
    
    post_i <- c()
    for (t in 1:length(Ti)) {
      # Patient has 1 timepoint
      if (length(Ti) == 1) {
        logprior <- log(c(1 - prior_i[t], prior_i[t]))
        logprobs <- logprior + c(
          dmvnrm_arma_fast(t(Xi[t, ]), t(mus[[1]][, t]), Sigma, TRUE),
          dmvnrm_arma_fast(t(Xi[t, ]), t(mus[[2]][, t]), Sigma, TRUE)
        )
        probs <- exp(logprobs - max(logprobs))
        post_i <- c(post_i, probs[2] / sum(probs), rep(0, 4))
      }
      
      # Currently on patient's first timepoint
      else if (t == 1) {
        X_comp <- c(Xi[t, ], Xi[t + 1, ])
        prior_it <- c(1 - prior_i[t], prior_i[t])
        if (Ti[t + 1] - Ti[t] == 1) {
          a <- a1
          Sigma_cond <- Sigma_cond1
        }
        else {
          a <- a2
          Sigma_cond <- Sigma_cond2
        }
        
        logprobs <- matrix(NA, 2, 2)
        for (yn in 1:2) {
          for (yc in 1:2) {
            p <- expit(transCoefs %*% c(1, yc - 1, Hi, Ti[t + 1], Tlogi[t + 1], as.integer(Ti[t + 1] - Ti[t] == 1), as.integer(Ti[t + 1] - Ti[t] == 1) * (yc - 1)))
            logprior <- log(prior_it[yc] * c(1 - p, p)[yn])
            muC <- mus[[yc]][, t]
            muN <- mus[[yn]][, t + 1] + a[1 + (yc == yn), ] * (Xi[t, ] - muC)
            loglike <- dmvnrm_arma_fast(t(Xi[t, ]), t(muC), Sigma, TRUE) + dmvnrm_arma_fast(t(Xi[t + 1, ]), t(muN), Sigma_cond, TRUE)
            logprobs[yc, yn] <- logprior + loglike
          }
        }
        probs <- exp(logprobs - max(logprobs))
        post_i <- c(post_i, sum(probs[2, ]) / sum(probs), rep(0, 4))
      }
      
      # Currently on patient's last timepoint
      else if (t == length(Ti)) {
        X_comp <- c(Xi[t - 1, ], Xi[t, ])
        prior_it <- c(1 - prior_i[t - 1], prior_i[t - 1])
        if (Ti[t] - Ti[t - 1] == 1) {
          a <- a1
          Sigma_cond <- Sigma_cond1
        }
        else {
          a <- a2
          Sigma_cond <- Sigma_cond2
        }
        
        logprobs <- matrix(NA, 2, 2)
        for (yn in 1:2) {
          for (yc in 1:2) {
            p <- expit(transCoefs %*% c(1, yc - 1, Hi, Ti[t], log(Ti[t] + 1), as.integer(Ti[t] - Ti[t - 1] == 1), as.integer(Ti[t] - Ti[t - 1] == 1) * (yc - 1)))
            logprior <- log(prior_i[yc] * c(1 - p, p)[yn])
            muC <- mus[[yc]][, t - 1]
            muN <- mus[[yn]][, t] + a[1 + (yc == yn), ] * (Xi[t - 1, ] - muC)
            loglike <- dmvnrm_arma_fast(t(Xi[t - 1, ]), t(muC), Sigma, TRUE) + dmvnrm_arma_fast(t(Xi[t, ]), t(muN), Sigma_cond, TRUE)
            logprobs[yc, yn] <- logprior + loglike
          }
        }
        probs <- exp(logprobs - max(logprobs))
        post_i <- c(post_i, sum(probs[, 2]) / sum(probs), probs / sum(probs))
      }
      
      # Currently on intermediate timepoint
      else {
        X_comp <- c(Xi[t - 1, ], Xi[t, ], Xi[t + 1, ])
        prior_it <- c(1 - prior_i[t - 1], prior_i[t - 1])
        if (Ti[t] - Ti[t - 1] == 1) {
          aL <- a1
          Sigma_condL <- Sigma_cond1
        }
        else {
          aL <- a2
          Sigma_condL <- Sigma_cond2
        }
        if (Ti[t + 1] - Ti[t] == 1) {
          aN <- a1
          Sigma_condN <- Sigma_cond1
        }
        else {
          aN <- a2
          Sigma_condN <- Sigma_cond2
        }
        
        logprobs <- array(NA, dim = c(2, 2, 2))
        for (yn in 1:2) {
          for (yc in 1:2) {
            for (yl in 1:2) {
              p1 <- expit(transCoefs %*% c(1, yl - 1, Hi, Ti[t], Tlogi[t], as.integer(Ti[t] - Ti[t - 1] == 1), as.integer(Ti[t] - Ti[t - 1] == 1) * (yl - 1)))
              p2 <- expit(transCoefs %*% c(1, yc - 1, Hi, Ti[t + 1], Tlogi[t + 1], as.integer(Ti[t + 1] - Ti[t] == 1), as.integer(Ti[t + 1] - Ti[t] == 1) * (yc - 1)))
              logprior <- log(prior_it[yl] * c(1 - p1, p1)[yc] * c(1 - p2, p2)[yn])
              muL <- mus[[yl]][, t - 1]
              muC <- mus[[yc]][, t] + aL[1 + (yl == yc), ] * (Xi[t - 1, ] - muL)
              muN <- mus[[yn]][, t + 1] + aN[1 + (yc == yn), ] * (Xi[t, ] - muC)
              loglike <- dmvnrm_arma_fast(t(Xi[t - 1, ]), t(muL), Sigma, TRUE) + 
                dmvnrm_arma_fast(t(Xi[t, ]), t(muC), Sigma_condL, TRUE) + 
                dmvnrm_arma_fast(t(Xi[t + 1, ]), t(muN), Sigma_condN, TRUE)
              logprobs[yl, yc, yn] <- logprior + loglike
            }
          }
        }
        probs <- exp(logprobs - max(logprobs))
        post_i <- c(post_i, sum(probs[, 2, ]) / sum(probs), probs[, , 1] / sum(probs[, , 1]))
      }
    }
    
    post <- c(post, post_i)
  }
  
  # Reformat and save output
  output <- post[seq(1, length(post), 5)]
  attr(output, "post2") <- array(post[-seq(1, length(post), 5)], dim = c(2, 2, length(output)))
  output
}


Mstep <- function(train, train_interTemp = NULL, nX = 10, r = 1, survival = FALSE) {
  priorModel <- glm(Y ~ Hlog + T + Tlog, weights = pY, data = train, family = "quasibinomial")

  if (is.null(train_interTemp)) {
    repInds <- which(duplicated(train$ID))
    transCoefs <- glm(train$Y[repInds] ~ train$Y[repInds - 1] + train$Hlog[repInds] + train$T[repInds] + train$Tlog[repInds] +
      as.factor(train$T[repInds] - train$T[repInds - 1] == 1) + as.factor(train$T[repInds] - train$T[repInds - 1] == 1):train$Y[repInds - 1],
    family = "quasibinomial"
    )$coefficients
  }
  else {
    transCoefs <- glm(Ycurr ~ Yprev + H + T + Tlog + as.factor(T - Tprev == 1) + Yprev:as.factor(T - Tprev == 1),
      weights = pY, data = train_interTemp, family = "quasibinomial")$coefficients
  }
  
  if (survival){
    transCoefs[2] <- 1000
  }

  likeModel <- fitGLS(train, nX = nX, r = r)

  return(list("priorModel" = priorModel, "transCoefs" = transCoefs, "likeModel" = likeModel))
}


# Full Markov implementation
Estep_full <- function(dat, trained, nX = 10) {
  expit <- function(x) {
    1 / (1 + exp(-x))
  }

  priorModel <- trained$priorModel
  likeModel <- trained$likeModel
  transCoefs <- trained$transCoefs
  a1 <- likeModel$autoCoefs[3:4, ]
  a2 <- likeModel$autoCoefs[1:2, ]
  amin1 <- apply(a1, 2, min)
  amin2 <- apply(a2, 2, min)
  prior_fitted <- predict(priorModel, dat, type = "response")

  post1 <- unlist(sapply(unique(dat$ID), function(patient) {
    encounters <- which(dat$ID == patient)
    Ti <- dat$T[encounters]
    Tlogi <- dat$Tlog[encounters]
    Hi <- unique(dat$H[encounters])
    Hlogi <- unique(dat$Hlog[encounters])
    prior_pat <- prior_fitted[encounters]

    mu0 <- likeModel$beta %*% rbind(1, 0, Hlogi, Ti, Tlogi, 0, 0, 0)
    mu1 <- likeModel$beta %*% rbind(1, 1, Hlogi, Ti, Tlogi, Hlogi, Ti, Tlogi)
    mus <- list(mu0, mu1)

    sigsq_base <- (likeModel$sigma * Hi^likeModel$alpha)^2
    sigsq_prior <- sapply(1:nX, function(i) {
      c(1, 1, Hlogi, mean(Ti), mean(Tlogi), Hlogi, mean(Ti), mean(Tlogi)) %*%
        likeModel$std.errors[, , i] %*% c(1, 1, Hlogi, mean(Ti), mean(Tlogi), Hlogi, mean(Ti), mean(Tlogi))
    })
    sig <- sqrt(sigsq_base + sigsq_prior)
    Sigma <- likeModel$corrMat * (sig %*% t(sig))
    A1 <- sqrt(1 - amin1^2) %*% t(sqrt(1 - amin1^2))
    A2 <- sqrt(1 - amin2^2) %*% t(sqrt(1 - amin2^2))
    Sigma_cond1 <- A1 * Sigma
    Sigma_cond2 <- A2 * Sigma

    Xi <- as.matrix(dat[encounters, paste0("X", 1:nX)])

    fwd <- matrix(0, 2, length(encounters))
    f_prev <- c(1 - prior_pat[1], prior_pat[1])
    for (i in 1:length(encounters)) {
      if (i == 1) {
        loglike <- c(
          dmvnrm_arma_fast(t(Xi[i, ]), t(mu0[, i]), Sigma),
          dmvnrm_arma_fast(t(Xi[i, ]), t(mu1[, i]), Sigma)
        )
        fwd[, i] <- f_prev * exp(loglike - max(loglike))
      }
      else {
        if (Ti[i] - Ti[i - 1] == 1) {
          a <- a1
          Sigma_cond <- Sigma_cond1
        }
        else {
          a <- a2
          Sigma_cond <- Sigma_cond2
        }
        p <- expit(cbind(1, c(0, 1), Hi, Ti[i], log(Ti[i] + 1), as.integer(Ti[i] - Ti[i - 1] == 1), as.integer(Ti[i] - Ti[i - 1] == 1) * c(0, 1)) %*% transCoefs)
        transition <- cbind(1 - p, p)
        probs <- (t(transition) * rbind(f_prev, f_prev))
        loglike <- matrix(c(
          dmvnrm_arma_fast(t(Xi[i, ]), t(mu0[, i] + a[2, ] * (Xi[i - 1, ] - mu0[, i - 1])), Sigma_cond, TRUE),
          dmvnrm_arma_fast(t(Xi[i, ]), t(mu1[, i] + a[1, ] * (Xi[i - 1, ] - mu0[, i - 1])), Sigma_cond, TRUE),
          dmvnrm_arma_fast(t(Xi[i, ]), t(mu0[, i] + a[1, ] * (Xi[i - 1, ] - mu1[, i - 1])), Sigma_cond, TRUE),
          dmvnrm_arma_fast(t(Xi[i, ]), t(mu1[, i] + a[2, ] * (Xi[i - 1, ] - mu1[, i - 1])), Sigma_cond, TRUE)
        ), 2, 2)
        probs <- probs * exp(loglike - max(loglike))
        fwd[, i] <- probs %*% c(1, 1)
      }

      fwd[, i] <- fwd[, i] / sum(fwd[, i])
      f_prev <- fwd[, i]
    }

    bkw <- matrix(0, 2, length(encounters))
    b_next <- bkw[, length(encounters)] <- c(1, 1)
    if (length(encounters) > 1) {
      for (i in (length(encounters) - 1):1) {
        if (Ti[i + 1] - Ti[i] == 1) {
          a <- a1
          Sigma_cond <- Sigma_cond1
        }
        else {
          a <- a2
          Sigma_cond <- Sigma_cond2
        }
        p <- expit(cbind(1, c(0, 1), Hlogi, Ti[i + 1], Tlogi[i + 1], as.integer(Ti[i + 1] - Ti[i] == 1), as.integer(Ti[i + 1] - Ti[i] == 1) * c(0, 1)) %*% transCoefs)
        transition <- cbind(1 - p, p)

        probs <- (transition * rbind(b_next, b_next))
        loglike <- matrix(c(
          dmvnrm_arma_fast(t(Xi[i + 1, ]), t(mu0[, i + 1] + a[2, ] * (Xi[i, ] - mu0[, i])), Sigma_cond, TRUE),
          dmvnrm_arma_fast(t(Xi[i + 1, ]), t(mu1[, i + 1] + a[1, ] * (Xi[i, ] - mu0[, i])), Sigma_cond, TRUE),
          dmvnrm_arma_fast(t(Xi[i + 1, ]), t(mu0[, i + 1] + a[1, ] * (Xi[i, ] - mu1[, i])), Sigma_cond, TRUE),
          dmvnrm_arma_fast(t(Xi[i + 1, ]), t(mu1[, i + 1] + a[2, ] * (Xi[i, ] - mu1[, i])), Sigma_cond, TRUE)
        ), 2, 2)
        probs <- probs * exp(loglike - max(loglike))

        bkw[, i] <- probs %*% c(1, 1)
        bkw[, i] <- bkw[, i] / sum(bkw[, i])
        b_next <- bkw[, i]
      }
    }

    post <- fwd[2, ] * bkw[2, ] / colSums(fwd * bkw)
    post[colSums(fwd * bkw) == 0] <- prior_pat[colSums(fwd * bkw) == 0]
    post
  }))
  attr(post1, "post2") <- array(cbind(0, sapply(2:length(post1), function(i) {
    c(1 - post1[i - 1], post1[i - 1]) %*% t(c(1 - post1[i], post1[i]))
  })), dim = c(2, 2, length(post1)))
  post1
}


cumulative_prob <- function(dat, trained, nX = 10) {
  expit <- function(x) {
    1 / (1 + exp(-x))
  }

  priorModel <- trained$priorModel
  likeModel <- trained$likeModel
  transCoefs <- trained$transCoefs
  a1 <- likeModel$autoCoefs[3:4, ]
  a2 <- likeModel$autoCoefs[1:2, ]
  amin1 <- apply(a1, 2, min)
  amin2 <- apply(a2, 2, min)
  prior_fitted <- predict(priorModel, dat, type = "response")

  cumprob <- unlist(sapply(unique(dat$ID), function(patient) {
    encounters <- which(dat$ID == patient)
    Ti <- dat$T[encounters]
    Tlogi <- dat$Tlog[encounters]
    Hi <- unique(dat$H[encounters])
    Hlogi <- unique(dat$Hlog[encounters])
    prior_pat <- prior_fitted[encounters]

    mu0 <- likeModel$beta %*% rbind(1, 0, Hlogi, Ti, Tlogi, 0, 0, 0)
    mu1 <- likeModel$beta %*% rbind(1, 1, Hlogi, Ti, Tlogi, Hlogi, Ti, Tlogi)
    mus <- list(mu0, mu1)

    sigsq_base <- (likeModel$sigma * Hi^likeModel$alpha)^2
    sigsq_prior <- sapply(1:nX, function(i) {
      c(1, 1, Hlogi, mean(Ti), mean(Tlogi), Hlogi, mean(Ti), mean(Tlogi)) %*%
        likeModel$std.errors[, , i] %*% c(1, 1, Hlogi, mean(Ti), mean(Tlogi), Hlogi, mean(Ti), mean(Tlogi))
    })
    sig <- sqrt(sigsq_base + sigsq_prior)
    Sigma <- likeModel$corrMat * (sig %*% t(sig))
    A1 <- sqrt(1 - amin1^2) %*% t(sqrt(1 - amin1^2))
    A2 <- sqrt(1 - amin2^2) %*% t(sqrt(1 - amin2^2))
    Sigma_cond1 <- A1 * Sigma
    Sigma_cond2 <- A2 * Sigma

    Xi <- as.matrix(dat[encounters, paste0("X", 1:nX)])

    fwd <- matrix(0, 2, length(encounters))
    f_prev <- c(1 - prior_pat[1], prior_pat[1])
    survival <- rep(0, length(encounters))
    for (i in 1:length(encounters)) {
      if (i == 1) {
        loglike <- c(
          dmvnrm_arma_fast(t(Xi[i, ]), t(mu0[, i]), Sigma),
          dmvnrm_arma_fast(t(Xi[i, ]), t(mu1[, i]), Sigma)
        )
        fwd[, i] <- f_prev * exp(loglike - max(loglike))
        survival[i] <- fwd[1, i]
      }
      else {
        if (Ti[i] - Ti[i - 1] == 1) {
          a <- a1
          Sigma_cond <- Sigma_cond1
        }
        else {
          a <- a2
          Sigma_cond <- Sigma_cond2
        }
        p <- expit(cbind(1, c(0, 1), Hi, Ti[i], log(Ti[i] + 1), as.integer(Ti[i] - Ti[i - 1] == 1), as.integer(Ti[i] - Ti[i - 1] == 1) * c(0, 1)) %*% transCoefs)
        transition <- cbind(1 - p, p)
        probs <- (t(transition) * rbind(f_prev, f_prev))
        loglike <- matrix(c(
          dmvnrm_arma_fast(t(Xi[i, ]), t(mu0[, i] + a[2, ] * (Xi[i - 1, ] - mu0[, i - 1])), Sigma_cond, TRUE),
          dmvnrm_arma_fast(t(Xi[i, ]), t(mu1[, i] + a[1, ] * (Xi[i - 1, ] - mu0[, i - 1])), Sigma_cond, TRUE),
          dmvnrm_arma_fast(t(Xi[i, ]), t(mu0[, i] + a[1, ] * (Xi[i - 1, ] - mu1[, i - 1])), Sigma_cond, TRUE),
          dmvnrm_arma_fast(t(Xi[i, ]), t(mu1[, i] + a[2, ] * (Xi[i - 1, ] - mu1[, i - 1])), Sigma_cond, TRUE)
        ), 2, 2)
        like <- exp(loglike - max(loglike))
        probs <- probs * like
        fwd[, i] <- probs %*% c(1, 1)
        survival[i] <- survival[i - 1] * transition[1, 1] * like[1, 1]
      }

      survival[i] <- survival[i] / sum(fwd[, i])
      fwd[, i] <- fwd[, i] / sum(fwd[, i])
      f_prev <- fwd[, i]
    }

    bkw <- matrix(0, 2, length(encounters))
    b_next <- bkw[, length(encounters)] <- c(1, 1)
    if (length(encounters) > 1) {
      for (i in (length(encounters) - 1):1) {
        if (Ti[i + 1] - Ti[i] == 1) {
          a <- a1
          Sigma_cond <- Sigma_cond1
        }
        else {
          a <- a2
          Sigma_cond <- Sigma_cond2
        }
        p <- expit(cbind(1, c(0, 1), Hi, Ti[i], log(Ti[i] + 1), as.integer(Ti[i + 1] - Ti[i] == 1), as.integer(Ti[i + 1] - Ti[i] == 1) * c(0, 1)) %*% transCoefs)
        transition <- cbind(1 - p, p)

        probs <- (transition * rbind(b_next, b_next))
        loglike <- matrix(c(
          dmvnrm_arma_fast(t(Xi[i + 1, ]), t(mu0[, i + 1] + a[2, ] * (Xi[i, ] - mu0[, i])), Sigma_cond, TRUE),
          dmvnrm_arma_fast(t(Xi[i + 1, ]), t(mu1[, i + 1] + a[1, ] * (Xi[i, ] - mu0[, i])), Sigma_cond, TRUE),
          dmvnrm_arma_fast(t(Xi[i + 1, ]), t(mu0[, i + 1] + a[1, ] * (Xi[i, ] - mu1[, i])), Sigma_cond, TRUE),
          dmvnrm_arma_fast(t(Xi[i + 1, ]), t(mu1[, i + 1] + a[2, ] * (Xi[i, ] - mu1[, i])), Sigma_cond, TRUE)
        ), 2, 2)
        probs <- probs * exp(loglike - max(loglike))

        bkw[, i] <- probs %*% c(1, 1)
        bkw[, i] <- bkw[, i] / sum(bkw[, i])
        b_next <- bkw[, i]
      }
    }

    surv <- survival * bkw[1, ] / colSums(fwd * bkw)
    1 - surv
  }))

  cumprob
}


EM <- function(train, observedPats, test = NULL, maxIt = 5, r = 0.8, tol = 0.01, Estep = Estep_partial, nX = 10, survival=FALSE) {
  observedIndices <- which(train$ID %in% observedPats)
  unobservedIndices <- setdiff(seq(nrow(train)), observedIndices)
  train$pY <- train$pInv <- 1
  supPred <- semisupPred <- prediction <- aucs <- NULL
  lastY <- train$Y[unobservedIndices]

  trained_sup <- trained_semisup <- Mstep(train[observedIndices, ], r = r, nX = nX, survival = survival)
  if (!is.null(test)) {
    supPred <- Estep(test, trained_sup, nX = nX)
    suppressMessages({aucs <- c(pROC::auc(test$Y, supPred), rep(0, maxIt))})
  }

  for (it in 1:maxIt) {

    # E-step
    prediction <- Estep(train[unobservedIndices, ], trained_semisup, nX = nX)
    transProbs <- attr(prediction, "post2")

    train_augmented <- rbind(train[observedIndices, ], train[unobservedIndices, ], train[unobservedIndices, ])
    train_augmented$pY <- c(rep(1, length(observedIndices)), 1 - prediction, prediction)
    train_augmented$pInv <- 1 / train_augmented$pY
    train_augmented$Y <- c(train$Y[observedIndices], rep(0, length(unobservedIndices)), rep(1, length(unobservedIndices)))
    train_augmented$ID <- c(train$ID[observedIndices], train$ID[unobservedIndices], paste0(train$ID[unobservedIndices], ".2"))
    train_augmented <- train_augmented[train_augmented$pInv < 1e8, ]

    train_interTemp <- rbind(
      train[observedIndices[-1], ], train[unobservedIndices[-1], ], train[unobservedIndices[-1], ],
      train[unobservedIndices[-1], ], train[unobservedIndices[-1], ]
    )
    train_interTemp$Ycurr <- c(
      train$Y[observedIndices[-1]], rep(0, length(unobservedIndices) - 1), rep(0, length(unobservedIndices) - 1),
      rep(1, length(unobservedIndices) - 1), rep(1, length(unobservedIndices) - 1)
    )
    train_interTemp$Yprev <- c(
      train$Y[observedIndices[-length(observedIndices)]], rep(0, length(unobservedIndices) - 1),
      rep(1, length(unobservedIndices) - 1), rep(0, length(unobservedIndices) - 1), rep(1, length(unobservedIndices) - 1)
    )
    train_interTemp$pY <- c(
      as.integer(train$ID[observedIndices[-1]] == train$ID[observedIndices[-length(observedIndices)]]),
      transProbs[1, 1, -1], transProbs[2, 1, -1], transProbs[1, 2, -1], transProbs[2, 2, -1]
    )
    train_interTemp$Tprev <- c(
      train$T[observedIndices[-length(observedIndices)]],
      rep(train$T[unobservedIndices[-length(unobservedIndices)]], 4)
    )
    
    if (survival){
      train_interTemp$pY[train_interTemp$Yprev==1 & train_interTemp$Ycurr==1] <- 1
      train_interTemp$pY[train_interTemp$Yprev==1 & train_interTemp$Ycurr==0] <- 0
    }

    # M-step
    trained_old <- trained_semisup
    trained_semisup <- Mstep(train_augmented, train_interTemp, r = r, nX = nX, survival = survival)

    if (!is.null(test)) {
      oldPred <- semisupPred
      semisupPred <- Estep(test, trained_semisup, nX = nX)
      suppressMessages({aucs[it + 1] <- pROC::auc(test$Y, semisupPred)})
      
      if (it > 1 && aucs[it+1] < aucs[it]){
        trained_semisup <- trained_old
        semisupPred <- oldPred
        break
      }
    }
  
    if (all(abs(prediction - lastY) < tol)) {
      break
    }
    lastY <- prediction
  }

  return(list(
    "fitted_semisup" = trained_semisup, "fitted_sup" = trained_sup,
    "supPrediction" = supPred, "semisupPrediction" = semisupPred, "aucs" = aucs
  ))
}


lineSearch <- function(train, observedPats, test = NULL, nCrosses = 5, alphas = seq(0, 1, .1),
                       r = 0.8, Estep = Estep_partial, nX = 10, survival = FALSE) {
  if (length(alphas) == 1) {
    alpha <- alphas
  }
  else {
    n <- length(observedPats)
    observedPats <- sample(observedPats)
    valPats_overall <- lapply(1:nCrosses, function(i) {
      observedPats[round(n * (i - 1) / nCrosses + 1):round(n * i / nCrosses)]
    })

    alpha_results <- as.matrix(foreach(
      i = 1:nCrosses, .combine = cbind #, .packages = c("pROC", "nlme"), .noexport="dmvnrm_arma_fast",
      #.export = c("EM", "Estep", "Mstep", "fitGLS", "cumulative_prob", "abind")
    ) %do% {
      validatePats <- valPats_overall[[i]]
      trainPats <- setdiff(observedPats, validatePats)
      validation <- train[train$ID %in% validatePats,]

      tryCatch(
        {
          em <- EM(train, trainPats, validation, r = r, Estep = Estep, nX = nX, survival = survival)
          supervised <- Estep(validation, em$fitted_sup, nX = nX)
          semisupervised <- Estep(validation, em$fitted_semisup, nX = nX)

          sapply(alphas, function(alpha) {
            mixture <- alpha * semisupervised + (1 - alpha) * supervised
            suppressMessages({pROC::auc(validation$Y, mixture)})
          })
        },
        error = function(e) {
          rep(NA, length(alphas))
        }
      )
    })

    alpha <- alphas[which.max(rowMeans(alpha_results, na.rm = T))]
  }

  if (!is.null(test)) {
    em <- EM(train, observedPats, test, r = r, Estep = Estep, nX = nX, survival=survival)
    supervised <- em$supPrediction # Estep(test, em$fitted_sup, nX = nX)
    semisupervised <- em$semisupPrediction # Estep(test, em$fitted_semisup, nX = nX)
    mixture <- alpha * semisupervised + (1 - alpha) * supervised

    cumSup <- cumulative_prob(test, em$fitted_sup, nX = nX)
    cumSemisup <- cumulative_prob(test, em$fitted_semisup, nX = nX)
    cumMixture <- alpha * cumSemisup + (1 - alpha) * cumSup
  }
  else {
    supervised <- semisupervised <- mixture <- cumSup <- cumSemisup <- cumMixture <- NULL
  }

  return(list(
    "alpha" = alpha, "prediction" = mixture,
    "margSup" = supervised, "margSemisup" = semisupervised, "margMix" = mixture,
    "cumSup" = cumSup, "cumSemisup" = cumSemisup, "cumMix" = cumMixture
  ))
}


cv.r <- function(train, observedPats, nCrosses = 5, rs = seq(0, 1, .1), Estep = Estep_partial, nX = 10, survival=FALSE) {
  n <- length(observedPats)
  observedPats <- sample(observedPats)
  valPats_overall <- lapply(1:nCrosses, function(i) {
    observedPats[round(n * (i - 1) / nCrosses + 1):round(n * i / nCrosses)]
  })

  suppressWarnings({
    grid <- foreach(
      i = 1:nCrosses, .combine = cbind#, .export = c("Mstep", "Estep", "fitGLS", "abind"),
      #.packages = c("pROC", "nlme", "foreach", "parallel", "doParallel"), .noexport="dmvnrm_arma_fast"
    ) %do% {
      validatePats <- valPats_overall[[i]]
      trainPats <- setdiff(observedPats, validatePats)

      foreach(r = rs, .combine = c#, .export = c("Mstep", "Estep", "fitGLS", "abind"),
              #.packages = c("pROC", "nlme"), .noexport="dmvnrm_arma_fast"
      ) %do% {
        fitted_M <- Mstep(train[train$ID %in% trainPats, ], r = r, nX = nX, survival = survival)
        supervised <- Estep(train[train$ID %in% validatePats, ], fitted_M, nX = nX)
        suppressMessages({pROC::auc(train$Y[train$ID %in% validatePats], supervised)})
      }
    }
  })
  means <- rowMeans(grid)

  return(list("results" = grid, "r_opt" = rs[which.max(means)]))
}


objective_w <- function(w, args, lambda = 0) {
  C <- args$C
  y <- args$y
  V <- args$V
  X <- as.matrix(C) %*% (w * V)
  X0 <- X[y == 0, ]
  X1 <- X[y == 1, ]
  mu0 <- colMeans(X0)
  mu1 <- colMeans(X1)
  epsilon0 <- t(X0) - mu0
  epsilon1 <- t(X1) - mu1
  Sigma <- 1 / (nrow(X) - 1) * (epsilon0 %*% t(epsilon0) + epsilon1 %*% t(epsilon1) + 1e-6 * diag(length(mu0)))

  return(c(-1 * (mu1 - mu0) %*% solve(Sigma) %*% (mu1 - mu0) + lambda * sum(abs(w))))
}


numericGradientDescent <- function(x0, f, args = NULL, constIndex = 1, alphas = c(1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, .01, .02, .05, .1, .2, .5, 1),
                                   lambda = 0, maxIt = 100, tol = 1e-4) {
  for (it in 1:maxIt) {
    grad <- nloptr::nl.grad(x0, f, heps = .Machine$double.eps^(1 / 3), args, lambda)
    if (any(x0 == 0)) {
      grad[x0 == 0] <- sign(grad[x0 == 0]) * sapply(abs(grad[x0 == 0]) - lambda, function(xi) {
        max(xi, 0)
      })
    }
    grad[constIndex] <- 0
    grad <- grad / c(sqrt(grad %*% grad))

    alphaDeltas <- sapply(alphas, function(alpha) {
      deltas <- alpha * grad * sqrt(length(x0))
      x1 <- sapply(1:length(x0), function(i) {
        if (x0[i] > 0) {
          max(x0[i] - deltas[i], 0)
        }
        else if (x0[i] < 0) {
          min(x0[i] - deltas[i], 0)
        }
        else {
          x0[i] - deltas[i]
        }
      })
      f(x1, args, lambda)
    })

    if (all(alphaDeltas >= c(f(x0, args, lambda) - tol))) {
      break
    }
    else {
      alpha <- alphas[which.min(alphaDeltas)]
      deltas <- alpha * grad * sqrt(length(x0))
      x0 <- sapply(1:length(x0), function(i) {
        if (x0[i] > 0) {
          max(x0[i] - deltas[i], 0)
        }
        else if (x0[i] < 0) {
          min(x0[i] - deltas[i], 0)
        }
        else {
          x0[i] - deltas[i]
        }
      })
    }
  }

  return(x0)
}


cv.lambda <- function(C, y, V, w0 = NULL, nCrosses = 5, lambdas = NULL, surrIndex = 1) {
  if (is.null(w0)) {
    w0 <- glm(y ~ C, family = "quasibinomial")$coefficients[-1]
    w0[is.na(w0)] <- 0
  }
  if (is.null(lambdas)) {
    lambdas <- 10^seq(-3, -0.5, 0.5)
  }

  n <- nrow(C)
  observedPats <- sample(n)
  valPats_overall <- lapply(1:nCrosses, function(i) {
    observedPats[round(n * (i - 1) / nCrosses + 1):round(n * i / nCrosses)]
  })

  suppressWarnings({
    grid <- foreach(
      i = 1:nCrosses, .combine = cbind#, .export = c("numericGradientDescent", "objective_w"),
      #.packages = c("pROC", "foreach", "parallel", "doParallel"), .noexport="dmvnrm_arma_fast"
    ) %do% {
      validatePats <- valPats_overall[[i]]
      trainPats <- setdiff(observedPats, validatePats)

      foreach(lambda = lambdas, .combine = c) %do%{ #, .export = c("numericGradientDescent", "objective_w"), .noexport="dmvnrm_arma_fast") %do% {
        w_opt <- numericGradientDescent(w0, objective_w,
          args = list("C" = C[trainPats, ], "y" = y[trainPats], "V" = V),
          constIndex = surrIndex, lambda = lambda, maxIt = 50, tol = 1e-4
        )
        objective_w(w_opt, args = list("C" = C[validatePats, ], "y" = y[validatePats], "V" = V), lambda = 0)
      }
    }
  })
  means <- rowMeans(grid)

  lambda_opt <- lambdas[which.min(means)]

  return(list("results" = grid, "lambda_opt" = lambda_opt))
}


#' Semi-supervised Adaptive Markov Gaussian Process (SAMGEP)
#' 
#' @param dat_train (optional if Xtrain is supplied) Raw training data set, including patient IDs (ID), healthcare utilization feature (H) and censoring time (C)
#' @param dat_test (optional) Raw testing data set, including patient IDs (ID), a healthcare utilization feature (H) and censoring time (C)
#' @param Cindices (optional if Xtrain is supplied) Column indices of EHR feature counts in dat_train/dat_test
#' @param w (optional if Xtrain is supplied) Pre-optimized EHR feature weights
#' @param w0 (optional if Xtrain is supplied) Initial (i.e. partially optimized) EHR feature weights
#' @param V (optional if Xtrain is supplied) nFeatures x nEmbeddings embeddings matrix
#' @param observed (optional if Xtrain is supplied) IDs of patients with observed outcome labels
#' @param nX Number of embedding features (defaults to 10)
#' @param Estep E-step function to use (Estep_partial or Estep_full; defaults to Estep_partial)
#' @param Xtrain (optional) Embedded training data set, including patient IDs (ID), healthcare utilization feature (H) and censoring time (C)
#' @param Xtest (optional) Embedded testing data set, including patient IDs (ID), healthcare utilization feature (H) and censoring time (C)
#' @param alpha (optional) Relative weight of semi-supervised to supervised MGP predictors in SAMGEP ensemble
#' @param r (optional) Scaling factor of inter-temporal correlation
#' @param lambda (optional) L1 regularization hyperparameter for feature weight (w) optimization
#' @param surrIndex (optional) Index (within Cindices) of primary surrogate index for outcome event
#' @param nCores Number of cores to use for parallelization (defaults to 1)
#' 
#' @return w_opt Optimized feature weights (w)
#' @return r_opt Optimized inter-temporal correlation scaling factor (r)
#' @return alpha_opt Optimized semi-supservised:supervised relative weight (alpha)
#' @return lambda_opt Optiized L1 regularization hyperparameter (lambda)
#' @return margSup Posterior probability predictions of supervised model (MGP Supervised)
#' @return margSemisup Posterior probability predictions of semi-supervised model (MGP Semi-supervised)
#' @return margMix Posterior probability predictions of SAMGEP
#' @return cumSup Cumulative probability predictions of supervised model (MGP Supervised)
#' @return cumSemisup Cumulative probability predictions of semi-supervised model (MGP Semi-supervised)
#' @return cumMix Cumulative probability predictions of SAMGEP
#'
#' @export
samgep <- function(dat_train = NULL, dat_test = NULL, Cindices = NULL, w = NULL, w0 = NULL, V = NULL, observed = NULL, nX = 10, covs=NULL, survival=FALSE,
                   Estep = Estep_partial, Xtrain = NULL, Xtest = NULL, alpha = NULL, r = NULL, lambda = NULL, surrIndex = NULL, nCores = 1) {
  if (is.null(observed)) {
    observed <- unique(dat_train$ID)
  }

  if (nCores > 1) {
    logfile <- "SAMGEP.log"
    writeLines(c(""), file(logfile, "w"))
    clust <- parallel::makeCluster(nCores, outfile = logfile)
    doParallel::registerDoParallel(clust)
  }

  if (is.null(Xtrain)) {
    Ctrain <- as.matrix(dat_train[, Cindices])
    if (!is.null(dat_test)) {
      Ctest <- as.matrix(dat_test[, Cindices])
    }
    observedIndices <- which(dat_train$ID %in% observed)

    # Optimize w
    if (is.null(w)) {
      message("Fitting feature weights")
      if (is.null(w0)) {
        w0 <- glm(dat_train$Y[observedIndices] ~ Ctrain[observedIndices, ], family = "quasibinomial")$coefficients[-1]
        w0[is.na(w0)] <- 0
        if (any(abs(w0) > 1000)){
          w0 <- rep(1,ncol(Ctrain))
          if (is.null(surrIndex)){surrIndex <- 1}
        }
        else if (is.null(surrIndex)){
          surrIndex <- which.max(w0)
          w0 <- w0 / abs(w0[surrIndex])
        }
      }
      w0 <- numericGradientDescent(w0, objective_w,
        args = list("C" = Ctrain[observedIndices, ], "y" = dat_train$Y[observedIndices], "V" = V),
        constIndex = surrIndex, lambda = 0, maxIt = 100, tol = 1e-3
      )

      # Optimize lambda
      if (is.null(lambda)) {
        message("Cross-validating lambda")
        lambda <- cv.lambda(Ctrain[observedIndices, ], dat_train$Y[observedIndices], V, w0, surrIndex = surrIndex)$lambda_opt
        save(lambda,file='LC_lambdaopt.RData')
      }

      w <- numericGradientDescent(w0, objective_w,
        args = list("C" = Ctrain[observedIndices, ], "y" = dat_train$Y[observedIndices], "V" = V),
        constIndex = surrIndex, lambda = lambda, maxIt = 200, tol = 1e-4
      )
    }

    # Define Xtrain, Xtest
    CWVtrain <- Ctrain %*% (w * V)
    nonZeros <- which(colSums(CWVtrain)!=0)
    CWVtrain <- CWVtrain[,nonZeros]
    nX <- ncol(CWVtrain)
    Xtrain <- data.frame(ID = dat_train$ID, Y = dat_train$Y, T = dat_train$T, Tlog = log(dat_train$T + 1), H = dat_train$H, Hlog = log(dat_train$H + 1))
    Xtrain <- cbind(Xtrain, CWVtrain)
    colnames(Xtrain)[-c(1:6)] <- paste0("X", seq(nX))
    Xtrain$pY <- Xtrain$pInv <- 1
    if (!is.null(dat_test)) {
      CWVtest <- Ctest %*% (w * V)
      CWVtest <- CWVtest[,nonZeros]
      Xtest <- data.frame(ID = dat_test$ID, Y = dat_test$Y, T = dat_test$T, Tlog = log(dat_test$T + 1), H = dat_test$H, Hlog = log(dat_test$H + 1))
      Xtest <- cbind(Xtest, CWVtest)
      colnames(Xtest)[-c(1:6)] <- paste0("X", seq(nX))
    }
  }

  # Optimize r
  if (is.null(r)) {
    message("Cross-validating r")
    r <- cv.r(Xtrain, observed, Estep = Estep, nX = nX, survival = survival)$r_opt
  }

  # Optimize alpha and predict
  message("Fitting MGP")
  if (is.null(alpha)) {
    result <- lineSearch(Xtrain, observed, Xtest, r = r, Estep = Estep, nX = nX, survival=survival)
  }
  else {
    result <- lineSearch(Xtrain, observed, Xtest, alphas = alpha, r = r, Estep = Estep, nX = nX, survival=survival)
  }
  alpha <- result$alpha
  
  if (nCores > 1){
    parallel::stopCluster(clust) 
  }

  # Result
  return(list(
    "w_opt" = w, "r_opt" = r, "alpha_opt" = alpha, "lambda_opt" = lambda,
    "margSup" = result$margSup, "margSemisup" = result$margSemisup, "margMix" = result$margMix,
    "cumSup" = result$cumSup, "cumSemisup" = result$cumSemisup, "cumMix" = result$cumMix
  ))
}
