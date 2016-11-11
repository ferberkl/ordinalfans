#### Ordinal FANS: Approach #1 ############

### Model Fitting

binFANS <- function(x, ...) UseMethod("binFANS")  # Generic function

#Function for running ordinal FANS algorithm
binFANS.default <- function(x, y, newx=NULL, niter=1, seed=2468, 
                            scale=TRUE, parallel=FALSE) { 

  if (class(y)[1] != "ordered") {
    stop("y must be an ordered factor.")
  }
  orig.y <- y
  y <- NULL
  # For each feature, add a small amount of noise to the duplicate values.
  # Avoids issues estimating densities.
  orig.x <- apply(x, 2, function(d) {
    set.seed(seed)
    d[duplicated(d)] <- jitter(d[duplicated(d)], factor=0.5)
    return(d)
  })
  if (scale) {
    trans <- preProcess(as.data.frame(orig.x), method=c("center", "scale")) 
    if (!is.null(newx)) newx <- predict(trans, as.data.frame(newx))
    orig.x <- predict(trans, as.data.frame(orig.x))
  } else {
    trans <-NULL
  }
  x <- NULL
  K <- length(unique(orig.y))
  n <- dim(orig.x)[1]
  levels <- unique(orig.y)
  set.seed(seed)
  partitions  <-  createDataPartition(orig.y, p=0.5, 
                                      times=ceiling(niter / 2)) 
  call <- match.call()
  matrix.sum <- function(mat1, mat2) {
    return(as.matrix(mat1) + as.matrix(mat2))
  }
  # Part of the algorithm that is repeated niter times
  if (parallel) {
    `%fun%` <- `%dopar%`
  } else {
    `%fun%` <- `%do%`    
  }
  fitted.model <- foreach(i=1:niter, .combine='rbind', 
                          .packages=c('sm', 'caret')) %fun% { 
    #Create the data partitions, i.e. (D1, D1^c), (D2, D2^c), ..., (DL, DL^c) 
    if(i%%2!=0) { 
      D.train.x <- orig.x[partitions[[(i - 1) / 2 + 1]], ]  #x in D1, D3, ...
      D.train.y <- orig.y[partitions[[(i - 1) / 2 + 1]]]  #y in D1, D3, ...
      D.test.x <- orig.x[-partitions[[(i - 1) / 2 + 1]], ]  #x in D1^c, D3^c, ...
      D.test.y <- orig.y[-partitions[[(i - 1) / 2 + 1]]]  #y in D1^c, D3^c, ...
    } else {
      D.train.x <- orig.x[-partitions[[i / 2]], ] #x in D2, D4, D6, ...
      D.train.y <- orig.y[-partitions[[i / 2]]] #y in D2, D4, D6, ...
      D.test.x <- orig.x[partitions[[i / 2]], ] #x in D2^c, D4^c, D6^c, ...
      D.test.y <- orig.y[partitions[[i / 2]]]   #y in D2^c, D4^c, D6^c, ...
    }
    D.train.y.bin <- matrix(nrow=length(D.train.y), ncol=K - 1)
    D.test.y.bin <- matrix(nrow=length(D.test.y), ncol=K - 1)
    for(k in 1:(K - 1)) {
      D.train.y.bin[, k] <- ifelse(as.numeric(D.train.y) <= k, 1, 0)
      D.test.y.bin[, k] <- ifelse(as.numeric(D.test.y) <= k, 1, 0)
    }
    f.marginals <- array(dim=c(dim(D.test.x)[1], dim(orig.x)[2], K - 1)) 
    g.marginals <- array(dim=c(dim(D.test.x)[1], dim(orig.x)[2], K - 1)) 
    #if newx is specified, fit model and predict outcome of new observations
    if(!is.null(newx)) {
      if(class(newx) == "numeric") {newx <- t(as.matrix(newx))}
      test.marginals.f <- array(dim=c(dim(newx)[1], dim(orig.x)[2], K - 1))
      test.marginals.g <- array(dim=c(dim(newx)[1], dim(orig.x)[2], K - 1))
      for(k in 1:(K - 1)) {
        for(j in 1:dim(orig.x)[2]) {
          f.estimates <- sm.density(D.train.x[as.numeric(D.train.y) <= k, j], 
                                    eval.points=c(D.test.x[, j], newx[, j]), 
                                    display="none")$estimate
          g.estimates <- sm.density(D.train.x[as.numeric(D.train.y) > k, j], 
                                    eval.points=c(D.test.x[, j], newx[, j]),
                                    display="none")$estimate
          #Marginal density estimates for observations in D^c
          f.marginals[, j, k] <- f.estimates[1:length(D.test.x[, j])]
          g.marginals[, j, k] <- g.estimates[1:length(D.test.x[, j])]
          #Marginal density estimates for newx
          test.marginals.f[, j, k] <- f.estimates[(length(D.test.x[, j]) + 1):
                                                  length(f.estimates)]
          test.marginals.g[, j, k] <- g.estimates[(length(D.test.x[, j]) + 1):
                                                  length(g.estimates)]
        }   
      }
      f.marginals[f.marginals < 0.001] <- 0.001
      g.marginals[g.marginals < 0.001] <- 0.001
      test.marginals.g[test.marginals.g < 0.001] <- 0.001
      test.marginals.f[test.marginals.f < 0.001] <- 0.001
      Z <- log(f.marginals / g.marginals) #array of augmented features. 
      Z.test <- log(test.marginals.f / test.marginals.g) 
      for (g in 1:dim(Z)[3]) {
        colnames(Z[, , g]) <- paste(colnames(orig.x), g, sep=".")
      }
      #Fit the K - 1 binary response models 
      p.mat <- array(dim=c(K, dim(Z.test)[1], K - 1))
      create.mat <- function(w) {
        c(rep(w, k), rep(1 - w, K - k))
      }
      coefs <- matrix(nrow=K - 1, ncol=dim(orig.x)[2] + 1)
      colnames(coefs) <- c("Intercept", colnames(orig.x))
      rownames(coefs) <- sapply(as.character(1:nrow(coefs)), 
                                function(x) {
                                  paste("Modeling logit[P(Y <= ", ")]", sep=x)
                                })
      fits <- list()
      min.AIC <- c()
      for(k in 1:(K - 1)) {
        fits[[k]] <- glmpath(x=Z[, , k], y=D.test.y.bin[, k],
                                       family=binomial)
        min.AIC[k] <- as.numeric(gsub("Step ", "", 
                                 rownames(summary(fits[[k]])
                                 [which.min(summary(fits[[k]])$AIC), ])))
        p.hat <- predict(fits[[k]], newx=Z.test[, , k],
                         s=min.AIC[k], type="response") #predict y<=k or y>k  
        coefs[k, ] <- predict(fits[[k]], s=min.AIC[k], type="coefficients") 
        p.mat[, , k] <- sapply(p.hat, create.mat)
      }
      scores <- apply(p.mat, c(1, 2), sum) #aggregate binary predictions
      return(list(fits=fits, minAIC=min.AIC, scores=t(scores), coefs=coefs,
                  D.train.x=D.train.x, D.train.y=D.train.y, trans=trans))   
    #if newx is NULL, just fit K - 1 binary response models
    } else {
      for(k in 1:(K - 1)) {
        for(j in 1:dim(orig.x)[2]) {
          #Marginal density estimates for observations in D^c
          f.marginals[, j, k] <- sm.density(D.train.x[as.numeric(D.train.y) <= k, j], 
                                            eval.points=c(D.test.x[, j]), 
                                            display="none")$estimate
          g.marginals[, j, k] <- sm.density(D.train.x[as.numeric(D.train.y) > k, j], 
                                            eval.points=c(D.test.x[, j]), 
                                            display="none")$estimate
        }   
      }
      f.marginals[f.marginals < 0.001] <- 0.001 
      g.marginals[g.marginals < 0.001] <- 0.001
      Z <- log(f.marginals / g.marginals) 
      #Fit the K - 1 binary response models 
      fits <- list()
      coefs <- matrix(nrow=K - 1, ncol=dim(orig.x)[2] + 1)
      colnames(coefs) <- c("Intercept", colnames(orig.x))
      rownames(coefs) <- sapply(as.character(1:nrow(coefs)), 
                                function(x) {
                                  paste("Modeling logit[P(Y <= ", ")]", sep=x)
                                })
      min.AIC <- c()
      for(k in 1:(K - 1)) {
        fits[[k]] <- glmpath(x=Z[, , k], y=D.test.y.bin[, k], family=binomial)
        min.AIC[k] <- as.numeric(gsub("Step ", "", 
                                 rownames(summary(fits[[k]])
                                 [which.min(summary(fits[[k]])$AIC), ])))
        coefs[k, ] <- predict(fits[[k]], s=min.AIC[k], type="coefficients") 
      }
      return(list(fits=fits, x=orig.x, D.train.x=D.train.x, coefs=coefs,
                  D.train.y=D.train.y, K=K, minAIC=min.AIC, niter=niter, 
                  call=call, trans=trans))
    }
  }
  if (niter == 1) {
    if(!is.null(newx)) {
      output <- list(coefs=fitted.model$coefs, niter=niter, K=K, x=orig.x,
                     fits=fitted.model$fits, minAIC=fitted.model$minAIC,
                     D.train.x=fitted.model$D.train.x, y=orig.y,
                     D.train.y=fitted.model$D.train.y, call=call,
                     trans=trans, newx=newx, scores=fitted.model$scores, 
                     pred.class=ordered(unique(orig.y)[apply(fitted.model$scores,
                                                       1, which.max)],
                                        levels=levels(orig.y)))
    } else {
      output <- list(coefs=fitted.model$coefs, niter=niter, K=K, x=orig.x,
                     fits=fitted.model$fits, minAIC=fitted.model$minAIC,
                     D.train.x=fitted.model$D.train.x, y=orig.y,
                     D.train.y=fitted.model$D.train.y, call=call,
                     trans=trans)
    }
  } else {
    if(!is.null(newx)) {
      fits <- lapply(fitted.model[, 1], as.matrix)
      minAIC <- fitted.model[, 2]
      scores <- lapply(fitted.model[, 3], as.matrix)
      agg.scores <- Reduce('+',lapply(fitted.model[, 3], as.matrix))
      coefs <- lapply(fitted.model[, 4], as.matrix)
      D.train.x <- lapply(fitted.model[, 5], as.matrix)
      D.train.y <- lapply(fitted.model[, 6], unlist)
      pred.class <- ordered(unique(orig.y)[apply(agg.scores, 1, which.max)],
                            levels=levels(orig.y))
      output <- list(coefs=coefs, niter=niter, K=K, x=orig.x, fits=fits, 
                  minAIC=minAIC, D.train.x=D.train.x, y=orig.y, 
                  D.train.y=D.train.y, call=call, trans=trans, newx=newx, 
                  scores=scores, agg.scores=agg.scores, pred.class=pred.class)
    } else {
      fits <- lapply(fitted.model[, 1], as.matrix)
      D.train.x <- lapply(fitted.model[, 3], as.matrix)
      D.train.y <- lapply(fitted.model[, 5], unlist)
      minAIC <- fitted.model[, 7]
      coefs <- lapply(fitted.model[, 4], as.matrix)
      output <- list(coefs=coefs, niter=niter, K=K, x=orig.x, fits=fits, 
                     minAIC=minAIC, D.train.x=D.train.x, y=orig.y, 
                     D.train.y=D.train.y, call=call, trans=trans)
    }
  }
  class(output) <- "binFANS"
  return(output)
}
