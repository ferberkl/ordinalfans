#### Ordinal FANS: Approach #2 #######################################

### Model Fitting

sum.across.L <- function(Li, Lj) {
  # Sums the results of the L fitted models
  #
  # Args:
  #    Li: list of results of a fitted model
  #    Lj: list of results of another fitted model
  #
  # Returns:
  #   A list composed of the sums of the elements of the L
  #      lists produced from the L fitted models.
  tmp<-rbind(Li, Lj)
  apply(tmp, 2, function(x) Reduce('+', x))
}

ordinalFANS <- function(x, ...) UseMethod("ordinalFANS")  # Generic function

ordinalFANS.default<-function(x, y, newx=NULL, eps=0.1, niter=1, 
                              seed=2468, mstop=100, scale=TRUE,
                              parallel=FALSE) {
  # Runs the ordinal FANS algorithm
  #
  # Args:
  #   x: the n x p design matrix
  #   y: the response vector
  #   newx: (optional) design matrix with which to predict outcomes
  #   niter: number of times to repeat the algorithm with new data partitions
  #   seed: sets the seed for the initial partitioning of the data
  #   mstop: Stopping iteration
  #
  call <- match.call()
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
  # Stratified random sampling
  partitions <- createDataPartition(orig.y, p=0.5, 
                                    times=ceiling(niter / 2))

  # Part of the algorithm that is repeated niter times
  if (parallel) {
    `%fun%` <- `%dopar%`
  } else {
    `%fun%` <- `%do%`    
  }
  if (niter > 1) {
    fitted.model <- foreach(i=1:niter, .combine='rbind', 
                            .packages=c('sm', 'caret')) %fun% {
      # Create the data partitions, (D1, D1^c), (D2, D2^c), ..., (DL, DL^c)
      if (i%%2 != 0) {
        D.train.x <- orig.x[partitions[[(i - 1) / 2 + 1]], ]  
        D.train.y <- orig.y[partitions[[(i - 1) /  2 + 1]]]   
        D.test.x <- orig.x[-partitions[[(i - 1) /  2 + 1]], ] 
        D.test.y <- orig.y[-partitions[[(i - 1) /  2 + 1]]]   
      } else {
        D.train.x <- orig.x[-partitions[[i / 2]], ]
        D.train.y <- orig.y[-partitions[[i / 2]]]  
        D.test.x <- orig.x[partitions[[i / 2]], ]  
        D.test.y <- orig.y[partitions[[i / 2]]]    
      }
      # Estimate marginals using data in Di, evaluate using data in Di^c,
      #   and calculate augmented features (log ratios)
      # If newx is supplied, fit model, then find predicted values for new obs
      if (!is.null(newx)) {
        f.estimates <- array(dim=c(dim(D.test.x)[1] + dim(newx)[1],
                                 K - 1, dim(orig.x)[2]))
        g.estimates <- array(dim=c(dim(D.test.x)[1] + dim(newx)[1],
                                 K - 1, dim(orig.x)[2]))
        for (k in 1:(K - 1)) {
          for (j in 1:dim(orig.x)[2]) {
            f.estimates[, k, j] <-
              sm.density(D.train.x[which(as.numeric(D.train.y) <= k), j],
                         eval.points=c(D.test.x[, j], newx[, j]),
                         display="none")$estimate
            g.estimates[, k, j] <-
              sm.density(D.train.x[which(as.numeric(D.train.y) > k), j],
                         eval.points=c(D.test.x[, j], newx[, j]),
                         display="none")$estimate
          }
        }
        # Marginal density estimates for observations in D^c
        f.marginals <- f.estimates[1:dim(D.test.x)[1], , ]
        g.marginals <- g.estimates[1:dim(D.test.x)[1], , ]
        # Marginal density estimates for newx
        test.marginals.f <- f.estimates[(dim(D.test.x)[1] + 1):
                                        dim(f.estimates)[1], , ]
        test.marginals.g <- g.estimates[(dim(D.test.x)[1] + 1):
                                        dim(g.estimates)[1], , ]

        # Winsorization to improve stability of estimates as 
        # suggested in FANS manuscript
        f.marginals[f.marginals < 0.001] <- 
          g.marginals[g.marginals < 0.001] <-
          test.marginals.g[test.marginals.g < 0.001] <-
          test.marginals.f[test.marginals.f < 0.001] <- 0.001
        # Array of augmented features for x. Within the array: 
        #   n x (K - 1) matrices of augmented features .
        #   array dimensions are n x (K - 1) x p
        Z <- log(f.marginals / g.marginals)
        full <- data.frame(Z)
        colnames(full) <- paste("z", paste(rep(1:dim(Z)[3], each=K - 1),
                                    rep(1:(K - 1), times=dim(Z)[3]), 
                                    sep="."), sep="")
        # Array of augmented features for newx. Within the array: 
        #   n x (K - 1) matrices of augmented features .
        #   array dimensions are n x (K - 1) x p
        Z.test <- log(test.marginals.f / test.marginals.g)
        full.test <- data.frame(Z.test)
        colnames(full.test) <- paste("z", paste(rep(1:dim(Z)[3], each=K - 1),
                                     rep(1:(K - 1), 
                                     times=dim(f.marginals)[3]),
                                     sep="."), sep="")
        orig.vars.index <- rep(1:dim(Z)[3], each=K - 1)
        # Fit proportional odds boosting model
        fit <- PO.boost(x=full, y=D.test.y, eps=eps, mstop=mstop)
        # Variable importance measure = absolute value of sum of the (K - 1) 
        #   augmented features coefficient estimates 
        var.importance <- apply(fit$coefs, 1, function(x) {
          tapply(x, orig.vars.index, function(y) {
            abs(sum(y))
          })
        })
        var.importance <- data.frame(t(var.importance))
        colnames(var.importance) <- colnames(orig.x)
        # function estimates for newx
        offset <- as.numeric(-1 * (matrix(apply(full, 2, mean), nrow=1) %*%
                  matrix(fit$coefs[mstop, ], ncol=1)))
        f.newx <- offset + (as.matrix(full.test) %*% fit$coefs[mstop, ])  
        # estimated posterior probabilities for newx
        post.probs <- response(f.newx, fit$theta[mstop, ])  
        risk <- apply(fit$loss, 1, sum)
        return(list(coefs=fit$coefs, theta=fit$theta, risk=risk,
                    post.probs=post.probs, D.train.x=D.train.x, 
                    D.train.y=D.train.y, trans=trans, 
                    var.importance=var.importance,
                    aug.newX=as.matrix(full.test), aug.X=as.matrix(full)))


        # If newx is NOT supplied 
      } else {
        # Initialize array of f marginals
        f.marginals <- array(dim=c(dim(D.test.x)[1], 
                                   K - 1, dim(orig.x)[2]))
        # Initialize array of g marginals
        g.marginals <- array(dim=c(dim(D.test.x)[1],
                                   K - 1, dim(orig.x)[2]))
        for (k in 1:(K - 1)) {
          for (j in 1:dim(orig.x)[2]) {
            # Marginal density estimates for observations in D^c
            f.marginals[, k, j] <-
              sm.density(D.train.x[which(as.numeric(D.train.y) <= k), j],
                         eval.points=D.test.x[, j],
                         display="none")$estimate
            g.marginals[, k, j] <-
              sm.density(D.train.x[which(as.numeric(D.train.y) > k), j],
                         eval.points=D.test.x[, j],
                         display="none")$estimate
          }
        }
        # Winsorization to improve stability of estimates 
        # as suggested in FANS
        f.marginals[f.marginals < 0.001] <- 
        g.marginals[g.marginals < 0.001] <- 0.001
        # Array of augmented features for x. Within the array: 
        #   n x (K - 1) matrices of augmented features.
        # Array dimensions are n x (K - 1) x p
        Z <- log(f.marginals / g.marginals)
        full <- data.frame(Z)
        colnames(full) <- paste("z",
                                paste(rep(1:dim(Z)[3], each=K - 1),
                                      rep(1:(K - 1), times=dim(Z)[3]),
                                      sep="."), sep="")
        orig.vars.index <- rep(1:dim(Z)[3], each=K - 1)
        # Fit proportional odds boosting model
        fit <- PO.boost(x=full, y=D.test.y, mstop=mstop)
        var.importance <- apply(fit$coefs, 1, function(x) {
          tapply(x, orig.vars.index, function(y) {
            abs(sum(y))
          })
        })
        var.importance <- data.frame(t(var.importance))
        colnames(var.importance) <- colnames(orig.x)
        risk <- apply(fit$loss, 1, sum)
        return(list(coefs=fit$coefs, theta=fit$theta, risk=risk,
                    D.train.x=D.train.x, D.train.y=D.train.y, trans=trans,
                    var.importance=var.importance, aug.X=as.matrix(full)))
      }
    }
    if (!is.null(newx))  {
      coefs <- lapply(fitted.model[, 1], as.matrix)
      mean.coefs <- Reduce('+', fitted.model[, 1]) / niter
      theta <- lapply(fitted.model[, 2], as.matrix)
      mean.theta <- Reduce('+', fitted.model[, 2]) / niter
      risk <- lapply(fitted.model[, 3], as.numeric)
      mean.risk <- as.numeric(Reduce('+', risk) / niter)
      if (dim(newx)[1] == 1) {
        posteriors <- lapply(fitted.model[, 4], function(q) t(as.matrix(q)))
      } else {
        posteriors <- lapply(fitted.model[, 4], as.matrix)
      }
      mean.posteriors <- Reduce('+', posteriors) / niter
      D.train.x <- lapply(fitted.model[, 5], as.matrix)
      D.train.y <- lapply(fitted.model[, 6], unlist)
      pred.class <-
        factor(levels(orig.y)[apply(mean.posteriors, 1, which.max)],
               levels=levels(orig.y), ordered=TRUE)
      var.importance <- lapply(fitted.model[, 8], data.frame)
      mean.var.importance <-
        data.frame(Reduce('+', fitted.model[, 8]) / niter)
      aug.newX <- lapply(fitted.model[, 9], as.matrix)
      aug.X <- lapply(fitted.model[, 10], as.matrix)
      output <- list(call=call, eps=eps, mstop=mstop,
                     niter=niter, K=K, coefs=coefs,
                     mean.coefs=mean.coefs, theta=theta,
                     mean.theta=mean.theta,
                     risk=risk, mean.risk=mean.risk, y=orig.y,
                     D.train.x=D.train.x, D.train.y=D.train.y,
                     trans=trans, var.importance=var.importance,
                     aug.X=aug.X,
                     mean.var.importance=mean.var.importance, 
                     newx=newx, aug.newX=aug.newX,
                     posteriors=posteriors,
                     mean.posteriors=mean.posteriors,
                     pred.class=pred.class)
    }
    else {
      coefs <- lapply(fitted.model[, 1], as.matrix)
      mean.coefs <- Reduce('+', fitted.model[, 1]) / niter
      theta <- lapply(fitted.model[, 2], as.matrix)
      mean.theta <- Reduce('+', fitted.model[, 2]) / niter
      risk <- lapply(fitted.model[, 3], as.numeric)
      mean.risk <- as.numeric(Reduce('+', risk) / niter)
      D.train.x <- lapply(fitted.model[, 4], as.matrix)
      D.train.y <- lapply(fitted.model[, 5], unlist)
      var.importance <- lapply(fitted.model[, 7], data.frame)
      aug.X <- lapply(fitted.model[, 8], as.matrix)
      mean.var.importance <-
        data.frame(Reduce('+', fitted.model[, 7]) / niter)
      output <- list(call=call, eps=eps, mstop=mstop,
                     niter=niter, K=K, coefs=coefs,
                     mean.coefs=mean.coefs, theta=theta,
                     mean.theta=mean.theta,
                     risk=risk, mean.risk=mean.risk, y=orig.y,
                     D.train.x=D.train.x, D.train.y=D.train.y,
                     trans=trans, var.importance=var.importance,
                     aug.X=aug.X, mean.var.importance=mean.var.importance)
    }


  # if niter=1
  } else {
    # Create the data partitions, i.e. (D1, D1^c)
    D.train.x <- orig.x[partitions[[1]], ]  #x in D1
    D.train.y <- orig.y[partitions[[1]]]    #y in D1
    D.test.x <- orig.x[-partitions[[1]], ]  #x in D1^c
    D.test.y <- orig.y[-partitions[[1]]]    #y in D1^c
    # Estimate marginals using data in Di, evaluate using data in Di^c,
    #    and calculate augmented features (log ratios)
    # If newx is supplied, fit model, then find predicted values for new obs
    if (!is.null(newx)) {
      f.estimates <- array(dim=c(dim(D.test.x)[1] + dim(newx)[1],
                                 K - 1, dim(orig.x)[2]))
      g.estimates <- array(dim=c(dim(D.test.x)[1] + dim(newx)[1],
                                 K - 1, dim(orig.x)[2]))
      for (k in 1:(K - 1)) {
        for (j in 1:dim(orig.x)[2]) {
          f.estimates[, k, j] <-
            sm.density(D.train.x[which(as.numeric(D.train.y) <= k), j],
                       eval.points=c(D.test.x[, j], newx[, j]),
                       display="none")$estimate
          g.estimates[, k, j] <-
            sm.density(D.train.x[which(as.numeric(D.train.y) > k), j],
                       eval.points=c(D.test.x[, j], newx[, j]),
                       display="none")$estimate
        }
      }
      # Marginal density estimates for observations in D^c
      f.marginals <- f.estimates[1:dim(D.test.x)[1], , ]
      g.marginals <- g.estimates[1:dim(D.test.x)[1], , ]
      # Marginal density estimates for newx
      test.marginals.f <- f.estimates[(dim(D.test.x)[1] + 1):
                                      dim(f.estimates)[1], , ]
      test.marginals.g <- g.estimates[(dim(D.test.x)[1] + 1):
                                      dim(g.estimates)[1], , ]
      # Winsorization to improve stability of estimates 
      # as suggested in FANS manuscript
      f.marginals[f.marginals < 0.001] <-
        g.marginals[g.marginals < 0.001] <-
        test.marginals.g[test.marginals.g < 0.001] <-
        test.marginals.f[test.marginals.f < 0.001] <- 0.001
        # Array of augmented features for x. Within the array: 
        #   n x (K - 1) matrices of augmented features.
        # Array dimensions are n x (K - 1) x p
      Z <- log(f.marginals / g.marginals)
      full <- data.frame(Z)
      colnames(full) <- paste("z",
                              paste(rep(1:dim(Z)[3], each=K - 1),
                                    rep(1:(K - 1), times=dim(Z)[3]), 
                                    sep="."), sep="")
        # Array of augmented features for newx. Within the array: 
        #   n x (K - 1) matrices of augmented features .
        #   array dimensions are n x (K - 1) x p
      Z.test <- log(test.marginals.f / test.marginals.g)
      full.test <- data.frame(Z.test)
      colnames(full.test) <- paste("z",
                                   paste(rep(1:dim(Z)[3], each=K - 1),
                                   rep(1:(K - 1), 
                                   times=dim(f.marginals)[3]),
                                   sep="."), sep="")
      orig.vars.index <- rep(1:dim(Z)[3], each=K - 1)
      # Fit proportional odds boosting model
      fit <- PO.boost(x=full, y=D.test.y, eps=eps, mstop=mstop)
      # Variable importance measure = absolute value of sum of the (K - 1) 
      #   augmented features coefficient estimates 
      var.importance <- apply(fit$coefs, 1, function(x) {
        tapply(x, orig.vars.index, function(y) {
          abs(sum(y))
        })
      })
      var.importance <- data.frame(t(var.importance))
      colnames(var.importance) <- colnames(orig.x)
      # function estimates for newx
      offset <- as.numeric(-1 * (matrix(apply(full, 2, mean), nrow=1) %*%
                matrix(fit$coefs[mstop, ], ncol=1)))
      f.newx <- offset + (as.matrix(full.test) %*% fit$coefs[mstop, ])  
      # estimated posterior probabilities for newx 
      post.probs <- response(f.newx, fit$theta[mstop, ])  
      risk <- apply(fit$loss, 1, sum)


      # If newx is NOT supplied 
    } else {
      # Initialize array of f marginals
      f.marginals <- array(dim=c(dim(D.test.x)[1],
                           K - 1, dim(orig.x)[2]))
      # Initialize array of g marginals
      g.marginals <- array(dim=c(dim(D.test.x)[1],
                           K - 1, dim(orig.x)[2]))
      for (k in 1:(K - 1)) {
        for (j in 1:dim(orig.x)[2]) {
          # Marginal density estimates for observations in D^c
          f.marginals[, k, j] <-
            sm.density(D.train.x[which(as.numeric(D.train.y) <= k), j],
                                  eval.points=D.test.x[, j],
                                  display="none")$estimate
          g.marginals[, k, j] <-
            sm.density(D.train.x[which(as.numeric(D.train.y) > k), j],
                                  eval.points=D.test.x[, j],
                                  display="none")$estimate
        }
      }
      # Winsorization to improve stability of estimates
      #  as suggested in the FANS manuscript
      f.marginals[f.marginals < 0.001] <- 
      g.marginals[g.marginals < 0.001] <- 0.001
      # Array of augmented features for x. Within the array: 
      #   n x (K - 1) matrices of augmented features .
      #   array dimensions are n x (K - 1) x p
      Z <- log(f.marginals / g.marginals)
      full <- data.frame(Z)
      colnames(full) <- paste("z",
                              paste(rep(1:dim(Z)[3], each=K - 1),
                                    rep(1:(K - 1), times=dim(Z)[3]),
                                    sep="."), sep="")
      orig.vars.index <- rep(1:dim(Z)[3], each=K - 1)
      # Fit proportional odds boosting model
      fit <- PO.boost(x=full, y=D.test.y, mstop=mstop)
      var.importance <- apply(fit$coefs, 1, function(x) {
        tapply(x, orig.vars.index, function(y) {
          abs(sum(y))
        })
      })
      var.importance <- data.frame(t(var.importance))
      colnames(var.importance) <- colnames(orig.x)
      risk <- apply(fit$loss, 1, sum)
    }
    if (!is.null(newx))  {
      coefs <- as.matrix(fit$coefs)
      theta <- as.matrix(fit$theta)
      risk <- as.numeric(risk)
      if (dim(newx)[1] == 1) {
        posteriors <- t(as.matrix(post.probs))
      } else {
        posteriors <- as.matrix(post.probs)
      }
      pred.class <-
        factor(levels(orig.y)[apply(posteriors, 1, which.max)],
               levels=levels(orig.y), ordered=TRUE)
      D.train.x <- as.matrix(D.train.x)
      D.train.y <- as.matrix(D.train.y)
      var.importance <- data.frame(var.importance)
      aug.newX <- as.matrix(full.test)
      output <- list(call=call, eps=eps, mstop=mstop, niter=niter,
                     K=K, coefs=coefs, theta=theta, risk=risk,
                     y=orig.y, D.train.x=D.train.x,
                     D.train.y=D.train.y, trans=trans, 
                     var.importance=var.importance, aug.X=full,
                     newx=newx, aug.newX=aug.newX,
                     posteriors=posteriors, pred.class=pred.class)
    } else {
      coefs <- as.matrix(fit$coefs)
      theta <- as.matrix(fit$theta)
      risk <- as.numeric(risk)
      D.train.x <- as.matrix(D.train.x)
      D.train.y <- as.matrix(D.train.y)
      var.importance <- data.frame(var.importance)
      output <- list(call=call, eps=eps, mstop=mstop,
                     niter=niter, K=K, coefs=coefs, 
                     theta=theta, risk=risk, y=orig.y,
                     D.train.x=D.train.x, 
                     D.train.y=D.train.y, trans=trans, 
                     var.importance=var.importance, aug.X=full)
    }
  }
  class(output) <- "ordinalFANS"
  return(output)
}
