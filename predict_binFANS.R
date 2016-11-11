### Predict the outcome of a new observation

predict.binFANS <- function(object, newx) {
  newx <- predict(object$trans, newx)  # Scale newx with mean and sd of X
  if (!is.null(object$newx)) {
    if (!is.null(object$pred.class) & all.equal(as.matrix(object$newx),
                                                newx)[1]) {
      #if (type=="class") {
      #  return(object$pred.class)
      #} else {
      #  return(object$scores)
      #}
      if (object$niter == 1) {
        return(list(scores=object$scores, pred.class=object$pred.class))
      } else {
        return(list(scores=object$scores, agg.scores=object$agg.scores,
                    pred.class=object$pred.class))
      }
    }
  }
  K <- object$K
  x <- object$x
  niter <- object$niter
  predicted <- foreach(i=1:niter, .combine='rbind', 
                       .packages=c('sm', 'caret')) %dopar% {
    if (niter==1) {
      D.train.x <- object$D.train.x
      D.train.y <- object$D.train.y
      fits <- object$fits
      min.AIC <- object$minAIC
    } else {
      D.train.x <- object$D.train.x[[i]]
      D.train.y <- object$D.train.y[[i]]
      fits <- object$fits[[i]]
      min.AIC <- object$minAIC[[i]]
    }
    test.marginals.f <- array(dim=c(dim(newx)[1],
                              dim(newx)[2], K - 1))
    test.marginals.g <- array(dim=c(dim(newx)[1],
                              dim(newx)[2], K - 1))
    for(k in 1:(K - 1)) {
      for(j in 1:dim(x)[2]) {
          test.marginals.f[, j, k] <-
            sm.density(D.train.x[as.numeric(D.train.y) <= k, j],
                                            eval.points=c(newx[, j]), 
                                            display="none")$estimate
          test.marginals.g[, j, k] <-
            sm.density(D.train.x[as.numeric(D.train.y) > k, j],
                                            eval.points=c(newx[, j]),
                                            display="none")$estimate
      }   
    }
    test.marginals.g[test.marginals.g < 0.001] <- 0.001
    test.marginals.f[test.marginals.f < 0.001] <- 0.001
    Z.test <- log(test.marginals.f / test.marginals.g) 
    p.mat <- array(dim=c(K, dim(Z.test)[1], K - 1))
    create.mat <- function(w) {
      c(rep(w, k), rep(1 - w, K - k))
    }
    for(k in 1:(K - 1)) {
      p.hat <- predict(fits[[k]], newx=Z.test[, , k], 
                       s=min.AIC[k], type="response") #predict y<=k or y>k
      p.mat[, , k] <- sapply(p.hat, create.mat)
    }
    scores <- t(apply(p.mat, c(1, 2), sum)) #aggregate binary predictions
    class <- apply(scores, 1, which.max)
    return(list(scores=scores, class=class))
  }
  if(niter > 1) {
    scores <- lapply(predicted[, 1], as.matrix)
    agg.scores <- Reduce('+', scores)
   # if (type=="class") {
      class <- apply(agg.scores, 1, which.max)
   #   return(pred.class=class)
   # } else {
    return(list(scores=scores,
                agg.scores=agg.scores,
                pred.class=levels(object$y)[class]))
    #}
  } else {
    #if (type=="class") {
    #  return(pred.class=predicted$class)
    #} else {
    return(list(scores=predicted$scores,
                pred.class=levels(object$y)[predicted$class]))
    #}
  }
}


