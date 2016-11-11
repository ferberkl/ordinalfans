### Predict the outcome of a new observation

predict.ordinalFANS <- function(object, newx, mstop=NULL) {
  newx <- predict(object$trans, newx)  # Scale newx with mean and sd of X
  
  if (is.null(mstop)) mstop <- object$mstop
  if (!is.null(object$newx)) {
    if (object$niter > 1) {
      if (!is.null(object$mean.posteriors) & 
          all.equal(as.matrix(object$newx), as.matrix(newx))) {
        return(list=list(mean.posteriors=object$mean.posteriors,
               pred.class=object$pred.class))
      }
    } else {
      if (!is.null(object$posteriors) &
          all.equal(as.matrix(object$newx), as.matrix(newx))) {
        return(list=list(posteriors=object$posteriors,
               pred.class=object$pred.class))
      }
    }
  }
  K <- object$K
  niter <- object$niter

  niter.preds <- foreach(l = 1:niter, .combine='rbind', 
                     .packages=c('sm', 'caret')) %do% {
    if (object$niter > 1) {
      D.train.x <- object$D.train.x[[l]]
      D.train.y <- object$D.train.y[[l]]
      coefs <- object$coefs[[l]]
      theta <- object$theta[[l]]
      aug.X <- object$aug.X[[l]]
    } else {
      D.train.x <- object$D.train.x
      D.train.y <- object$D.train.y
      coefs <- object$coefs
      theta <- object$theta
      aug.X <- object$aug.X
    }
    test.marginals.f <- array(dim=c(dim(newx)[1], K - 1, dim(newx)[2]))
    test.marginals.g <- array(dim=c(dim(newx)[1], K - 1, dim(newx)[2]))
    for (k in 1:(K - 1)) {
      for (j in 1:dim(newx)[2]) {
        test.marginals.f[, k, j] <- sm.density(D.train.x[D.train.y <= k, j],
                                  eval.points=newx[, j],
                                  display="none")$estimate
        test.marginals.g[, k, j] <- sm.density(D.train.x[D.train.y > k, j],
                                  eval.points=newx[, j],
                                  display="none")$estimate
      }
    }
    # Winsorization to improve stability of estimates as suggested in FANS
    # manuscript
    test.marginals.g[test.marginals.g < 0.001] <- 
    test.marginals.f[test.marginals.f < 0.001] <- 0.001
    # Array of augmented features for newx. Within the array: p matrices of 
    # augmented features of dimension  n x (K - 1)
    Z.test <- log(test.marginals.f / test.marginals.g)
    full.test <- data.frame(Z.test)
    colnames(full.test) <- paste("z", paste(rep(1:dim(Z.test)[3], each=K - 1), 
                                            rep(1:(K - 1),
                                            times=dim(test.marginals.f)[3]), 
                                            sep="."),
                                 sep="")
    # function estimates for newx
    offset <- as.numeric(-1 * (matrix(apply(aug.X, 2, mean), nrow=1) %*%
                         matrix(coefs[mstop, ], ncol=1)))
    f.newx <- offset + (as.matrix(full.test) %*% as.matrix(coefs[mstop, ]))
    # estimated posterior probabilities for newx
    post.probs <- response(f.newx, theta[mstop, ])  
    return(list=list(post.probs))
  }
  if (object$niter > 1) {
    if (ncol(as.matrix(niter.preds[[1]])) == 1) {
      mean.posteriors <- Reduce('+', niter.preds) / niter
      pred.class <- factor(levels(object$y)[apply(t(as.matrix(mean.posteriors)), 1, which.max)],
                           levels=levels(object$y), ordered=TRUE)
      posteriors <- lapply(1:niter, function(w) t(as.matrix(niter.preds))[[w]])
    } else {
      mean.posteriors <- Reduce('+', niter.preds) / niter
      pred.class <- factor(levels(object$y)[apply(mean.posteriors, 1, which.max)],
                           levels=levels(object$y), ordered=TRUE)
      posteriors <- lapply(1:niter, function(w) niter.preds[[w]])
    }
    return(list=list(posteriors=posteriors, 
                     mean.posteriors=mean.posteriors, pred.class=pred.class))
  } else {
    if (ncol(as.matrix(niter.preds[[1]])) == 1) {
      posteriors <- t(as.matrix(niter.preds[[1]]))
    } else {
      posteriors <- as.matrix(niter.preds[[1]])
    }
    pred.class <- factor(levels(object$y)[apply(posteriors, 1, which.max)],
                         levels=levels(object$y), ordered=TRUE)
    return(list=list(posteriors=posteriors, pred.class=pred.class))
  }
}
