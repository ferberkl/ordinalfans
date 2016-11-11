
### Extract coefficient estimates

coef.binFANS <- function(object, model.select="mean", only.nonzero=FALSE)  {
  if(is.numeric(model.select) & model.select > object$niter) {
    stop(paste("model.select cannot be greater than ", object$niter, sep=""))
  }
  if (object$niter > 1) {
    if (model.select=="mean") {
      mean.coefs <- Reduce('+', object$coefs) / object$niter
      if (only.nonzero) {
        beta <- mean.coefs[, which(apply(mean.coefs, 2, 
                                         function(x) abs(sum(x)))!=0)]
      } 
    } else if (is.numeric(model.select)) {
      beta <- object$coefs[[model.select]]
      if (only.nonzero) {
        beta <- beta[, which(apply(beta, 2, function(x) abs(sum(x)))!=0)]
      }
    } else if (model.select=="all") {
      beta <- object$coefs
      if (only.nonzero) {
        nonzero <- function(x) {
          x[, which(apply(x, 2, function(y) abs(sum(y)))!=0)]
        }
        beta <- lapply(beta, function(l) nonzero(l))
      }
    }
  } else {
    beta <- object$coefs
    if (only.nonzero) {
      beta <- beta[, which(apply(beta, 2, function(x) abs(sum(x)))!=0)]
    } 
  }
  beta
}

