### Extract Model Coefficients

coef.ordinalFANS <- function(object, model.select="mean", m=NULL, 
                             only.nonzero=FALSE) {
  if (object$niter > 1) {
    if (is.null(m)) {
      m <- dim(object$coefs[[1]])[1]
    }
    if (model.select=="mean") {
        beta <- object$mean.coefs[m, ]
        theta <- object$mean.theta[m, ]
        if (only.nonzero) {
          beta <- object$mean.coefs[m, ][object$mean.coefs[m, ]!=0]
        }
    } else if (is.numeric(model.select)) {
        beta <- object$coefs[[model.select]][m, ]
        theta <- object$theta[[model.select]][m, ]
        if (only.nonzero) {
          beta <- beta[beta!=0]
          theta <- theta[theta!=0]
        }
    } else if (model.select=="all") {
        beta <- lapply(object$coefs, function(x) {x[m, ]})
        theta <- lapply(object$theta, function(x) {x[m, ]})
        if (only.nonzero) {
          beta <- lapply(beta, function(x) {x[x!=0]})
          theta <- lapply(theta, function(x) {x[x!=0]})
        }
    }
  } else {
      if (is.null(m)) {
        m <- dim(object$coefs)[1]
      }
      beta <- object$coefs[m, ]
      theta <- object$theta[m, ]
      if (only.nonzero) {
        beta <- beta[beta!=0]
      }
    }
  c(theta, beta)
}
