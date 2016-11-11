### Plot the model output

plot.ordinalFANS <- function(object, type="coefficients", xlab=NULL, 
                             ylab=NULL, main=NULL) {
  if (is.null(xlab)) xlab="Step"
  if (is.null(ylab)) {
    if (type == "coefficients") {
      ylab <- expression(hat(beta))
      if (is.null(main)) {
       main <- "Coefficient Path"
      }
    } else if (type == "risk") {
      ylab <- "Empirical Risk"
      if (is.null(main)) {
        main <- "Empirical Risk Path"
      }
    } else if (type=="var.importance") {
      ylab <- "Variable Importance"
      if (is.null(main)) {
        main <- "Variable Importance Path"
      }
    }
  }
  if (type=="coefficients") {
    if (object$niter > 1) {
      mean.beta <- object$mean.coefs
    } else {
      mean.beta <- object$coefs
    }
    last.row <- nrow(mean.beta)
    y.positions <- mean.beta[last.row, ][which(mean.beta[last.row, ] != 0)]
    varorder <- names(y.positions)
    matplot(mean.beta, type="l", xlab=xlab, ylab=ylab, main=main)
    axis(side=4, labels=varorder, at=y.positions, cex.axis=0.69)
  } else if (type=="risk") {
      if (object$niter > 1) {
        mean.risk <- object$mean.risk
      } else {
        mean.risk <- object$risk
      }
      plot(mean.risk, type="b", xlab=xlab, ylab=ylab, main=main)
  } else if (type=="var.importance") {
      if (object$niter > 1) {
        mean.importance <- object$mean.var.importance
      } else {
        mean.importance <- object$var.importance
      }
      last.row <- nrow(mean.importance)
      y.positions <- mean.importance[last.row, ][
                      which(mean.importance[last.row, ] != 0)]
      varorder <- names(y.positions)
      matplot(mean.importance, type="l", xlab=xlab, ylab=ylab, main=main)
      axis(side=4, labels=varorder, at=y.positions, cex.axis=0.69)
    }
}