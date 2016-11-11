### Plot coefficient paths of individual binary response models

plot.binFANS <- function(object, type="coefficients", xlab=NULL, 
                         ylab=NULL, main=NULL, ...) {
  num.plots <- object$niter * (object$K - 1)
  if (object$niter > 1) {
    for (l in 1:object$niter) {
      for (k in 1:(object$K - 1)) {
        if(is.null(main))  {
          title <- paste(paste("FANS Iteration: ", l), 
                         paste("Modeling Augmented Feature: ", "", 
                         sep=as.character(k)), sep="\n")
        }
        x11()
        par(oma=c(0,0,1,1))
        plot(object$fits[[l]][[k]], main="", breaks=FALSE, type=type, ...)
        title(main=title)
      }
    }
  } else {
    for (k in 1:num.plots) {
      if (is.null(main)) {
        title <- paste(paste("FANS Iteration: ", object$niter), 
                       paste("Modeling Augmented Feature: ", "", 
                       sep=as.character(k)), sep="\n")
      }
      x11()
      par(oma=c(0,0,1,1))
      plot(object$fits[[k]], main="", breaks=FALSE, type=type, ...)
      title(main=title)
    }
  }
}
