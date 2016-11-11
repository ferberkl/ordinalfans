### Summarize the fitted model

summary.ordinalFANS <- function(object, m=NULL) {
  if (object$L > 1) {
    if (is.null(m)) {
      m <- dim(object$coefs[[1]])[1]
    }
    rownames(object$mean.var.importance) <- NULL
    var.importance <- object$mean.var.importance[m,
                                    which(object$mean.var.importance[m, ]!=0)]
  } else {
      if (is.null(m)) {
        m <- dim(object$coefs)[1]
      }
      rownames(object$var.importance) <- NULL
      var.importance <- object$var.importance[m,
                                         which(object$var.importance[m, ]!=0)]
  } 
  return(list(object=object, var.importance=var.importance))
}