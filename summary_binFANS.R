### Summary of the model fitting

summary.binFANS <- function(object) {
  features <- colnames(coef.binFANS(object, model.select="mean",
                       only.nonzero=TRUE))[-1]
  return(list(object=object, features=features))
}

