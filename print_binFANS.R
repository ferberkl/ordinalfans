### Print a summary of the model fitting

print.binFANS <- function(object, ...) {
  cat("Call:\n")
  print(object$call)
  cat("\nNumber of FANS iterations: niter = ", object$niter)
  cat("\n")
}

