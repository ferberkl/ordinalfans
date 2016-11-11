### Print the model output

print.ordinalFANS <- function(object, ...) {
  cat("Call:\n")
  print(object$call)
  cat("\nNumber of boosting iterations: mstop = ", object$mstop)
  cat("\nNumber of FANS iterations: niter = ", object$niter)
  cat("\nStep size: ", object$eps, "\n")
}