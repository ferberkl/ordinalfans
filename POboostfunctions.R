
### Functions used in P/O Boosting

d2t <- function(delta) {
  # Used for constraining the threshold estimates to be nondecreasing
  #
  # Args:
  #   delta: the vector of unconstrained threshold values
  #
  # Returns:
  #   Theta, the vector of constrained threshold values
  delta[1] + cumsum(c(0, exp(delta[-1])))
}

plloss <- function(theta, y, f, w = w) {
  # Calculates the value of the loss (negative logL)
  #
  # Args:
  #    theta: the vector of constrained threshold values
  #    y:      the response vector
  #    f:      the current function estimate
  #    w:      the vector of weights
  #
  # Returns:
  #   The value of the loss function at the current step
  if (length(f) == 1) f <- rep(f, length(y))
  tmp <- lapply(1:(length(theta) + 1), function(i) {
    if (i == 1) return(1 + exp(f - theta[i]))
    if (i == (length(theta) + 1)) {
        return(1 - 1 / (1 + exp(f - theta[i - 1])))
    }
    return(1 / (1 + exp(f - theta[i])) -
             1 / (1 + exp(f - theta[i - 1])))
  })
  loss <- log(tmp[[1]]) * (y == levels(y)[1])
  for (i in 2:nlevels(y)) {
    loss <- loss - log(tmp[[i]]) * (y == levels(y)[i])
  }
  return(loss)
}

riskS <- function(delta, y, fit, w = w) {
  # Calculates the value of the empirical risk.
  #   Converts delta to theta (constrained), then sums Loss over i=1, ..., n.
  #
  # Args:
  #    delta: the vector of unconstrained threshold values
  #    y:      the response vector
  #    fit:    the current function estimate
  #    w:      the vector of weights
  #
  # Returns:
  #   The value of the risk function at the current step
  sum(w * plloss(y = y, f = fit, theta = d2t(delta)))
}

neg.grad <- function(y, f, theta = theta, w = w) {
  # Calculates the negative gradient at current estimate of f and theta
  #
  # Args:
  #    y:      the response vector
  #    f:      the current function estimate
  #    theta: the vector of constrained threshold values
  #    w:      the vector of weights
  #
  # Returns:
  #   The value of the risk function at the current step
  if (length(f) == 1) f <- rep(f, length(y))
  # Calculates negative gradient for each subject, then sums over i=1...n
  ng <- sapply(1:(length(theta) + 1), function(i) {
    if (i > 1 & i < (length(theta) + 1)) {
      ret <- (1 - exp(2 * f - theta[i - 1] - theta[i]))  /
             (1 + exp(f - theta[i - 1])  +
                  exp(f - theta[i])  +
                  exp(2 * f - theta[i - 1] - theta[i]))
    } else {
        if (i == 1) {
          ret <- -1 / (1 + exp(theta[i] - f))
        } else {
          ret <- 1 / (1 + exp(f - theta[i - 1]))
        }
    }
    return(ret * (y == levels(y)[i]))
  })
  rowSums(ng)
 }

#Empirical risk function, treat theta as fixed, optimize over f - used to initialize f
# risk <- function(y, f, w = w)
#     sum(w * plloss(y = y, f = f, theta = theta))

#Initialize delta / theta and f (maybe just initialize f to 0)
#initial.f <- function(y, w = w) {
#  delta<-log(cumsum(pi.0) / (1 - cumsum(pi.0)))[1:(K - 1)]
#  optimize(risk, interval = c(-5, 5), y = y, w = w)$minimum
# }

response <- function(f, theta) {
   # Calculates posterior probabilities
  #
  # Args:
  #    f: the current function estimate
  #
  # Returns:
  #   Vector of posterior probabilities (of length K)
  ret <- sapply(1:(length(theta) + 1), function(i) {
    if (i == 1) return(1 / (1 + exp(f - theta[i])))  #P(Y=1|X)
    if (i == (length(theta) + 1)) {
      return(1 - 1 / (1 + exp(f - theta[i - 1])))    #P(Y=K|X)
    }
    return(1 / (1 + exp(f - theta[i])) -      #P(Y=2or3or...orK-1|X)
       1 / (1 + exp(f - theta[i - 1])))
    })
    ret
   }
