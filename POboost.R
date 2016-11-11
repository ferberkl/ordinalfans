### Proportional Odds Boosting

PO.boost<-function(x, y, eps=0.1, w=1, mstop) {
  # Proportional odds (P/O) boosting algorithm
  #
  # Args:
  #    x:      Design matrix
  #    y:      Response vector
  #    eps:    A real-valued step length factor
  #    w:      Vector of weights
  #    K:      Number of levels in the response
  #    mstop:  Stopping iteration
  # Returns:
  #   coefs: Vector of coefficient estimates
  #   f.hat: Vector of function estimate
  #   theta: Vector of threshold estimates
  ###
  # Step 1: Initialize the n-dimensional vector f_hat[0] and the 
  #  K−1 threshold parameter estimates with offset values.
  #y <- factor(y, ordered=TRUE)
  K <- length(unique(y))
  p <- dim(x)[2] / (K - 1)
  x <- scale(x, center=T, scale=F)
  f.hat <- rep(0, length(y))
  pi.0  <-  table(y) / length(y)
  theta <- matrix(0, nrow=mstop, ncol=K - 1)
  theta[1, ] <- delta <- log(cumsum(pi.0) / (1 - cumsum(pi.0)))[1:(K - 1)]
  colnames(theta) <- paste("theta", 1:(K - 1), sep="")
  # Matrix of coefficient estimates (in all m steps)
  coefs <- matrix(0, nrow=mstop, ncol=dim(x)[2])  
  # Matrix of latest, updated vector of coefficient estimates
  coefs.latest <- rep(0, dim(x)[2])  
  colnames(coefs) <- names(coefs.latest) <- colnames(x)
  all.loss <- matrix(0, nrow=mstop, ncol=length(y)) 
  # Step 2: Specify base-learners, set m = 0.
  vars <- list()
  for (w in 1:p) {
    vars[[w]] <- colnames(x)[(w + (w-1)):(w + (w-1) + (K-2))]
  }
  m <- 0
  # Iterate mstop times.
  for (j in 1:mstop) {
    # Step 3: Increase m by 1.
    m <- m + 1
    # Step 4a: Calculate negative gradient vector with current estimate of
    # theta and f
    U <- neg.grad(y=y, f=f.hat, theta = theta[m, ])
    # Step 4b: Fit the negative gradient vector U[m] using each of the p 
    # base learners. This yields p vectors of predicted values, where each 
    # vector is an estimate of the negative gradient vector U[m]. 
    # Use R^2 to determine best base learner.
    fit.rsq <- function(aug.vars) {
      # Fit a linear model for a given base learner
      #
      # Args:
      #    aug.vars: Names of the augmented features used to fit the model
      # Returns:
      #   R^2 for given model.
      # fit.i <- glmnet(x=data.matrix(x[, aug.vars]), y=U, 
      #                 family="gaussian", alpha = 0)
      fit.i <- lm(U ~ 0 + ., data=data.frame(U, x[, aug.vars]))
      return(summary(fit.i)$r.squared)
      #return(max(fit.i$dev.ratio)) # Return max R^2 of all lambda values used
    }
    r.sq <- vapply(vars, fit.rsq, 1) # R^2 for all p base learners
    # Step 4c: Select the base-learner that fits U[m] best according to the Rsq 
    # goodness-of-fit criterion. Set U_hat[m] equal to the fitted values of the 
    # best model.
    max.rsq <- which.max(r.sq) # Base learner with largest R^2
    fit.maxrsq <- lm(U ~ 0 + ., data=data.frame(U, x[, vars[[max.rsq]]]))
    u.hat <- predict(fit.maxrsq, newdata = data.frame(x[, vars[[max.rsq]]]),
                     interval="none")
    vars.update <- names(fit.maxrsq$coefficients)
    coefs.latest[vars.update] <- coefs.latest[vars.update] + 
                                 eps * fit.maxrsq$coefficients
    coefs[m, ] <- coefs.latest
    # Step 4d:  Update f_hat[m] ← f_hat[m−1] + eps * U_hat[m], where 0 < eps ≤ 1
    #   is a real-valued step length factor.
    f.hat <- f.hat + eps * u.hat
    # Step 5: Plug f_hat[m] into the empirical risk function and minimize the
    #   empirical risk over θ. Set θ_hat[m] equal to the newly obtained estimate 
    #   of θ. (Convert delta (unconstrained) to theta, minimize risk over theta,
    #   return optimum delta).
    delta <- optim(par = delta, fn = riskS, y = y,
              fit = f.hat, w = w, method = "BFGS")$par
    # Convert optimum delta (above) to theta (constrained).
    theta[m, ] <- d2t(delta)  
    all.loss[m, ] <-  t(plloss(theta=theta[m, ], y=y, f=f.hat))
  }
  return(list=list(coefs=coefs, f.hat=f.hat, theta=theta, loss=all.loss))
}
