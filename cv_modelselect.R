### Select stopping iteration by cross validation

cv.modelselect <- function(x, y, parallel=TRUE, num.folds=5, seed=2468, 
                           mstop.seq=floor(seq(from=10, to=100, by=10))) {
  set.seed(seed)
  folds <- createFolds(y=y, k=num.folds, list=F)
  if (parallel) {
    `%fun%` <- `%dopar%`
  } else {
    `%fun%` <- `%do%`    
  }
    modelselect <- foreach (m=mstop.seq, .combine=rbind) %fun% {
      dxy <- rep(0, num.folds)
      for (i in 1:num.folds) {
        fit <- ordinalFANS(x=x[-which(folds == i), ], y=y[-which(folds == i)],
                           newx=x[which(folds == i), ], niter=1, mstop=m)
        dxy[i] <- rcorr.cens(as.numeric(fit$pred.class),
                             y[which(folds == i)])[2]
      }
      CV.estimate <- sum((table(folds) / length(y)) * dxy)
      dxy.sd <- sqrt(sum((table(folds) / length(y)) * (dxy - mean(dxy))^2))
      return(c(m, CV.estimate, dxy.sd))
    }
  colnames(modelselect) <- c("mstop", "CV Somer's Dxy Estimate", "SE")
  return(modelselect)
}
