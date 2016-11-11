### Does the same thing as predict()

fitted.ordinalFANS <- function(object, newx, mstop=NULL) {
  predict.ordinalFANS(object=object, newx=newx, mstop=mstop)
}
