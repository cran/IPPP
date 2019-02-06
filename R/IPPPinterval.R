IPPPinterval <-
function(from, to, xrate,yrate,no=1,expsamplesize=NULL){
  #make sure the input "from" and "to" are within the range of xrate
  adjinterval=pmin(pmax(c(from,to),min(xrate)),max(xrate))
  #get all necessary function values, then sort them and delete duplicates
  yrate=approx(xrate,yrate,xout=c(adjinterval,xrate))
  so=sort(yrate[[1]],index.return=TRUE)
  yrate=(yrate[[2]][so$ix])[!duplicated(so$x)]
  xrate=so$x[!duplicated(so$x)]
  #integrate using the trapezoidal rule
  Fval0=c(0,cumsum( (0.5*(c(0,yrate)+c(yrate,0))*(c(xrate,0)-c(0,xrate)))[2:length(xrate)] ))
  #get rescaling factors to adjust for expsamplesize if necessary
  l=1
  if (is.null(expsamplesize)==FALSE){ l=expsamplesize/(Fval0[length(xrate)]) }
  #result are "no" Poisson distributed random numbers with parameter according to the integral from "from" to "to", (potentially rescaled by "l")
  return(rpois(no,lambda=l*(Fval0[which(xrate==adjinterval[2])]-Fval0[which(xrate==adjinterval[1])])))
}
