IPPPnthpointdens <-
function(x,n,samplesize, xrate,yrate){
  yrate=approx(xrate,yrate,xout=c(pmin(pmax(x,min(xrate)),max(xrate)),xrate),yleft=0,yright=0)
  so=sort(yrate[[1]],index.return=TRUE)
  yrate=(yrate[[2]][so$ix])[!duplicated(so$x)]
  xrate=so$x[!duplicated(so$x)]
  
  Fval0=c(0,cumsum( (0.5*(c(0,yrate)+c(yrate,0))*(c(xrate,0)-c(0,xrate)))[2:length(xrate)] ))
  res=n*choose(samplesize,n)*((1/max(Fval0))*approx(xrate,Fval0,xout=x,yleft=0,yright=0)[[2]])^(n-1)*(1-(1/max(Fval0))*approx(xrate,Fval0,xout=x,yleft=0,yright=0)[[2]])^(samplesize-n)*(1/max(Fval0))*approx(xrate,yrate,xout=x,yleft=0,yright=0)[[2]]
  return(list(x=x,densval=res))
}
