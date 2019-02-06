IPPPcond <-
function(samplesize,xrate,yrate,rnpr=100){
  rvals=numeric()
  #rejection method
  while (length(rvals)<samplesize){ #if not enough random numbers have been generated
    lx=runif(min=min(xrate),max=max(xrate),rnpr)     #get nrp uniform random x values to minimize no. of runs through the loop
    rvals=c(rvals,lx[ runif(min=0,max=max(yrate),rnpr) <= approx(xrate,yrate,xout=lx)[[2]] ])    #get the lx values for which the corresponding ly values are below f(lx), reject the others 
  }
  #return the desired number of random numbers
  return(rvals[0:samplesize])
}
