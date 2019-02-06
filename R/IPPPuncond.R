IPPPuncond <-
function(xrate,yrate,expsamplesize=NULL,rnpr=100){
#Determine sample size of output
  if (is.null(expsamplesize)==TRUE){ #If no expected sample size is given
    expsamplesize = sum( (0.5*(c(0,yrate)+c(yrate,0))*(c(xrate,0)-c(0,xrate)))[2:length(xrate)] ); #get expected sample size as the integral over the rate function. integration is done using the trapezoidal rule
  }
  samplesize=rpois(1,lambda=expsamplesize) #determine the samplesize
  #Determine location of the samples for the output
  rvals=numeric()
  #rejection method
  while (length(rvals)<samplesize){ #if not enough random numbers have been generated
      lx=runif(min=min(xrate),max=max(xrate),rnpr)     #get nrp uniform random x values to minimize no. of runs through the loop
      rvals=c(rvals,lx[ runif(min=0,max=max(yrate),rnpr) <= approx(xrate,yrate,xout=lx)[[2]] ])    #get the lx values for which the corresponding ly values are below f(lx), reject the others 
  }
return(rvals[0:samplesize])
}
