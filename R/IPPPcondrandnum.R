IPPPcondrandnum <-
function(pointlocation,xrate,yrate,no=1,nthpoint=1,mode='forward'){
  if      (mode=='backward'){ l=-1  }
  else if (mode=='forward') { l=1   }
  else{ stop("Invalid mode. Either 'forward' or 'backward'.") }
  
  yrate = approx(xrate,yrate,xout=c(pointlocation,max(xrate)+100*diff(range(xrate)),min(xrate)-100*diff(range(xrate)),xrate),yleft=yrate[1],yright=yrate[length(xrate)])
  so = sort(yrate[[1]],index.return=TRUE)
  yrate = (yrate[[2]][so$ix])[!duplicated(so$x)]
  xrate = so$x[!duplicated(so$x)]
  
  Fval0=c(0,cumsum( (0.5*(c(0,yrate)+c(yrate,0))*(c(xrate,0)-c(0,xrate)))[2:length(xrate)] ))
  a = (0.5*(c(yrate,0)-c(0,yrate))/(c(xrate,0)-c(0,xrate)))[2:length(xrate)];#determine the (interval wise) leading coefficients of F_0 (which are 0.5 times the gradient of f_0)
  b = yrate[1:length(xrate)-1]-2*a*xrate[1:length(xrate)-1];                     #determine coefficients of linear parts of F_0
  c = -1*a*xrate[1:length(xrate)-1]^2-b*xrate[1:length(xrate)-1]+Fval0[1:length(xrate)-1]
  
  #determine sample relative to the transformed pointlocation
  sam = Fval0[which(xrate==pointlocation)]+l*rgamma(no,shape=nthpoint,rate=1)
  sieve = rep(0,length(sam))
  for (j in 1:length(xrate)){sieve = sieve + as.numeric( sam >= Fval0[j] )  }

  root1 = (-1*b[sieve]+sqrt(b[sieve]^2-4*a[sieve]*(c[sieve] - sam)))/(2*a[sieve])          #solve quadratic equations
  mysample = (-1*b[sieve]-sqrt(b[sieve]^2-4*a[sieve]*(c[sieve] - sam)))/(2*a[sieve])
  mysample[(a[sieve]==0)] = -1*((c[sieve] - sam)[(a[sieve]==0)])/((b[sieve])[a[sieve]==0])         #when an interval of f_0 is constant, F_0 is linear and the inversion can be done analytically
  mysample[(root1>=xrate[sieve])==TRUE & (root1<=xrate[sieve+1])==TRUE] = root1[(root1>=xrate[sieve])==TRUE & (root1<=xrate[sieve+1])==TRUE]
  return(mysample)
}
