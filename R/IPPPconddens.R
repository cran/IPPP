IPPPconddens <-
function(x,pointlocation,xrate,yrate,nthpoint=1,mode="forward"){
  #determine whether density for events below ("backward") or above ("forward") the given point "pointlocation" is determined
  if      (mode=="backward"){ l=-1  }
  else if (mode=="forward") { l=1   }
  else{ stop("Invalid mode. Either 'forward' or 'backward'.") }
  
  #prepare values for integration: get all function values, then sort them and delete duplicates
  yrate=approx(xrate,yrate,xout=c(x,pointlocation,xrate),yleft=yrate[which(xrate==min(xrate))],yright=yrate[which(xrate==max(xrate))])
  so=sort(yrate[[1]],index.return=TRUE)
  yrate=(yrate[[2]][so$ix])[!duplicated(so$x)]
  xrate=so$x[!duplicated(so$x)]
  
  #integrate
  Fval0=c(0,cumsum( (0.5*(c(0,yrate)+c(yrate,0))*(c(xrate,0)-c(0,xrate)))[2:length(xrate)] ))
  
  #determine density at the points x
  res=approx(xrate,yrate,xout=x,yleft=yrate[1],yright=yrate[length(xrate)])[[2]]*dgamma(l*(approx(xrate,Fval0,xout=x)[[2]]-approx(xrate,Fval0,xout=pointlocation)[[2]]),shape=nthpoint,rate=1)
  return(list(x=x,densval=res))
}
