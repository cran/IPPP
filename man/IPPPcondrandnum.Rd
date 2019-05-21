\name{IPPPcondrandnum}
\alias{IPPPcondrandnum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Events Relative to a Known Event of an Inhomogeneous Poisson Point Process
}
\description{
Given knowledge of the location of one event from an inhomogeneous Poisson point process (IPPP), this function determines random numbers corresponding to the location of the n-th event above/below the known event.

The IPPP is described by a rate function r.
}
\usage{IPPPcondrandnum(pointlocation, xrate, yrate, no = 1, nthpoint = 1, mode = "forward")}

\arguments{
  \item{pointlocation}{Location of the known event}
  \item{xrate}{Vector of (strictly increasing) real numbers}
  \item{yrate}{Vector of strictly positive real numbers of the same length as
xrate. The vectors xrate and yrate form the rate function r in the sense
that r=approxfun(xrate,yrate)}
  \item{no}{OPTIONAL, default is 1. The number of random events to generate}
  \item{nthpoint}{OPTIONAL, default is 1. Setting this to a value i will simulate the location of the ith event above/below the event at the position pointlocation}
  \item{mode}{OPTIONAL, default is "forward". Determines whether the pdf for points above(right) the known event is determined (mode="forward""), or whether the pdf for points below(left) the known event is determined (mode="backward"")}
}
\details{The value of rate function r below min(xrate) is set to r(min(xrate)), and for values above max(xrate) it is set to r(max(xrate)). Both r(min(xrate)) and r(max(xrate)) need to be strictly positive for the results to be correct.}
\value{
Vector of length no, containing the locations of the simulated events.
}

\author{Niklas Hohmann}

\references{
Hohmann, Niklas. "Conditional Densities and Simulations of Inhomogeneous Poisson Point Processes: The R package "IPPP"" arXiv 2019. <arXiv:1901.10754>
}

\seealso{
\code{vignette("IPPP")} for an overview of the features of the IPPP package and some background.
}

\examples{
#define rate function
xrate=seq(0,2*pi,length.out=1000)
yrate=sin(xrate)+1.01
plot(xrate,yrate,type="l",main="Rate Function")
#define known event:
pointlocation=3
lines(c(3,3),c(0,3),lwd=3)
#values where the density is calulated:
x=xrate
#simulate the next point above the known event:
p1=IPPPcondrandnum(pointlocation, xrate,yrate)
points(p1,0.5,pch=1,cex=2)
#simulate 10 times the next event below the known location:
p2=IPPPcondrandnum(pointlocation, xrate,yrate,no=10,mode="backward")
points(p2,rep(1,length(p2)),cex=2,pch=2)
#simulate the second point above the known event:
p3=IPPPcondrandnum(pointlocation, xrate,yrate,no=1,nthpoint=2)
points(p3,1.5,cex=2,pch=3) #might be out of the plot
}

\keyword{ Poisson }
\keyword{ Poisson point process}
\keyword{ Point process }
\keyword{ IPPP }
