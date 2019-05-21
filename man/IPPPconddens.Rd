\name{IPPPconddens}
\alias{IPPPconddens}
\title{Conditional Probability Density Given the Location of one Event}
\description{Given knowledge of the location of one event, this function determines the probability density function (pdf) of the distribution of the n-th point above/below
the known event.}
\usage{IPPPconddens(x, pointlocation, xrate, yrate, nthpoint = 1, mode = "forward")}

\arguments{
  \item{x}{Points at which the values of the pdf are determined}
  \item{pointlocation}{Location of the known event}
  \item{xrate}{Vector of (strictly increasing) real numbers}
  \item{yrate}{Vector of strictly positive real numbers of the same length as xrate. The vectors xrate and yrate form the rate function r in the sense that r=approxfun(xrate,yrate)}
  \item{nthpoint}{OPTIONAL, default is 1. Setting this to a value i will determine the pdf of the ith event above/below the event at the position pointlocation.}
  \item{mode}{OPTIONAL, default is "forward". Determines whether the pdf for points above(right) the known event is determined (mode="forward"), or whether the pdf for points below(left) the known event is determined (mode="backward")}
}
\details{The value of the rate function r below min(xrate) is set to r(min(xrate)), and for values above max(xrate) it is set to r(max(xrate)). Both r(min(xrate)) and r(max(xrate)) need to be strictly positive for the results to be correct.}
\value{Returns a list containing two entries:
\item{ x }{Duplicate of the input of the same name}
\item{ densval }{Vector consisting of the values of the conditional probability density, evaluated at x}
}

\author{Niklas Hohmann}

\references{
Hohmann, Niklas. "Conditional Densities and Simulations of Inhomogeneous Poisson Point Processes: The R package "IPPP"" arXiv 2019. <arXiv:1901.10754>
}

\seealso{
\code{\link{IPPPnthpointdens}} for the probability density of the location of the of the n-th event, given that a fixed number of events occur.

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
#simulate the pdf of the event above the known event:
r1=IPPPconddens(x,pointlocation,xrate,yrate)
lines(r1$x,r1$densval,col="red",lwd=3)
#simulate the pdf of the event below the known event:
r2=IPPPconddens(x,pointlocation,xrate,yrate,mode="backward")
lines(r2$x,r2$densval,col="blue",lwd=3)
#simulate the pdf of the third event below the known event:
r2=IPPPconddens(x,pointlocation,xrate,yrate,mode="backward",nthpoint=3)
lines(r2$x,r2$densval,col="green",lwd=3)
}

\keyword{ Poisson }
\keyword{ Poisson point process}
\keyword{ Point process }
\keyword{ IPPP }
\keyword{ density }
