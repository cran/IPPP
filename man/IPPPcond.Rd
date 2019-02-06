\name{IPPPcond}
\alias{IPPPcond}

\title{Simulate Location of Events for a Fixed Number of Events}
\description{
Returns a vector of random numbers that correspond to the location of events of an inhomogeneous Poisson point process (IPPP) in the case that the number of events is fixed.

The IPPP is described by a rate function r.
}
\usage{IPPPcond(samplesize, xrate, yrate, rnpr = 100)}

\arguments{
  \item{samplesize}{Natural number. The number of events occurring}
  \item{xrate}{Vector of (strictly increasing) real numbers}
  \item{yrate}{Vector of positive real numbers of the same length as xrate. The vectors xrate and yrate form the rate function r in the sense that r=approxfun(xrate,yrate)}
  \item{rnpr}{OPTIONAL, default is 100. The number of random numbers used per run in the loop of the rejection method. For details see the corresponding preprint}
}
\details{
Below min(xrate) and above max(xrate), the rate function r is assumed to have the value zero.
}
\value{A vector of length samplesize, containing the locations of the simulated events}

\author{Niklas Hohmann}
\references{
Hohmann, Niklas. "Conditional Densities and Simulations of Inhomogeneous Poisson Point Processes." arXiv 2019. <arXiv:1901.10754>}

\seealso{
\code{\link{IPPPuncond}} for the unconditioned version with a random number of events occurring
}

\examples{
#define rate function
sx=1:5
sy=c(0,1,1,3,0)
sm=c(1,0,1,0,-1)
xrate=seq(1,5,length.out=100)
yrate=splinefunH(sx,sy,sm)(xrate)
#plot rate function
plot(xrate,yrate,type='l',xlim=c(0.5,5.5), main='Rate Function')
#simulate location of 5 events:
samplesize=5
p5=IPPPcond(samplesize,xrate,yrate) 
points(p5,rep(0,samplesize),cex=2) #plot results
#simulate location of 30 events
samplesize=30
p30=IPPPcond(samplesize,xrate,yrate) 
points(p30,rep(1,samplesize),cex=2,pch=3) #compare with former results
}

\keyword{ Poisson }
\keyword{ IPPP }