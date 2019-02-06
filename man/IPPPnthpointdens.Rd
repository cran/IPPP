\name{IPPPnthpointdens}
\alias{IPPPnthpointdens}

\title{ Probability Density of Events in an Inhomogeneous Poisson Point Process}
\description{
Determines values of the probability density function of the n-th point of an inhomogeneous Poisson point process (IPPP), given that a fixed number of events occur.
The IPPP is described by a rate function r.
}
\usage{
IPPPnthpointdens(x, n, samplesize, xrate, yrate)
}

\arguments{
  \item{x}{Vector of real numbers where the value of the probability density is determined}
  \item{n}{Natural number smaller than samplesize. Determines the event for which the probability density is determined. 1=first event from the left, i=i-th event from the left.}
  \item{samplesize}{Natural number. The overall number of events occurring}
  \item{xrate}{Vector of (strictly increasing) real numbers}
  \item{yrate}{Vector of positive real numbers of the same length as xrate. The vectors xrate and yrate form the rate function r in the sense that r=approxfun(xrate,yrate)}
}
\details{ Below min(xrate) and above max(xrate), the rate function r is assumed to have the value zero.}
\value{
A list containing two entries:
\item{ x }{A duplicate of the input of the same name}
\item{ densval }{A vector consisting of the values of the probability density, evaluated at x}
}

\author{Niklas Hohmann}

\references{
Hohmann, Niklas. "Conditional Densities and Simulations of Inhomogeneous Poisson Point Processes." arXiv 2019. <arXiv:1901.10754>}


\seealso{ \code{\link{IPPPconddens}} for the probability density of the n-th event above/below the location of some given event. }
\examples{
sx=c(1,2,3,4,5)
sy=c(0,1,1,3,0)
sm=c(1,0,1,0,-1)
xrate=seq(1,5,length.out=100)
yrate=splinefunH(sx,sy,sm)(xrate)
#plot rate function
plot(xrate,yrate,type='l',xlim=c(0.5,5.5), main='Rate Function') 

x=seq(0,6,length.out = 1000) #point where the pdf is determined
n=1          #get the pdf of the first point from the left ...
samplesize=5  #... out of a sample of five
ll=IPPPnthpointdens(x,n,samplesize,xrate,yrate)
plot(ll$x,ll$densval,type='l') #plot the resulting pdf
legend('topleft',legend=paste('pdf of point no.',as.character(n),
'out of ',as.character(samplesize)), 'points' ,lty=1)
}

\keyword{ Poisson }
\keyword{ IPPP }
