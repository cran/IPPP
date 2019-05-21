\name{IPPPuncond}
\alias{IPPPuncond}

\title{Simulate Events according to an Inhomogeneous Poisson Point Process}
\description{
Generates random numbers corresponding to the locations of events of an inhomogeneous Poisson point process (IPPP). The IPPP is described by a rate function r.
}
\usage{IPPPuncond(xrate, yrate, expsamplesize = NULL, rnpr = 100)}

\arguments{
  \item{xrate}{Vector of (strictly increasing) real numbers}
  \item{yrate}{Vector of positive real numbers of the same length as xrate. The vectors xrate and yrate form the rate function r in the sense that r=approxfun(xrate,yrate)}
  \item{expsamplesize}{OPTIONAL, default is NULL. If this variable is set to a numeric value, it determines the average number of events occurring}
  \item{rnpr}{OPTIONAL, default is 100. The number of random numbers used per run in the loop of the rejection method. For more details see the corresponding preprint}
}
\details{Below min(xrate) and above max(xrate), the rate function r is assumed to have
the value zero.}
\value{A vector of variable length, containing the generated locations of the events. If no events occur, the output is numeric(0)}

\author{Niklas Hohmann}

\references{
Hohmann, Niklas. "Conditional Densities and Simulations of Inhomogeneous Poisson Point Processes: The R package "IPPP"" arXiv 2019. <arXiv:1901.10754>
}

\seealso{
\code{\link{IPPPuncond}} for a version with a fixed number of events occurring.

\code{vignette("IPPP")} for an overview of the features of the IPPP package and some background.
}

\examples{
#Define rate function
sx=1:5
sy=c(0,1,1,3,0)
sm=c(1,0,1,0,-1)
xrate=seq(1,5,length.out=100)
yrate=splinefunH(sx,sy,sm)(xrate)
#plot rate function
plot(xrate,yrate,type='l',xlim=c(0.5,5.5), main='Rate Function') 

p=IPPPuncond(xrate,yrate) #simulate one set of events
points(p,rep(0,length(p)),cex=2) #plot results
#simulate location of events with the expected number of events being 20
expsamplesize=20
pp=IPPPuncond(xrate,yrate,expsamplesize) 
length(pp) #in most cases, the result is not exactely 20
points(pp,rep(1,length(pp)),cex=2,pch=3) #compare with former results
}

\keyword{ Poisson }
\keyword{ Poisson point process}
\keyword{ Point process }
\keyword{ IPPP }