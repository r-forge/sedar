\name{mspa}
\alias{mspa}
\alias{local.rtest}
\title{Global and local tests}
\description{
 
}
\usage{
mspa()
\method{scatter}{mspa}()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{}{}
  \item{listw}{a list of weights of class \code{listw}. Can be obtained
    easily using the function \code{chooseCN}.}
  \item{k}{integer: the number of highest \eqn{R^2}{R^2} summed to form the test statistics}
  \item{nperm}{integer: the number of randomisations to be performed.}
}
\details{
 These tests were originally part of the adegenet package for R.
}
\value{
An object of class \code{randtest}.
}
\references{
  Jombart T, Dray S, and Dufour, A-B. (2009) Finding essential scales of
  spatial variation in ecological data: a multivariate
  approach. \emph{Ecography} \blod{32}: 161-168.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\seealso{\code{\link{chooseCN}} to obtain a list of spatial weights.}
\examples{

}
\keyword{multivariate}
\keyword{spatial}
