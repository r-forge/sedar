\name{xytrend}
\alias{xytrend}
\title{Trend Surface for Spatial Coordinates}
\description{
Makes trend surface using third-degree polynomial of spatial terms.
}
\usage{
xytrend(xy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xy}{a matrix with 2 columns for x and y coordinates.}
}
\value{
A matrix with.
}
\details{
The x and y coordinates are standardized (scaled to zero mean and unit variance) by the function \code{\link[vegan]{decostand}}. Polynomial terms are: $x$, $y$, $x^2$, $xy$, $y^2$, $x^3$, $x^2*y$, $x*y^2$, $y^3$, according to Legendre (1993).
}
\references{
Legendre, P. 1993. Spatial autocorrelation: trouble or new paradigm? Ecology 74, 1659--1673.
}
\author{Peter Solymos, \email{Solymos.Peter@aotk.szie.hu}, \url{http://www.univet.hu/users/psolymos/personal/}}
\seealso{
  \code{\link[vegan]{decostand}}
}
\examples{
(xy <- cbind(runif(20,0,100), runif(20,0,100)))
(xy2 <- xytrend(xy))
}
\keyword{ manip }
\keyword{ spatial }

