\encoding{UTF-8}
\name{xytrend}
\alias{xytrend}
\title{Trend Surface for Spatial Coordinates}
\description{
Makes trend surface using third-degree polynomial of spatial terms.
}
\usage{
xytrend(xy)
}
\arguments{
  \item{xy}{a matrix with 2 columns for x and y coordinates.}
}
\value{
A matrix with 9 columns.
}
\details{
The x and y coordinates are standardized (scaled to zero mean and unit variance) by the function \code{\link[vegan]{decostand}}. Polynomial terms are: $x$, $y$, $x^2$, $xy$, $y^2$, $x^3$, $x^2*y$, $x*y^2$, $y^3$, according to Legendre (1993).
}
\references{
Legendre, P. 1993. Spatial autocorrelation: trouble or new paradigm? \emph{Ecology} 74, 1659--1673.
}
\author{
\enc{Péter Sólymos}{Peter Solymos}, \email{solymos@ualberta.ca}
}
\seealso{
  \code{\link[vegan]{decostand}}
}
\examples{
xy <- cbind(runif(20,0,100), runif(20,0,100))
xy
xytrend(xy)
}
\keyword{ manip }
\keyword{ spatial }

