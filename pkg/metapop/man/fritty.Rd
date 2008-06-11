\name{fritty}
\alias{fritty}
\docType{data}
\title{ Forged Data for Metapopulation Analysis}
\description{
 Data are fictitious. They were generated using
 \code{\link{metacycle}}. However, the patch network is made to resemble
 the example in Hanski et al. (1994).
}
\usage{data(fritty)}
\format{
  A data frame with 50 observations on the following 5 variables.
  \describe{
    \item{\code{x.crd, y.crd}}{ Coordinates of the patches (locations)}
    \item{\code{A}}{Areas of patches}
    \item{\code{p}}{Occupancy (0 or 1) in the patches.}
    \item{\code{p2}}{Occupancy in the same patches next year.}
  }
}

\references{
Hanski, I., Kuussaari, M. & Nieminen, M. 1994. Metapopulation structure
and migration in butterfly \emph{Melitaea cinxia}. \emph{Ecology} 75: 747--762.
}
\examples{
data(fritty)
str(fritty)
## Plot
with(fritty,
plot(x.crd, y.crd, asp = 1, xlab = "", ylab = "",
cex = sqrt(A *+ 5), pch = 21, col = p + 1, bg = 5 * p)
)
}
\keyword{datasets}
