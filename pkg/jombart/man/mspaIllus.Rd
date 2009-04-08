\encoding{UTF-8}
\name{mspaIllus}
\alias{mspaIllus}
\docType{data}
\title{Simulated data illustrating MSPA}
\description{
  These data simulate multi-scale spatial patterns, obtained by linear
  combination of Moran's eigenvectors. See the \code{example} to produce
  analogous dataset. See \code{\link{mspa}} manpage for an analysis of
  this dataset.
}
\usage{
data(mspaIllus)
}
\format{
  \code{mspaIllus} is data.frame with 100 observations (rows) and 35
  variables (columns). Observations are considered to be spatially
  distributed on a regular grid. As an extra attribute, it contains a list of
  spatial weights derived from this grid (\code{attr(mspaIllus,"lw")}).
}
\author{
  Thibaut Jombart \email{t.jombart@imperial.ac.uk}
}
\seealso{ \code{\link{mspa}} }
\references{
  Jombart T, Dray S, and Dufour, A-B. (2009) Finding essential scales of
  spatial variation in ecological data: a multivariate
  approach. \emph{Ecography} \bold{32}: 161-168.
}
\examples{
### This is how data were obtained ###
##
## output object is 'X'

## spatial coordinates
xy <- expand.grid(1:10,1:10)

## list of spatial weights
lw <- nb2listw(cell2nb(10,10))

## Moran's eigenvectors
U <- orthobasis.listw(lw)

## data.frame of simulated variables
V <- list()
# large scale
V$V1 <- 0.5*U[,1] + 0.5*U[,2] + 0.5*U[,3] + rnorm(100)
V$V2 <- 0.5*U[,1] - 0.8*U[,2] + 0.5*U[,3] + rnorm(100)
V$V3 <- U[,1] - U[,2] + 0.5*U[,3] + rnorm(100)

# medium scale
V$V4 <- 0.6*U[,44] + U[,45] +  0.8*U[,46] + rnorm(100)

# fine scale
V$V5 <-  0.5*U[,97] + U[,98] + U[,99]  + rnorm(100)
V$V6 <- 0.5*U[,97] + 0.5*U[,98] - U[,99] + rnorm(100)
V$V7 <- 0.6*U[,97] + 0.6*U[,98] + 0.8*U[,99]  + rnorm(100)

V <- as.data.frame(V)
M <- matrix(rnorm(100*28),nrow=100)
X <- cbind.data.frame(V,M)
colnames(X) <- paste("V",1:ncol(X),sep='')
attr(X,"lw") <- lw
}
\keyword{datasets}
\keyword{spatial}

