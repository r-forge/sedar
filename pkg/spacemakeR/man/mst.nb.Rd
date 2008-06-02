\name{mst.nb}
\alias{mst.nb}
\title{ Function to compute neighborhood based on the minimum spanning tree. }
\description{
  Returns an object of the class nb (see spdep package).
}
\usage{
mst.nb(dxy)
}

\arguments{
  \item{dxy}{ A distance matrix based on spatial coordinates of samples }
}

\value{
  An object of the class nb.
}

\author{ Stephane Dray }
\seealso{ \code{\link[spdep]{tri2nb}}, \code{\link[spdep]{graph2nb}}, \code{\link{give.thresh}} }
\examples{
xy <- matrix(rnorm(60),30,2)
dxy <- dist(xy)
th <- give.thresh(dxy)
nb1 <- mst.nb(dxy)
nb1
wh1 <- which(as.matrix(dxy)==th,arr.ind=TRUE)
plot(nb1,xy,pch=20,cex=2,lty=3)
lines(xy[wh1[1,],1],xy[wh1[1,],2],lwd=2)
title(main="Maximum distance of the minimum spanning tree in bold")

}
\keyword{ spatial }
