\name{give.thresh}
\alias{give.thresh}
\title{ Function that returns the maximum distance of the minimum spanning tree based on a distance matrix. }
\description{
  Can be useful to select a truncation value for the original PCNM approach. It returns the minimum value that keep all samples connected.
}
\usage{
give.thresh(distxy)
}

\arguments{
  \item{distxy}{ A distance matrix (i.e. obtained by dist() function) }
}
\value{
  The maximum value of the minimum spanning tree
}

\author{ Stephane Dray }
\examples{
xy <- matrix(rnorm(60),30,2)
dxy <- dist(xy)
th <- give.thresh(dxy)
nb1 <- mst.nb(dxy)
wh1 <- which(as.matrix(dxy)==th,arr.ind=TRUE)
plot(nb1,xy,pch=20,cex=2,lty=3)
lines(xy[wh1[1,],1],xy[wh1[1,],2],lwd=2)
title(main="Maximum distance of the minimum spanning tree in bold")


}
\keyword{ spatial }
