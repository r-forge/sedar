\name{scores.listw}
\alias{scores.listw}

\title{ Function to compute Moran's eigenvectors of a listw object. }
\description{
  This functions compute eigenvector's of a doubly centered spatial weighting matrix. Corresponding eigenvalues are linearly related to Moran's index of spatial autocorrelation.
}
\usage{
scores.listw(listw, echo = FALSE, MEM.autocor = c("all","positive", "negative"))
}

\arguments{
  \item{listw}{ An object of the class \code{listw} created by functions of the \code{spdep} package. }
  \item{echo}{ A logical value, if TRUE a message indicates which eigenvectors corresponding to null eigenvalues are removed. }
\item{MEM.autocor}{ A string indicating if all MEMs must be returned or
  only those corresponding to positive or negative autocorrelation. }
}
\details{
  In order to remove the effect of units of measurments (cm, km ...), testing the nullity of eigenvalues is based on E(i)/E(1) where E(i) is i-th eiegnvalue and E(1) is the maximum absolute value of eigenvalues.
}
\value{
  A list:
  \item{vectors }{A matrix of eigenvectors,}
  \item{values }{A vector of eigenvalues.}
  \item{call }{Call of the function.}
}
\references{ Griffith D. A. (1996) Spatial autocorrelation and eigenfunctions of the geographic weights matrix accompanying geo-referenced data. \emph{Canadian Geographer} \bold{40}, 351--367. }
\author{ Stephane Dray }



\seealso{ \code{\link[spdep]{nb2listw}} }
\examples{
data(oribatid)
fau <- sqrt(oribatid$fau/outer(apply(oribatid$fau,1,sum),rep(1,ncol(oribatid$fau)),"*"))
faudt <- resid(lm(as.matrix(fau)~as.matrix(oribatid$xy)))
nbtri <- tri2nb(as.matrix(oribatid$xy))
sc.tri <- scores.listw(nb2listw(nbtri,style="B"))
par(mfrow=c(1,3))
s.value(oribatid$xy,sc.tri$vectors[,1],cleg=0, sub="Vector 1",csub=3)
s.value(oribatid$xy,sc.tri$vectors[,2],cleg=0, sub="Vector 2",csub=3)
s.value(oribatid$xy,sc.tri$vectors[,3],cleg=0, sub="Vector 3",csub=3)


    
  }

\keyword{ spatial }
