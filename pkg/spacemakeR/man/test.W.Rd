\name{test.W}
\alias{test.W}
\title{ Function to compute and test eigenvectors of spatial weighting matrices. }
\description{
  This function is a user-friendly way to compute and test eigenvectors for various definitions of spatial weighting matrices.
  It combines calls to the functions \code{scores.listw} and \code{ortho.AIC}.
  It allows to test various definitions of the spatial weighting matrix and return results of \code{scores.listw} for the best one.
}
\usage{
test.W(Y, nb, xy, f = NULL, ...)
}

\arguments{
  \item{Y}{ A matrix with response variables (univariate or multivariate response). }
  \item{nb}{ An object of the class \code{nb} created by functions of the \code{spdep} package. }
  \item{xy}{ Coordinates of the samples, this argument is optional and is required only if the argument \code{f} is not null.}
  \item{f}{ A function of the distance that can be used as a weighting spatial function. This argument is optional.  }
  \item{\dots}{ Others arguments for the function \code{f}. It defines the range of parameters which will be tested. }
}
\details{
  This functions allows to test one binary spatial weighting matrix (if only Y and nb are provided).
  It allows also to test a weighting function based on distances (if f is provided) and a weighting function with different values of parameters if other arguments of \code{f} are provided.
}
\value{
  A list with the following elements:
  \item{all }{A data.frame where each row correspond to one spatial weighint matrix tested. It contains value of parameteres tested and corrected AIC and number of orthogonal vectors for the best model.}
  \item{best }{A list containing results of scores.listw and ortho.AIC of the best spatial weighting matrix according to corrected AIC.}
  
}
\references{ Dray S., Legendre P. and Peres-Neto P. R. (submitted)
  Spatial modeling: a comprehensive framework for principal coordinate
  analysis of neighbor matrices (PCNM). \emph{Ecological Modelling}.
 }
\author{ Stephane Dray }
\seealso{  \code{\link{ortho.AIC}}, \code{\link{scores.listw}} }
\examples{

\dontrun{
data(oribatid)
# Hellinger transformation
fau <- sqrt(oribatid$fau/outer(apply(oribatid$fau,1,sum),rep(1,ncol(oribatid$fau)),"*"))
# remove gradient effect
faudt <- resid(lm(as.matrix(fau)~as.matrix(oribatid$xy)))

# test a binary spatial weighting matrix
nbtri <- tri2nb(as.matrix(oribatid$xy))
tri.res <- test.W(faudt,nbtri) 

maxi <- max(unlist(nbdists(nbtri,as.matrix(oribatid$xy))))

# test a simple spatial weighting function of the distance
f1 <- function(x) {1-(x)/(maxi)}
tri.f1 <- test.W(faudt,nbtri,f=f1,xy=as.matrix(oribatid$xy)) 

# test a spatial weighting function with various values of parameters
f2 <- function(x,dmax,y) {1-(x^y)/(dmax)^y}
tri.f2 <- test.W(faudt,nbtri, f=f2,y=2:10,dmax=maxi,xy=as.matrix(oribatid$xy)) 





}   


}
\keyword{ spatial }
