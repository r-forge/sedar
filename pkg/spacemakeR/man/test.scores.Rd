\name{test.scores}
\alias{test.scores}
\title{ Function to compute and test Moran's I for eigenvectors of spatial weighting matrices.}
\description{
  This function tests Moran's I for each eigenvector of a spatial
  weighting matrix. It uses the permutation procedure implemented in the
  \code{spdep} package.
}
\usage{
test.scores(scores,listw,nsim=999)
}

\arguments{
  \item{scores}{ An object created by the score.listw function. }
  \item{listw}{ An object of the class \code{listw} created by functions of the \code{spdep} package. }
  \item{nsim}{ Number of permutations.}
}
\details{Negative autocorrelation is tested for eigenvectors associated
  with negative eigenvalues. Positive autocorrelation is tested for eigenvectors associated
  with positive eigenvalues.}
\value{
  A \code{data.frame} with values of Moran's I and p-value
  
}

\author{ Stephane Dray }


\keyword{ spatial }
