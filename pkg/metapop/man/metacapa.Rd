\name{metacapa}
\alias{metacapa}
\alias{plot.metacapa}

\title{Metapopulation Capacity }
\description{
Metapopulation capacity is the ability of the landscape to sustain a
metapopulation.  
}
\usage{
metacapa(d, A, alpha = 1, ...)
\method{plot}{metacapa}(x, locations = NULL,  ...)
}

\arguments{
  \item{d}{ Distances among patches }
  \item{A}{ Sizes of patches }
  \item{alpha}{ Decay parameter of the exponential connectivity equation}
  \item{x}{\code{metacapa} result object}
  \item{locations}{Locations of patches. If missing, the locations are
    found from metric scaling of distances.}
  \item{\dots}{ =ther parameters to the functions }
}
\details{
  Metapopulation capacity is defined as the largest eigenvalue of matrix
  with non-diagonal elements \eqn{A_i \exp(-\alpha d_{ij} A_j} and zero
  diagonal. 
}
\value{
Function returns a list of class \code{"metacapa"} with following items:
  \item{capacity }{Metapopulation capacity}
  \item{x }{contributions of patches to the metapopulation capacity}
  \item{d}{distances among patches}
  \item{A}{size of patches}
}
\references{
Hanski, I. & Ovaskainen, O. 2000. The metapopulation capacity of
fragmented landscape. \emph{Nature} 404: 755--758.
}
\author{ Jari Oksanen }

\seealso{ nothing yet }
\examples{
data(fritty)
attach(fritty)
d <- dist(cbind(x.crd, y.crd))
sol <- metacapa(d, A)
plot(sol)
}

\keyword{ models }

