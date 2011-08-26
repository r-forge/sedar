\encoding{UTF-8}
\name{metacycle}
\alias{metacycle}
\alias{metastep}

\title{ Simulate a Metapopulation System }
\description{
  Function \code{metacycle} runs simulation of a metapopulation
  system. 
}
\usage{
metacycle(steps = 1, d, A, p = rep(1, length(A)), y = 1, x = 1, e = min(A), alpha = 1, locations = NULL, ...)
}

\arguments{
  \item{steps}{ Number of time steps }
  \item{d}{ Distances among patches  }
  \item{A}{ Areas of patches }
  \item{p}{ Initial occupancies of the patches, must be 0 or 1. }
  \item{y, x, e}{ Parameters of the metapopulation models }
  \item{alpha}{ Exponential decay rate in connectivity }
  \item{locations}{ Coordinates of the patches }
  \item{\dots}{ Other parameters of the model }
}
\details{
 For details of the model see the pdf document in the \code{doc} folder
 of the packages: Chapter 4 explains simulation, and Chapter 3 describes
 the estimeatin of the parameters.

 Function \code{metastep} performs one simulation step. It is a lower
 level function and trickier to parametrize. See the code for its
 working, or consult the pdf file in \code{doc} folder.
}
\value{
  Function \code{metacycle} returns an object of class
 \code{"metacycle"} with following items. 
  \item{p }{Occupancy matrix with a row for each location and column
    for each step}
  \item{d }{ connectivity from the exponenential decay model.}
  \item{A}{ Sizes of patches}
  \item{y, x, e, alpha}{Input parameters of the model}
  \item{locations}{Coordinates of the patches, found from
    \code{\link{cmdscale}} if not given. }
  \item{J.obs}{Observed incidences for each patch.}
  \item{P.obs}{Observed number of occupied site for each step.}
  \item{J.pot, S.pot, C.pot}{Expected incidence and parameters \eqn{S}
    and \eqn{C}.}
}
\references{
Hanski, I. 1999. Metapopulation Ecology. Oxford UP.
}
\author{ Jari Oksanen }

\seealso{ nothing yet }
\examples{
## Parameters similar to a Melitae cinxia network
data(fritty)
d <- with(fritty, dist(cbind(x.crd, y.crd)))
sim <- with(fritty, metacycle(100, d, A, p, x = 0.42, e = 0.061, y = 15.2))
plot(sim$P.obs,   type = "l", col = "blue", lwd = 2,
 xlab = "Time", ylab = "Population Size")
plot(sim$J.pot, sim$J.obs, xlab="Expected incidence", ylab="Simulated incidence")
}

\keyword{ models }

