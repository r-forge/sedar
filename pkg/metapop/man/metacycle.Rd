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
%- maybe also 'usage' for other objects documented here.
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
  \item{pmat }{Occupancy matrix with a row for each location and column
    for each step}
  \item{d }{ connectivity from the exponenential decay model.}
  \item{A}{ Sizes of patches}
  \item{y, x, e, alpha}{Input parameters of the model}
  \item{locations}{Coordinates of the patches}
}
\references{
Hanski, I. 1999. Metapopulation Ecology. Oxford UP.
}
\author{ Jari Oksanen }
\note{
  Awaits documetations
}
\seealso{ nothing yet }
\examples{
## Nothing yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

