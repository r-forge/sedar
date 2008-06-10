\name{rhubbell}
\alias{rhubbell}

\title{Simulates a community from Hubbell's species generator. }
\description{
  Simulates a community using Hubbell's species generator.
}
\usage{
rhubbell(theta, J)
}

\arguments{
  \item{theta}{Hubbell's ultimate diversity number \eqn{\theta}. }
  \item{J}{Size of generated community. }
}
\details{
  Hubbell's species generator \eqn{P = \theta/(\theta + j -1)} gives the
  probability that \eqn{j}th collected individual belongs to a new
  species.  If, by chance,
  this happens, a new species is added to the community with abundance
  1.  If the \eqn{j}th individual belongs to an old species, abundance 1 is
  added to some of the old species.  The probability of getting this new
  indivdiual is proportional to species abundances in community at the
  moment. 
}
\value{
  Returns a non-ranked vector of species abundances (integers).
}
\references{
  Hubbell, S.P. (2001). The Unified Neutral Theory of Biodiversity and
  Biogeography.  Princeton Univ. Press.
}
\author{ Jari Oksanen }
\note{ The \R code is crystal clear -- at least much clearer than
  Hubbell's description (p. 290--291, incomprehensible spaghetti diagram
  in Fig. 9.1).  Use the force, read the source.
}


\seealso{
  \code{\link[vegan]{radfit}} for displaying results,
  \code{\link{rbstick}} for an alternative model (but giving very
  similar results).
}

\examples{
## A nemoral forest
comm <- rhubbell(20, 200)
require(vegan)
plot(radfit(comm))
}
\keyword{ distribution }

