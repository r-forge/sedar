\name{hubbell}
\alias{hubbell}
\alias{plot.hubbell}
\alias{print.hubbell}
\alias{hubb2vect}

\title{ Plays a Hubbell game. }
\description{
  Hubbell game is a zero-sum game, where a number of individuals is
  killed during a timestep, and these are replaced either from within the
  remaining community or from the surrounding metacommunity. In either
  case, the replacement probabilities are given by species frequencies.
}
\usage{
hubbell(comm, D, m=0, P=NULL)
hubb2vect(vec)
\method{plot}{hubbell}(x, sp.max=64, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{comm, x}{A \code{hubbell} object giving the counts of each species
    in a community. }
  \item{D}{Number of individuals killed during time step. }
  \item{m}{Probability of immigration. }
  \item{P}{Relative frequencies of species in the metacommunity. }
  \item{vec}{A vector giving counts of each species.}
  \item{sp.max}{Species number parameter for plots.}
  \item{...}{Other parameters passed to functions.}
}
\details{
  Hubbell game is zero-sum game for a local community.  In each time
  step, \code{D} individuals are killed.  These may be replaced either
  by the remaining species in the local community, or if \code{m}
  \eqn{> 0}, they may be replaced by the surrounding metacommunity.
  Metacommunity frequencies \code{P} must be given by the user -- a
  Hubbellian way is to use \code{\link{rhubbell}}.

  The community \code{comm} is an object of class \code{hubbell}.  This
  is simply a vector giving the species id. number for each individual
  in the community.  Function \code{hubb2vect} can be used to produce
  the initial community.  Function 
  \code{plot} tries to displays the canopy using colour codes for
  species, together with a ranked abundance -- diversity (RAD)  curve.  The
  parameter \code{sp.max} determines the number of colours in plot and
  the scaling of horizontal axis in RAD plot.
}
\value{
  A \code{hubbell} object.
}
\references{
  Hubbell, S.P. (2001). The Unified Neutral Theory of Biodiversity and
  Biogeography.  Princeton Univ. Press.
}
\author{ Jari Oksanen }

\seealso{ \code{\link[vegan]{radfit}} for alternative visualization and
  analysis.  }

\examples{
## Island Biogeography: completely isolated community 
## starting with 40 species, each with 30 individuals (J=1200)
## and going for 100 time steps. Kill 240 individuals each year.
comm <- hubb2vect(rep(40,30))
comm
plot(comm, sp.max=30)
comm <- hubbell(comm, 240)
plot(comm, sp.max=30)
## Look for 100 time steps ahead
for (t in 1:100) comm <- hubbell(comm, D=240)
plot(comm, sp.max=30)
## Make a Hubbell metacommunity, and start immigration
metacomm <- rhubbell(20, 16000)
for (t in 1:100) comm <- hubbell(comm, D=240, m=0.1, P=metacomm)
comm  ## More species now
plot(comm)
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{models}
\keyword{ distribution }


