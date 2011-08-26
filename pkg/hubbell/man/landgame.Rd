\name{landgame}
\alias{landgame}
\alias{makeland}
\alias{mutate}
\alias{plot.landscape}
\alias{as.matrix.landscape}

\title{Build a Metacommunity from scratch. }
\description{
  Function defines a Metacommunity as a network of local
  communities. Each local community follows its own dynamics, plants may
  immigrate to the local community, and at birth, they may evolve to a
  new species.  
}

\usage{
landgame(land, D, m, theta)
makeland(nrow, ncol, J)
mutate(genes)
\method{plot}{landscape}(x,main, ...)
\method{as.matrix}{landscape}(x, ...)
}

\arguments{
  \item{land, x}{Object of class `landscape', giving the co-ordinates local
  communities and identities of each \eqn{J} plants of the local community. }
  \item{D}{Number of plants killed in a local community per timestep.}
  \item{m}{Probability of immigration from neighbouring local communities. }
  \item{theta}{Hubbell's fundamental biodiversity number \eqn{\theta}.}
  \item{nrow}{Number of rows of local communities in a Metacommunity.}
  \item{ncol}{Number of columns of local communities in a
    Metacommunity.}
  \item{J}{Size of each local community.}
  \item{genes}{Genes defining the species, as a hexidecimal RGB colour.}
  \item{main}{Main title for the plot.}
  \item{...}{Other parameters to the function.}
}
\details{
  Hubbell's Metacommunity game builds a metacommunity from scratch.
  Metacommunity consists of \code{nrow} rows and \code{ncol} columns of
  local communities, each with \code{J} plants. Function \code{makeland}
  creates a one-species Metacommunity.  The steps of one run are:
  \enumerate{
    \item
    Kill \code{D} individuals in each local community.
    \item
    Fill the vacated slots with immigrants from immediate siding
    neighbours at probability \code{m}. The species identities of the
    immigrants are randomly chosen from the surrounding species pool.
    \item
    Fill the remaining free slots with offspring of the remaining local
    community, again using probabilities proportional to
    post-disturbance abundances.
    \item
    Mutate new \code{D} individuals to a different species at the
    probability \eqn{\nu} derived from the Ultimate Diversity
    \eqn{\theta} and Metacommunity size \eqn{J_M}: \eqn{\nu =
      \theta/(2J_M)}.  Species are identified by RGB colours, and
    function \code{mutate} changes one of the RGB components into a random
    value, without checking whether the value really changed from the
    old one, or whether this species already exists. 
  }

  Function \code{plot} does not plot the Original Species, or the species
  with identity \code{#FFFFFF}, and other species are displayed by their
  colour.  For better visibility, it may be better to use grey
  background (see \code{bg} in \code{\link{par}}); grey species are late
  in evolution.

  Function \code{as.matrix} returns the \code{landgame} result as an
  ordinary community data matrix where values are counts of species in
  local communities. The species are named after their
  \code{\link{rgb}} colour codes.  
}

\value{
  Function returns an object of class \code{landscape}, which is simply
  a three dimensional array.
}
\references{
  Hubbell, S.P. (2001). The Unified Neutral Theory of Biodiversity and
  Biogeography.  Princeton Univ. Press.
}
\author{ Jari Oksanen}
\note{
  Unlike Hubbell's original program, this function uses a torus
  landscape so that opposite margins join. In this way it is possible to
  avoid edge effects and whole landscape can be inspected without buffer
  zone.

  The function is slow for large problems. For instance, Hubbell
  analysed a landscape of 101 \eqn{\times} 101 local communities, each
  with \eqn{J=16} and run 100,000 time steps.  With my crummy, old
  laptop (i586, 233MHz) this would take 19 days with isolated
  communities (\eqn{m=0.005)} and 38 days with open communities
  (\eqn{m=0.5}). 
}


\seealso{\code{\link{hubbell}}.}

\examples{
land <- makeland(21, 21, 16)
land
op <- par(bg="gray")
plot(land)
for (t in 1:20) land <- landgame(land, D=4, m=0.2, theta=14)
land
plot(land)
par(op)
}
\keyword{ models }

