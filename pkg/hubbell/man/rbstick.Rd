\name{rbstick}
\alias{rbstick}

\title{Generates a Random Brokenstick Community }

\description{Generates a random brokenstick community.}

\usage{
rbstick(n)
}

\arguments{
  \item{n}{Number of species. }
}
\details{
  Brokenstick has been a popular null-model for abundance relations
  since the days of MacArthur.  Brokenstick community is generated so
  that species simultaneously divide the available community (or resources) in
  random pieces.  When these pieces are arranged in rank
  abundance order, the result looks very similar to log-Normal.

  The current functin returns a random community. Expected brokenstick
  communities can be generated using \pkg{vegan} function
  \code{\link[vegan]{bstick}}. 
}
\value{
  Returns the ranked abundance proportions (sum = 1).
}
\references{
  Tokeshi, M. (1993). Species abundance patterns and community
  structure. Adv. Ecol. Res. 24, 111-186.

  Wilson, J.B. (1993). Would we recognize a Broken-Stick community if we
  found one? Oikos 67: 181-183.
}
\author{Jari Oksanen }

\seealso{ \code{\link[vegan]{radfit}} for displaying results and fitting
  a brokenstick model, and \code{\link[vegan]{bstick}} for the expected
  brokenstick distribution. Function
  \code{\link{rhubbell}} gives a Hubbellian alternative.
}

\examples{
## Generate and plot a broken-stick community
comm <- rbstick(30)
require(vegan)
(mod <- radfit(comm, family=gaussian))
plot(mod) 
}

\keyword{ distribution }




