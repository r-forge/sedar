\name{nbmc}
\alias{nbmc}
\title{Neighbourhood list based on metacommunity models}
\description{
The function uses grouping vectors or factors for regions and habitat classes to make neighbourhood list that can be used e.g. in the 'prabclus' package and in the \code{\link{permatsar}} function.
}
\usage{
nbmc(loc, reg, hab = NULL, process = c("dl", "ss", "ls", "me", "zd", "ud"), 
cont = rep(1, length(loc)))
}
\arguments{
  \item{loc}{a vector containing names of the locations.}
  \item{reg}{grouping vector for regions, factor levels or unique values are used.}
  \item{hab}{grouping vector for habitat classes, factor levels or unique values are used. If \code{NULL}, one level is used for as \code{rep(1, length(loc))}.}
  \item{process}{character, choose the metacommunity model to use: \code{dl} dispersal limitation, \code{ss}species sorting, \code{ls} dispersal limitation and species sorting, \code{me} mass effect, \code{zd} zero dispersal, \code{ud} unlimited dispersal.}
  \item{cont}{grouping vector for continents, factor levels or unique values are used. Continents are larger areas, and above that scale, dispersal is not possible.}
}
\value{
A list with neighborhood(s) of each element of \code{loc}.
}
\author{Peter Solymos, \email{Solymos.Peter@aotk.szie.hu}, \url{http://www.univet.hu/users/psolymos/personal/}}
\examples{
loc <- 1:8
reg <- c(1,1,1,1,2,2,2,2)
hab <- c(1,1,2,2,1,1,2,2)

nbmc(loc,reg,hab,"dl")
nbmc(loc,reg,hab,"ss")
nbmc(loc,reg,hab,"me")
}
\keyword{ manip }
