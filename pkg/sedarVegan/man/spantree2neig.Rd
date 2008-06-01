\name{spantree2neig}
\alias{spantree2neig}

\title{ vegan spantree to ade4 neighbourhood translation }
\description{
 The function changes \code{\link[vegan]{spantree}} object of
 \pkg{vegan} into \code{\link[ade4]{neig}} of \pkg{ade4}. 
}
\usage{
spantree2neig(tree)
}
\arguments{
  \item{tree}{Result object of \pkg{vegan} \code{\link[vegan]{spantree}}.}
}
\details{
  The \code{\link[ade4]{neig}} objects can be further changed with
  \code{\link[ade4]{neig2nb}} to neighbours list objects of
  \pkg{spdep} package.
}
\value{
  Returns an object of class \code{\link[ade4]{neig}}. 
}

\author{ Jari Oksanen }

\seealso{ \code{\link[ade4]{neig}}, \code{\link[ade4]{neig2nb}},
\code{\link[spdep]{summary.nb}}, \code{\link[spdep]{plot.nb}}. }

\examples{
## vegan usage
require(vegan)
xy <- matrix(runif(100), nrow=50)
tr <- spantree(dist(xy))
plot(xy, asp=1)
lines(tr, xy)
## ade4 translation
require(ade4)
ntr <- spantree2neig(tr)
ntr
summary(ntr)
## spdep translation
require(spdep)
sntr <- neig2nb(ntr)
sntr
summary(sntr)
plot(sntr, coord=xy)
} 
\keyword{ multivariate }
\keyword{ spatial }
