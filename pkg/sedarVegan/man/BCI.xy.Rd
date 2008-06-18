\name{BCI.xy}
\alias{BCI.xy}
\docType{data}
\title{ Spatial Coordinates of the Barro Colorado Rainforest Plots}
\description{
  UTM coordinates of 50 1-ha sample plots in the Barro Colorado
  Island, Panama.
}
\usage{data(BCI.xy)}
\format{
  A data frame with 50 observations on the following 2 variables.
  \describe{
    \item{\code{UTM.EW, UTM.NS}}{UTM coordinates}
  }
}
\details{
  The tree stem count data are given in the data set
  \code{\link[vegan]{BCI}} in \pkg{vegan}. 
}
\source{
  \url{http://www.sciencemag.org/cgi/content/full/295/5555/666/DC1}
}
\references{
   Condit, R, Pitman, N, Leigh, E.G., Chave, J., Terborgh, J., Foster,
  R.B., \enc{Nuñez}{Nunez}, P., Aguilar, S., Valencia, R., Villa, G., Muller-Landau,
  H.C., Losos, E. & Hubbell, S.P. (2002). Beta-diversity in tropical
  forest trees. \emph{Science} 295, 666--669.
}
\seealso{Data sets \code{\link[vegan]{BCI}} (\pkg{vegan}),
\code{\link[BiodiversityR]{BCI.env}} (\pkg{BiodiversityR}).}
\examples{
data(BCI.xy)
library(vegan)
data(BCI)
## Unconstrained ordination
ord <- cca(BCI)
ordisurf(BCI.xy, scores(ord, display="sites", choices=1),
  main="Unonstrained CA")
## Multiscale ordination and variogram
mord <- mso(ord, BCI.xy)
msoplot(mord)
## Constrain by pcnm (the version in sedarVegan)
nm <- sedarVegan::pcnm(BCI.xy)$vectors
## Take half of the pcnm vectors
take <- round(ncol(nm)/2)
ord <- vegan::cca(BCI, nm[,1:take])
## Use LC scores in plotting
ordisurf(BCI.xy, scores(ord, display="lc", choices=1), 
  main = "Constrained by PCNM")
msoplot(mso(ord, BCI.xy))
}
\keyword{datasets}
