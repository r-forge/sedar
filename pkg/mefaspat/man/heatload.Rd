\encoding{UTF-8}
\name{heatload}
\alias{heatload}
\title{Heat Load Index}
\description{
Calculates heat load index based on aspect according to McCune and Keon (2002). Converts 0--360 degrees to 0--1 interval for heat load. It handles character vectors, too, for aspect notation.
}
\usage{
heatload(aspect, adjust = 45, flat = 0)
}
\arguments{
  \item{aspect}{a numeric or character vector. If numeric, it is a vactor with aspect values in degrees (0--360). If character, the first characters are used (not case sensitive) for major points of the compass or indication flat (\code{c("n", "e", "s", "w", "f")}).}
  \item{adjust}{numeric, in degrees. Degrees from south, where the index should have its maximum. The default 45 degrees means south-west.}
  \item{flat}{if \code{angle} is character, a flat area can be coded as \code{"flat"}. This value gives the angle in degrees which replaces \code{"flat"} values.}
}
\value{
A numeric vector.}
\references{
McCune, B. and D. Keon. 2002. Equations for potential annual direct incident radiation and heat load. \emph{Journal of Vegetation Science} 13, 603--606. \url{http://oregonstate.edu/~mccuneb/JVSreprint13.603-606.pdf}
}
\author{
\enc{Péter Sólymos}{Peter Solymos}, \email{solymos@ualberta.ca}
}
\examples{
x <- round(runif(25, 1, 360))
heatload(x)
plot(heatload(x), heatload(x, adjust = 0))
y <- sample(c("n", "e", "s", "w", "f"), 25, replace = TRUE)
heatload(y)
data.frame(y=y, f0=heatload(y, adjust = 0), f45=heatload(y, adjust = 0, flat = 45))
}
\keyword{ manip }
