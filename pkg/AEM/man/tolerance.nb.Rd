\name{tolerance.nb}
\alias{tolerance.nb}
\alias{find.angles}

\title{ Function to construct edges based on a tolerance angle and a maximum distance }
\description{
  This function creates an object of class \code{nb} (defined in the library spdep) containing a connexion diagram. The edges between sites are based on a tolerance angle and a maximum distance. The angle is directional; its direction is always from the bottow to the top of the screen.
}
\usage{
tolerance.nb(coords, unit.angle = "degrees", max.dist, tolerance, rot.angle, plot.sites = TRUE)
}

\arguments{
  \item{coords}{ A matrix or a data frame containing the X and Y coordinates of the study sites. }
  \item{unit.angle}{ Character. The measurement units in which angles are defined: either "degrees" (default) or "radians". }
  \item{max.dist}{ Numeric. The maximum distance of an edge linking two sites together. }
  \item{tolerance}{ Numeric. The tolerance angle in which a site can influence another site. The angle is measured vertically and from bottom to top of the pictures after rotation of the points.}
  \item{rot.angle}{ Numeric, optional. An angle at which a set of coordinates should be rotated before creating the connexion diagram. The set of coordinates is rotated counterclockwise. Negative values will produce a clockwise rotation.}
  \item{plot.sites}{ Logical (\code{TRUE}, \code{FALSE}) determining if the site should be plotted in a graphic window. This graph allows one to make sure the points are rotated in a correct direction.}
}
\details{
  Even though this function creates a connexion diagram based on a tolerance angle going from the bottom to the top of the screen, the resulting object is symmetric, meaning that a site influences another and vice versa. The final object does not represent a directional connexion network. It may be used to construct a site-by-edge matrix representing a directional process through the function \code{\link{build.binary}}
}
\value{
The function returns an object of class \code{nb} with a list of integer vectors  corresponding to neighbour region numbers.
}
\keyword{ tolerance, angle, coordinate }
\author{ F. Guillaume Blanchet }

\seealso{ \code{\link{knn2nb}}, \code{\link{dnearneigh}}, \code{\link{cell2nb}}, \code{\link{graph2nb}}, \code{\link{tri2nb}}, \code{\link{knn2nb}} }
\examples{

data(oribatid)

### Construct object of class nb with a tolerance angle of 30 degrees and a maximum distance of 2 m.
nb.oribatid<-tolerance.nb(oribatid$xy, unit.angle = "degrees", max.dist=2, tolerance = 30)

### Construct object of class nb with a tolerance angle of 30 degrees and a maximum distance of 2 m. The coordinates are rotated at an angle of 45 degrees counterclockwise.
nb.oribatid.2<-tolerance.nb(oribatid$xy, unit.angle = "degrees", max.dist=2, tolerance = 30, rot.angle = 45)

### Construct object of class nb with a tolerance angle of pi/8 radians and a maximum distance of 1.5 m. The coordinates are rotated at an angle of pi/4 radians clockwise.
nb.oribatid.3<-tolerance.nb(oribatid$xy, unit.angle = "radians", max.dist=1.5, tolerance = pi/8,rot.angle = -pi*2/3)


par(mfrow=c(1,3))
plot(nb.oribatid,oribatid$xy)
plot(nb.oribatid.2,oribatid$xy)
plot(nb.oribatid.3,oribatid$xy)
}

