\encoding{UTF-8}
\name{chooseCN}
\alias{chooseCN}
\alias{plot.chooseCN}

\title{Function to choose a connection network}

\description{The function \code{chooseCN} is a simple interface to
  build a connection network (CN) from xy coordinates. The user chooses
  from seven types of graph. \code{chooseCN} calls functions from
  appropriate packages, handles non-unique coordinates and returns a
  connection network either with class \code{nb} or \code{listw} defined
  in \pkg{spdep} package.  
}
\usage{
chooseCN(xy, type = "delauney", result.type = "nb", d1 = NULL, 
    d2 = NULL, k = NULL,  a=NULL, dmin=NULL, support = c("vegan", "ade4")) 
}
\arguments{
  \item{xy}{an matrix or data.frame with two columns for x and y coordinates.}
  \item{type}{a text giving the type of graph (see details).}
  \item{result.type}{a character giving the class of the returned
    object. Either "nb" (default) or "listw", both from \pkg{spdep}
    package. See details.}
  \item{d1}{the minimum distance between any two neighbours. Used if
    \code{type "neighb.by.dist".}}
  \item{d2}{the maximum distance between any two neighbours. Used if
    \code{type = "neighb.by.dist"}.}
  \item{k}{the number of neighbours per point. Used if
    \code{type = "knearest"}.}
  \item{a}{the exponent of the inverse distance matrix. Used if
    \code{type = "inv.dist"}.}
  \item{dmin}{the minimum distance between any two distinct points. Used
    to avoid infinite spatial proximities (defined as the inversed
    spatial distances). Used if \code{type = "inv.dist"}.}
  \item{support}{Prefer \pkg{vegan} or \pkg{ade4} functions; only
    concerns \code{type = "mst"}.} 
}
\details{There are 7 kinds of graphs proposed: \cr
  Delaunay triangulation (type \code{"delauney"})\cr
  Gabriel graph (type \code{"gabriel"})\cr
  Relative neighbours (type \code{"rel.neighbours"})\cr
  Minimum spanning tree (type \code{"mst"})\cr
  Neighbourhood by distance (type \code{"neighb.by.dist"})\cr
  K nearests neighbours (type \code{"knearest"})\cr
  Inverse distances (type \code{"inv.dist"})

  The last option (\code{type = "inv.dist"}) is not a true neighbouring
  graph: all sites are neighbours, but the spatial weights are directly
  proportional to the inversed spatial distances.  Also note that in
  this case, the output of the function is always a \code{listw} object,
  even if \code{nb} was requested.


}
\value{Returns a connection network having the class \code{nb} or
  \code{listw}. The xy coordinates are passed as attribute to the
  created object.
}
\seealso{\code{\link[adegenet]{chooseCN}}
}
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr},
  veganified by Jari Oksanen. }

\note{This function was originally part of the \pkg{adegenet} package
  for R, and then ported to \code{sedar} project in
  \url{http://sedar.r-forge.r-project.org/}.  Jari Oksanen took the code
  from Thibaut Jombart's working package to see how to change the
  function for standard \pkg{vegan} style. The major changes were that
  now \code{chooseCN} works silently without asking anything from the
  user, and that the \code{plot.chooseCN} function is separate from the
  calculation. }

\examples{
if(require(spdep) && require(ade4)){

xy <- matrix(runif(100),ncol=2)
op <- par(mfrow=c(2,2))
plot(cn1 <- chooseCN(xy, type="del"))
plot(cn2 <- chooseCN(xy, type="gab"))
plot(cn3 <- chooseCN(xy, type="mst"))
plot(cn4 <- chooseCN(xy, type="nei", d1=0, d2=0.2))
par(op)
}
}
\keyword{spatial}
\keyword{utilities}