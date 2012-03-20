\name{PCNM-package}
\alias{PCNM-package}
\docType{package}
\title{ PCNM spatial eigenfunctions }
\description{
Function \code{\link{PCNM}} computes principal coordinates of neighbour matrices (PCNM) and function \code{\link{quickPCNM}} carries out spatial eigenfunction analysis. Function \code{\link{pcoa}} computes principal coordinate decomposition of a distance matrix and implements two correction methods for negative eigenvalues. Functions \code{\link{STImodels}} and  \code{\link{quickSTI}} allow testing space-time ANOVA models in the absence of replication.
}
\details{
\tabular{ll}{
Package: \tab PCNM\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2009-05-27\cr
License: \tab Unlimited\cr
LazyLoad: \tab yes\cr
}
From a geographic distance matrix, function \code{\link{PCNM}} computes principal coordinates of neighbour matrices (PCNM) or Moran's eigenvector maps (MEM) based on distances, in view of spatial eigenfunction analysis.

Function \code{\link{quickPCNM}} automatically computes principal coordinates of neighbour matrices (PCNM) and carries out the spatial eigenfunction analysis.

Function \code{\link{pcoa.all}} computes principal coordinate decomposition (also called classical multidimensional scaling) of a distance matrix D and outputs the eigenvectors corresponding to all eigenvalues, positive and negative.

Function \code{\link{pcoa}} computes principal coordinate decomposition (also called classical multidimensional scaling) of a distance matrix D. It implements two correction methods for negative eigenvalues.

Function \code{\link{STImodels}} performs two-way ANOVA with test of the space-time interaction without replication using one among a set of possible models. In many models, degrees of freedom are saved by coding space and/or time parsimoniously using principal coordinates of neighbour matrices (\code{\link{PCNM}}).

Function \code{\link{quickSTI}} performs two-way ANOVA with test of the space-time interaction without replication in a simplified way.

}

\author{
Pierre Legendre, Daniel Borcard, F. Guillaume Blanchet, Miquel De Caceres and Stephane Dray

Maintainer: Pierre Legendre <pierre.legendre@umontreal.ca>
}

\references{
Borcard, D. and P. Legendre. 1994. Environmental control and spatial structure in ecological communities: an example using Oribatid mites (Acari, Oribatei). Environmental and Ecological Statistics 1: 37-53.
     
Borcard, D. and P. Legendre. 2002. All-scale spatial analysis of ecological data by means of principal coordinates of neighbour matrices. Ecological Modelling 153: 51-68.

Borcard, D., P. Legendre, C. Avois-Jacquet and H. Tuomisto. 2004. Dissecting the spatial structure of ecological data at multiple scales. Ecology 85: 1826-1832 [Statistical report].

Dray, S., P. Legendre and P. R. Peres-Neto. 2006. Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). Ecological Modelling 196: 483-493.

Legendre, P., M. De Caceres and  D. Borcard (Submitted to Ecology). Community surveys through space and time to assess environmental changes: testing space-time interaction in the absence of replication. 
}

\examples{

# Oribatid mite data from Borcard and Legendre (1994)
library(vegan)
data(mite)      # 70 peat cores, 35 species
data(mite.xy)   # Core geographic coordinates

# Compute PCNM eigenfunctions
res.pcnm <- PCNM(dist(mite.xy), thresh=1.012)

# PCNM eigenfunction analysis: quickPCNM
mite.hel <- decostand(mite, "hellinger")
res.quick <- quickPCNM(mite.hel, mite.xy, thresh=1.012)

# Principal coordinate ordination
mite.log <- log1p(mite[1:30,-35])
mite.D <- vegdist(mite.log, "bray")
res <- pcoa(mite.D)
mite.log.st = apply(mite.log, 2, scale, center=TRUE, scale=TRUE)
biplot(res$vectors, mite.log.st)
}

\keyword{ package }
