\name{PCNM}
\alias{PCNM}
\title{ PCNM spatial eigenfunctions }
\description{
Function \code{\link{PCNM}} computes principal coordinates of neighbour matrices (PCNM) from a geographic distance matrix, in view of spatial eigenfunction analysis.

This version allows users to obtain either the classical PCNM eigenfunctions with positive eigenvalues only, or all PCNM eigenfunctions with positive and negative eigenvalues.
}
\usage{
PCNM(matdist, thresh=NULL, all=FALSE, include.zero=FALSE)
}

\arguments{
  \item{matdist}{ A distance matrix of class \code{dist} or \code{matrix}. }
  \item{thresh}{ A threshold value for truncation of the geographic distance matrix. }
  \item{all}{ If TRUE, the eigenvectors corresponding to all eigenvalues, positive and negative, are shown in the output list. Default value: \code{all=FALSE}. }
  \item{include.zero}{ If FALSE (default value), the zero eigenvalues as well as their eigenvectors are excluded from the output list. See details. }
}

\details{
If a truncation value is not provided, the largest distance in a minimum spanning tree linking all objects on the map is computed. That value is used as the truncation value.

When \code{all} is \code{FALSE}, the principal coordinate norm is the square root of the corresponding eigenvalue, as in Gower's principal coordinate analysis. When \code{all} is \code{TRUE}, the principal coordinates are normalized to length 1 to avoid creating complex principal coordinate axes for the negative eigenvalues.

When \code{include.zero} is \code{TRUE}, the null eigenvalue(s) and associated PCNM eigenfunction(s) are shown in the output list. These eigenfunctions are meaningless and should not be used for spatial modelling.
}

\value{

  \item{PCNM_eigenvalues }{The eigenvalues. }
  \item{rel.values }{The relative eigenvalues. }
  \item{rel.cum.values }{The cumulative relative eigenvalues. }
  \item{PCNM }{The PCNM eigenfunctions (principal coordinates of the truncated distance matrix). }
  \item{trace }{The trace of the truncated distance matrix. This is also the sum of all eigenvalues, positive and negative. }
}

\references{ 
Borcard, D. and P. Legendre. 1994. Environmental control and spatial structure in ecological communities: an example using Oribatid mites (Acari, Oribatei). Environmental and Ecological Statistics 1: 37-53.
     
Borcard, D. and P. Legendre. 2002. All-scale spatial analysis of ecological data by means of principal coordinates of neighbour matrices. Ecological Modelling 153: 51-68.

Borcard, D., P. Legendre, C. Avois-Jacquet and H. Tuomisto. 2004. Dissecting the spatial structure of ecological data at multiple scales. Ecology 85: 1826-1832 [Statistical report].

Dray, S., P. Legendre and P. R. Peres-Neto. 2006. Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). Ecological Modelling 196: 483-493.
}

\author{ Stephane Dray, Universite Lyon I, Pierre Legendre and Daniel Borcard, Universite de Montreal }

\seealso{\code{\link{pcoa.all}}, \code{\link{quickPCNM}} }

\examples{
# Oribatid mite data from Borcard and Legendre (1994)
library(vegan)
data(mite)      # 70 peat cores, 35 species
data(mite.xy)   # Core geographic coordinates

mite.D <- dist(mite.xy)   # Geographic distance matrix
pcnm1 <- PCNM(mite.D, thresh=1.012)
pcnm2 <- PCNM(mite.D, thresh=1.012, all=TRUE)
pcnm3 <- PCNM(mite.D, thresh=1.012, all=TRUE, include.zero=TRUE)

# Plot maps of the first 3 PCNM eigenfunctions of file pcnm1
par(mfrow=c(1,3))
s.value(mite.xy, pcnm1$vectors[,1], clegend=0, sub="PCNM 1", csub=2.5)
s.value(mite.xy, pcnm1$vectors[,2], clegend=0, sub="PCNM 2", csub=2.5)
s.value(mite.xy, pcnm1$vectors[,3], clegend=0, sub="PCNM 3", csub=2.5)

# Plot maps of PCNM eigenfunctions no. 1 (eigenvalue = 77.12495), 
# no. 44 (eigenvalue = 0), and no. 70 (eigenvalue = -18.92951) of file pcnm3
par(mfrow=c(1,3))
s.value(mite.xy, pcnm3$vectors[,1] , clegend=0, sub="PCNM 1" , csub=2.5)
s.value(mite.xy, pcnm3$vectors[,44], clegend=0, sub="PCNM 44", csub=2.5)
s.value(mite.xy, pcnm3$vectors[,70], clegend=0, sub="PCNM 70", csub=2.5)
}

\keyword{ multivariate }
\keyword{ spatial }
