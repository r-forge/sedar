\name{PCNM}
\alias{PCNM}
\title{PCNM spatial eigenfunctions}
\description{
Function \code{\link{PCNM}} computes principal coordinates of neighbour matrices (PCNM) from a geographic distance matrix, in view of spatial eigenfunction analysis.

This version allows users to obtain either the classical PCNM eigenfunctions with positive eigenvalues only, or all PCNM eigenfunctions with positive and negative eigenvalues.
}
\usage{
PCNM(matdist, thresh=NULL, dbMEM=FALSE, moran=NULL, all=FALSE, include.zero=FALSE)
}

\arguments{
  \item{matdist}{ A distance matrix of class \code{dist} or \code{matrix}. }
  \item{thresh}{ A threshold value for truncation of the geographic distance matrix. If \code{thresh=NULL}, the length of the longest edge of the minimum spanning tree will be used as the threshold. }
  \item{dbMEM}{ If FALSE, the diagonal of matrix D contains zeros; classical PCNM eigenfunctions are produced. If TRUE, the diagonal of D contains 4*thresh; the eigenfunctions produced are called distance-based MEM or db-MEM (Dray et al. 2006). db-MEM eigenfunctions with positive eigenvalues all have positive Moran's I. The db-MEMs with positive eigenvalues are fewer in number than the PCNMs with positive eigenvalues. }
  \item{moran}{ If TRUE, Moran's I are computed for all PCNMs. If FALSE, Moran's I are not computed. When \code{moran=NULL} (default value), Moran's I are computed for PCNM but not for dbMEM eigenfunctions because, for dbMEM, Moran's I are proportional to the eigenvalues. That is not the case for classical PCNMs. }
  \item{all}{ If TRUE, the eigenvectors corresponding to all eigenvalues, positive and negative, are shown in the output list. Default value: \code{all=FALSE}. }
  \item{include.zero}{ If FALSE (default value), the zero eigenvalues as well as their eigenvectors are excluded from the output list. See details. }
}

\details{
If a truncation value is not provided, the largest distance in a minimum spanning tree linking all sites on the map is computed. That value is used as the truncation value.

If the calculation of Moran's I is requested but no file of sites coordinates is provided, cartesian coordinates are recomputed by principal coordinate analysis from the distance matrix \code{matdist}. 

When \code{all} is \code{FALSE}, the principal coordinate norm is the square root of the corresponding eigenvalue, as in Gower's principal coordinate analysis. When \code{all} is \code{TRUE}, the principal coordinates are normalized to length 1 to avoid creating complex principal coordinate axes for the negative eigenvalues.

When \code{include.zero} is \code{TRUE}, the null eigenvalue(s) and associated PCNM eigenfunction(s) are shown in the output list. These eigenfunctions are meaningless and should not be used for spatial modelling.
}

\value{

  \item{values }{The PCNM eigenvalues. }
  \item{vectors }{The PCNM or db-MEM eigenfunctions (principal coordinates of the truncated distance matrix). }
  \item{Moran_I }{A 3-column table giving: \code{Moran} = Moran's I statistics, \code{p.value} = p-values (2-tailed parametric test), \code{Positive} = 1 for Moran's I larger than the expected value. In most instances, only the PCNM eigenfunctions displaying positive spatial correlation (i.e., those with Moran's I larger than the expected value) are used in PCNM modelling. Computed using function \code{moran.I.multi} of Guillaume Blanchet's AEM package. }
  \item{expected_Moran }{The expected value of Moran's I. }
  \item{spanning }{An R object containing the minimum spanning tree. }
  \item{thresh }{The truncation threshold. If a \code{thresh} value was provided by the user, that value is given here. If \code{thresh=NULL}, this is the length of the longest edge in the minimum spanning tree; the threshold is augmented by 0.000001 before printing to correct the cases where the threshold is rounded to a value slightly smaller than the true threshold. }
}

\author{ Stephane Dray, Universite Lyon I, Pierre Legendre and Daniel Borcard, Universite de Montreal, Guillaume Blanchet, University of Alberta }

\references{ 
Borcard, D. and P. Legendre. 1994. Environmental control and spatial structure in ecological communities: an example using Oribatid mites (Acari, Oribatei). Environmental and Ecological Statistics 1: 37-53.
     
Borcard, D. and P. Legendre. 2002. All-scale spatial analysis of ecological data by means of principal coordinates of neighbour matrices. Ecological Modelling 153: 51-68.

Borcard, D., P. Legendre, C. Avois-Jacquet and H. Tuomisto. 2004. Dissecting the spatial structure of ecological data at multiple scales. Ecology 85: 1826-1832 [Statistical report].

Dray, S., P. Legendre and P. R. Peres-Neto. 2006. Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). Ecological Modelling 196: 483-493.
}

\seealso{\code{\link{pcoa.all}}, \code{\link{quickPCNM}} }

\examples{
# Oribatid mite data from Borcard and Legendre (1994) using 
library(vegan)
data(mite)      # 70 peat cores, 35 species
data(mite.xy)   # Core geographic coordinates

mite.D <- dist(mite.xy)   # Geographic distance matrix

# thresh=1.012 is the value used in Borcard and Legendre (2002)
mite.pcnm1 <- PCNM(mite.D, thresh=1.012)
mite.pcnm2 <- PCNM(mite.D, thresh=1.012, all=TRUE)
mite.pcnm3 <- PCNM(mite.D, thresh=1.012, all=TRUE, include.zero=TRUE)

# Plot maps of the first 3 PCNM eigenfunctions of file mite.pcnm1
library(ade4)
par(mfrow=c(1,3))
s.value(mite.xy, mite.pcnm1$vectors[,1], clegend=0, sub="PCNM 1", csub=2.5)
s.value(mite.xy, mite.pcnm1$vectors[,2], clegend=0, sub="PCNM 2", csub=2.5)
s.value(mite.xy, mite.pcnm1$vectors[,3], clegend=0, sub="PCNM 3", csub=2.5)

# Plot maps of PCNM eigenfunctions no. 1 (eigenvalue = 77.12495), 
# no. 44 (eigenvalue = 0), and no. 70 (eigenvalue = -18.92951) of mite.pcnm3
library(ade4)
par(mfrow=c(1,3))
s.value(mite.xy, mite.pcnm3$vectors[,1] , clegend=0, sub="PCNM 1" , csub=2.5)
s.value(mite.xy, mite.pcnm3$vectors[,44], clegend=0, sub="PCNM 44", csub=2.5)
s.value(mite.xy, mite.pcnm3$vectors[,70], clegend=0, sub="PCNM 70", csub=2.5)

# Compute the db-MEM eigenfunctions. 
# By construction, they all have positive eigenvalues and Moran's I
mite.dbMEM1 <- PCNM(mite.D, dbMEM=TRUE)
dim(mite.dbMEM1$vectors)
plot(mite.dbMEM1$spanning, mite.xy)

# Compute all (n-1) db-MEM eigenfunctions, including those that 
# have negative eigenvalues and Moran's I modelling negative spatial correlation
mite.dbMEM2 <- PCNM(mite.D, dbMEM=TRUE, all=TRUE)
dim(mite.dbMEM2$vectors)
}

\keyword{ multivariate }
\keyword{ spatial }
