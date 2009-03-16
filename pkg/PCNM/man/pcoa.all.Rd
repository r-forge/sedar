\name{pcoa.all}
\alias{pcoa.all}
\title{ Principal Coordinate Analysis }
\description{
Function \code{\link{pcoa.all}} computes principal coordinate decomposition (also called classical scaling) of a distance matrix D (Gower 1966) and outputs the eigenvectors corresponding to all eigenvalues, positive and negative.
}
\usage{
pcoa.all(D, diagonal=0, all=FALSE, include.zero=FALSE, rn=NULL)
}

\arguments{
  \item{D}{ A distance matrix of class \code{dist} or \code{matrix}. }
  \item{diagonal}{ In an ordinary dissimilarity or distance matrix, the diagonal is 0. Other constant values can be forced on the diagonal of matrix D. For example, for the construction of db-MEM eigenfunctions, the diagonal contains 4*thresh. }
  \item{all}{ If TRUE, the eigenvectors corresponding to all eigenvalues, positive and negative, are shown in the output list. Default value: \code{all=FALSE}. }
  \item{include.zero}{ If FALSE (default value), the zero eigenvalues as well as their eigenvectors are excluded from the output list. See details. }
  \item{rn}{ An optional vector of row names, of length n, for the n objects. }
}

\details{
This function was developed to provide principal coordinate decomposition for the function \code{PCNM} of this library. It computes the eigenvectors corresponding to all eigenvalues, positive and negative, using matrix algebra. It is faster than \code{cmdscale}, which only outputs the positive eigenvalues and corresponding eigenvectors.

Another function in this library, \code{pcoa}, offers more options for principal coordinate analysis of ordinary distance matrices computed from community composition or other response data tables. Its companion for the production or ordination diagrams is \code{biplot.pcoa}.

In \code{pcoa.all}, when there are no negative eigenvalues, the eigenvectors (principal coordinates) are scaled to a length equal to the square root of the corresponding  eigenvalues. When negative eigenvalues are present, the eigenvectors are scaled to length 1 to prevent the production of complex engenvectors. Eigenvectors scaled in that way do not preserve the distances in matrix \code{D} among the objects. They should not be used, in particular, to draw ordination graphs.

When \code{include.zero} is \code{TRUE}, the null eigenvalue(s) and associated eigenvector(s) are shown in the output list. These eigenvectors are meaningless and should not be used for modelling.
}

\value{

  \item{values }{The eigenvalues. }
  \item{rel.values }{The relative eigenvalues. }
  \item{rel.cum.values }{The cumulative relative eigenvalues. }
  \item{vectors }{The principal coordinates. }
  \item{trace }{The trace of the distance matrix. This is also the sum of all eigenvalues, positive and negative. }
}

\references{ 
Gower, J. C. 1966. Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika. 53: 325-338.
     
Legendre, P. and L. Legendre. 1998. Numerical ecology, 2nd English edition. Elsevier Science BV, Amsterdam.
}

\author{ Pierre Legendre, Universite de Montreal }

\seealso{\code{\link{PCNM}}, \code{\link{pcoa}} }

\examples{
# Oribatid mite data from Borcard and Legendre (1994)
library(vegan)
data(mite)      # 70 peat cores, 35 species

mite.D <- vegdist(mite, "bray")
res <- pcoa.all(mite.D)
res <- pcoa.all(mite.D, all=TRUE)
}

\keyword{ multivariate }
