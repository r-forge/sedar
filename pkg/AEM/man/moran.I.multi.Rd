\name{moran.I.multi}
\alias{moran.I.multi}
\alias{moran.I.uni}
\alias{moran.I.basic}
\title{ Moran's I for symmetric and asymmetric weighting matrix }
\description{
  Calculate Moran's I for spatial variables out of the eigenfunction-based spatial filtering framework. The function also tests each Moran's I with a permutation or parametric test. The calculation of Moran's I and the test of significance are carried out by the function.

}

\usage{
moran.I.uni(x, mat.W, scaled = FALSE,normalize=FALSE, na.rm = FALSE, test.type="permutation",nperm=999,alternative = "greater")
moran.I.multi(eigenvector.mat,link,weight,scaled=FALSE,normalize=FALSE,na.rm = FALSE,test.type="permutation",nperm=999,plot.res=TRUE)
}

\arguments{
  \item{x}{ A numeric vector. }
  \item{mat.W}{ A matrix of weights. }
  \item{eigenvector.mat}{ Matrix. A set of orthogonal spatial variables, created by the function \code{aem}. }
  \item{link}{ A 2-column matrix describing the link edges. It has 2 columns (from, to) and as many rows as there are edges. The object names in the From-To list are the order numbers of the objects, not their names if the names differ from the order numbers. }
  \item{weight}{ A vector providing weights associated to the edges. If no weights are given, the function consider all edges to have the same weights.}
  \item{scaled}{ Logical (\code{TRUE}, \code{FALSE}) determining if the index should be scaled to allow comparisons between indices (default to \code{FALSE}). }
  \item{normalize}{ Logical (\code{TRUE}, \code{FALSE}) determining if the matrix of weights should be row-normalized to 1 (default to \code{FALSE}). }
  \item{na.rm}{ Logical (\code{TRUE}, \code{FALSE}) determining whether missing values should be removed (default to \code{FALSE}). }
  \item{alternative}{ A character string specifying the alternative hypothesis that is tested against the null hypothesis of no spatial autocorrelation; must be of one "two.sided", "less", or "greater", or any unambiguous abbrevation of these. }
  \item{test.type}{ A character string specifying the type of test to be carried out. Either "permutation", or "parametric", or any unambiguous abbrevation of these. Default is "permutation"}
  \item{nperm}{ Numeric. a number specifying the number of permutation to be carried out. This argument is inactive when performing a parametric test. }
  \item{plot.res}{ Logical (\code{TRUE}, \code{FALSE}) determining if the expected and observed values and the results of the test of significance of Moran's I, for each AEM eigenfunction, should be plotted (default to \code{TRUE}). }
}
\details{
  \code{moran.I.uni} is a modification of \code{\link[ape]{Moran.I}} (the code was heavily borrowed from library \pkg{ape}). In \code{\link[ape]{Moran.I}} the normalization is carried out by default, whereas in \code{moran.I.uni} the choice is given to the user.
  
  Other type of orthogonal spatial variables (MEM, PCNM) created through the Moran's eigenvector maps framework (result of \code{\link[spacemakeR]{scores.listw}} or \code{\link[spacemakeR]{pcnm}}) can also be tested with those functions. 
  
  Function \code{moran.I.basic} is simple function to compute Moran's I. It is used by both \code{moran.I.uni} and \code{moran.I.multi}.
}
\value{
  \item{ observed }{ The computed Moran's I}
  \item{ expected }{ Numeric. The expected Moran's I under the null hypothesis. }
  \item{ sd }{ Numeric. The standard deviation of the Moran's I under the null hypothesis. This value is calculated only during parametric tests.}
  \item{p.value}{ The P-value of the null hypothesis's test against the alternative hypothesis specified in \code{ alternative } }
  \item{ res.mat }{ A 2-column matrix. The first column gives the observed Moran's I value, the second column gives the associated p-value. }
}

\author{ F. Guillaume Blanchet }

\seealso{ \code{\link[ape]{Moran.I}} }

\examples{

### Construct AEM eigenfunctions
nb<-cell2nb(5,5,"queen")
xy <- cbind(1:25,expand.grid(1:5,1:5))
bin.mat <- build.binary(nb,xy)
aem1 <- aem(build.binary=bin.mat)

### Calculate and test Moran's I for each AEM eigenfunction
moran.I.multi(aem1$vectors, bin.mat[[2]])

#---------------------------
### Example using spacemakeR
#---------------------------
require(spacemakeR)

### Construct Moran's eigenvector maps
nb<-cell2nb(5,5,"queen")
sc.queen<-scores.listw(nb2listw(nb,style="B"))

### Calculate and test Moran's I for each MEM eigenfunction
moran.I.multi(sc.queen$vector,listw2sn(nb2listw(nb))[,1:2])

}
\keyword{ spatial }