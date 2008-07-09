\name{moran.I.multi}
\alias{moran.I.multi}
\alias{moran.I.uni}
\title{ Moran's I for symmetric and asymmetric weighting matrix }
\description{
  Calculates Moran's I for spatial variables out of the eigenfunction-based spatial filtering framework. The function also tests each Moran's I with a parametric test. The calculation of Moran's I and the test of significance are carried out by the function.

}

\usage{
moran.I.uni(x, mat.W, scaled = FALSE, normalize = FALSE, na.rm = FALSE, alternative = "two.sided")
moran.I.multi(eigenvector.mat, link, weight, scaled = FALSE, normalize = FALSE, alternative = "two.sided", plot.res = TRUE)
}

\arguments{
  \item{x}{ A numeric vector. }
  \item{mat.W}{ A matrix of weights. }
  \item{eigenvector.mat}{ Matrix. A set of orthogonal spatial variables, created by the function \code{aem}. }
  \item{link}{ A 2-column matrix describing the link edges. It has 2 columns (from, to) and as many rows as there are edges. }
  \item{weight}{ A vector providing weights associated to the edges. If no weights are given, the function consider all edges to have the same weights.}
  \item{scaled}{ Logical (\code{TRUE}, \code{FALSE}) determining if the index should be scaled to allow comparisons between indices (default to \code{FALSE}). }
  \item{normalize}{ Logical (\code{TRUE}, \code{FALSE}) determining if the matrix of weights should be row-normalized to 1 (default to \code{FALSE}). }
  \item{alternative}{ A character string specifying the alternative hypothesis that is tested against the null hypothesis of no spatial autocorrelation; must be of one "two.sided", "less", or "greater", or any unambiguous abbrevation of these. }
  \item{plot.res}{ Logical (\code{TRUE}, \code{FALSE}) determining if the expected and observed values and the results of the test of significance of Moran's I, for each AEM eigenfunction, should be plotted (default to \code{TRUE}). }
}
\details{
  \code{moran.I.uni} is a modification of \code{\link{Moran.I}}. In \code{\link{Moran.I}} the normalization is carried out by default, whereas in \code{moran.I.uni} the choice is given to the user.
  
  Other type of orthogonal spatial variables (MEM, PCNM) created through the Moran's eigenvector maps framework (result of \code{\link{scores.listw}} or \code{\link{pcnm}}) can also be tested with those functions. 
}
\value{
  \item{ observed }{ The computed Moran's I}
  \item{ expected }{ Numeric. The expected Moran's I under the null hypothesis. }
  \item{ sd }{ Numeric. The standard deviation of the Moran's I under the null hypothesis. }
  \item{p.value}{ The P-value of the null hypothesis's test against the alternative hypothesis specified in \code{ alternative } }
  \item{ res.mat }{ A 2-column matrix. The first column gives the observed Moran's I value, the second column gives the associated p-value. }
}

\author{ F. Guillaume Blanchet }

\seealso{ \code{\link{Moran.I}} }

\examples{

### Construct AEM eigenfunctions
nb<-cell2nb(5,5,"queen")
xy <- cbind(1:25,expand.grid(1:5,1:5))
bin.mat <- build.binary(nb,xy)
aem1 <- aem(build.binary=bin.mat)

### Calculate Moran's I for each AEM eigenfunction
moran.I.multi(aem1$vectors, bin.mat[[2]])

}
\keyword{ spatial, Moran }