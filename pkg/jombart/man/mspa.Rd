\name{mspa}
\alias{mspa}
\alias{scatter.mspa}
\title{Multi-Scale Pattern Analysis}
\description{
  The multi-scale pattern analysis (MSPA) investigates the main scales
  of spatial variation in a multivariate dataset. This implementation
  allows one to perform a MSPA using any multivariate analysis (stored
  as a \code{\link[pkg:ade4]{dudi}} object), and a list of spatial
  weights.\cr

  The \code{scatter} method is used for plotting the results. 
}
\usage{
mspa(dudi, lw, scannf = TRUE, nf = 2, centring=c("param","sim"),
     nperm=1000)

\method{scatter}{mspa}(x, xax = 1, yax = 2, clab.var = 0.75, clab.sca = 1, 
             posieig = "top", sub = NULL, ratio = 1/4, bary=TRUE, circle=TRUE, \dots)
}
\arguments{
  \item{dudi}{a duality diagram (i.e. a reduced space ordination)
    obtained by a \code{\link[pkg:ade4]{dudi}} function (for instance
    \code{\link[pkg:ade4]{dudi.pca}}).}
  \item{lw}{a list of weights of class \code{listw}. Can be obtained
    easily using the function \code{\link{chooseCN}}.}
  \item{scannf}{logical, indicating whether the screeplot should be
    displayed to choose the number or retained factors.}
  \item{nf}{the number of retained factors}
  \item{centring}{a character string indicating if parametric ("param")
    or non-parametric ("sim") centring should be used}
  \item{nperm}{an integer giving the number of permutations used to
    compute the theoretical coefficients of determination (1000 by
    default); used if centring="sim".}
  \item{x}{a mspa object.}
  \item{xax}{an integer indicating the x axis to be displayed.}
  \item{yax}{an integer indicating the y axis to be displayed.}
  \item{clab.var}{a numeric indicating the size of the labels of variables.}
  \item{clab.sca}{a numeric indicating the size of the labels of scales.}
  \item{posieig}{a character indicating the position of the screeplot
    (any of the four combination between "top", "bottom", "left" and
    "right").}
  \item{sub}{an optional character string added as legend to the plot.}
  \item{ratio}{a number between 0 and 1 indicating the size of the
    screeplot as a proportion of the plot.}
  \item{bary}{a logical indicating whether the barycenter of the variables should be displayed.}
  \item{circle}{a logical indicating whether a circle of radius one should be displayed.}
  \item{\dots}{}
}
\details{
 These tests were originally part of the adegenet package for R.
}
\value{
An object of class \code{randtest}.
}
\references{
  Jombart T, Dray S, and Dufour, A-B. (2009) Finding essential scales of
  spatial variation in ecological data: a multivariate
  approach. \emph{Ecography} \blod{32}: 161-168.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\seealso{\code{\link{chooseCN}} to obtain a list of spatial weights.}
\examples{

}
\keyword{multivariate}
\keyword{spatial}
