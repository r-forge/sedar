\name{betaols}
\alias{betaols}
\alias{print.betaols}
\alias{plot.betaols}
\alias{summary.betaols}
\title{
Permutation test for regression on distance matrices between two species groups
}
\description{
Permutation test for differences in slope and intercepts of OLS regression on distance matrices between two species groups.
}
\usage{
betaols(datamatrix, latlong, gr1, gr2 = NULL, x.loc = NULL, method = "sim", gcd = TRUE, times = 100)
\method{print}{betaols}(x, decimals = 4, \dots)
\method{plot}{betaols}(x, leg = TRUE, leg.pos = "topright", leg.text = c("Group 1", "Group 2"), crit=0.05, \dots)
\method{summary}{betaols}(object, decimals = 4, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamatrix}{a matrix or data frame containing the data.}
  \item{latlong}{a matrix or data frame with latitude/longitude (when \code{gcd=TRUE}) or x/y coordinates (when \code{gcd=FALSE}).}
  \item{gr1}{a vector for group membership (columns of \code{datamatrix}) of species belonging to group 1.}
  \item{gr2}{a vector for group membership (columns of \code{datamatrix}) of species belonging to group 2 (optional), \code{NULL} means all species that are not in group 1.}
  \item{x.loc}{locations along the geographic distance axis where intercepts are compared.}
  \item{method}{method for dissimilarity computation, one \code{index} of \code{\link[vegan]{betadiver}} (only character).}
  \item{gcd}{logical, if \code{TRUE}, great circle distances (in kms) are calculated, if \code{FALSE}, Euclidean distances are calculated instead.}
  \item{times}{number of permutations.}
  \item{x, object}{an object of class 'betaols'.}
  \item{decimals}{number of decimal places to passed to print method.}
  \item{leg}{logical, whether legend should be plotted.}
  \item{leg.pos}{position of legend, passed to \code{\link{legend}}.}
  \item{leg.text}{character vector with legend text, passed to \code{\link{legend}}.}
  \item{crit}{critical value for color plotting of locations for comparisons. If \code{NULL}, locations are not plotted.}
  \item{\dots}{other arguments.}
}
\details{
Using linear regression for distance matrices is not necessarily the best choice, but it is a convenient way when the question is to compare slopes and intercepts for parts of the data set (here two subsets defined by two non-overlaping ranges of columns). The permutation of columns among the two groups essentially mean, that one do a regression on a subsample of the original data set, thus expected slope and intercept should be equal to that of the other subset, and the whole data set as well (slightly modified by different sampling biases introduced by the subsetting). Thus simple comparison of confidence levels would yield the same results (under all assumptions met), but this implementation uses permutation test according to works of La Sorte et al (2006, 2008), that is more robust and non-parametric solution for distance matrices.

Summary method returns original coeffeicients for the two groups of species, and differences between slopes and intercept(s) for the two groups and corresponding p-values.

Plot shows scatter of beta index vs. geographic distances for the two groups of species, corresponding regression lines and locations for comparison  (in black if p-value < \code{crit}, else in grey).
}
\value{
  An object of class 'betaols' that is a list with elements:
  \item{call}{the function call.}
  \item{m}{data in \code{x}.}
  \item{gr1}{a vector for group membership in group 1.}
  \item{gr2}{a vector for group membership in group 2.}
  \item{x.loc}{locations along the geographic distance axis where intercepts are compared. If \code{x.loc=NULL}, only 0 is used for comparison. Only unique values are used.}
  \item{n.perm}{number of permutations.}
  \item{mat}{a list containing distance matrices \code{gcd1}, \code{gcd2}, \code{beta1} and \code{beta2} (geographic distances and dissimilarities for the two groups) for plotting and reuse.}
  \item{d.perm}{data frame with \code{times} rows and colums for slope and intercepts, this contains permutation results for further validation.}
  \item{betadiv}{dissimilarity index used by \code{method}.}
  \item{geogr}{geographic distance calculation according to \code{gcd}.}
  \item{coef.gr1}{coefficients for \code{lm} object fitted to group 1.}
  \item{coef.gr2}{coefficients for \code{lm} object fitted to group 2.}
  \item{test}{a matrix with test results (observed differences and p values).}
}
\references{
La Sorte, F.A., McKinney, M.L., Pysek, P., Klotz, S., Rapson, G.L. Celesti-Grapow, L., and Thompson, K. (2008) Distance decay in similarity among European urban floras: the impact of anthropogenic activities on beta diversity. Global Ecology and Biogeography, 17, 363-371.

La Sorte, F.A., and McKinney, M.L. (2006) Compositional similarity and the distribution of geographical range size for native and non-native species in urban floras. Diversity and Distributions, 12, 679-686.
}
\author{Peter Solymos}
\examples{
xy <- cbind(x <- runif(50), y <- runif(50))
gr1 <- c(1:10)
gr2 <- c(11:30)
m <- matrix(c(rpois(500,5), rpois(2000,1)), 50, 50)

(bols1 <- betaols(m, xy, gr1, gr2, gcd=FALSE))
(bols2 <- betaols(m, xy, gr1, x.loc = seq(0.1,1.1,0.2)))
plot(bols2)
hist(bols2$d.perm[,1])
}
\keyword{classes}
\keyword{methods}
\keyword{htest}
\keyword{regression}
\keyword{spatial}
