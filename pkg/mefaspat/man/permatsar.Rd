\encoding{UTF-8}
\name{permatsar}
\alias{permatsar}
\title{Matrix Permutation Algorithms or Presence-Absence and Count Data}
\description{
Individual (for count data) or incidence (for presence-absence data) based null models can be generated for community level simulations. This function uses Cristian Hennig's and Bernhad Hausdorf's \strong{prabclus} package to obtain parametric bootstrap matrices based on the estimation of disjunction probability and the SAR model.
}
\usage{
permatsar(m, nb, pd = NULL, mtype = "count", 
times = 100, round.sar = TRUE)
}
\arguments{
  \item{m}{a community data matrix with plots (samples) as rows and species (taxa) as columns.}
  \item{nb}{neighbourhood list.}
  \item{pd}{disjunction probability, if \code{NULL} \code{autoconst} function in package 'prablcus' is used to estimate it.}
  \item{mtype}{matrix data type, either \code{"count"} for count data, or \code{"prab"} for presence-absence data.}
  \item{times}{number of permuted matrices.}
  \item{round.sar}{logical, if \code{TRUE} the sum of the resulting matrices will be the same as that of the original. If \code{FALSE} the original estimated results are used, which may not be the same, because SAR model is fitted according to the log abundances.}
}
\value{
An object of class 'permat'.
  \item{call}{the function call.}
  \item{orig}{the original data matrix used for permutations.}
  \item{perm}{a list of permuted matrices with length \code{times}.}
  \item{specs}{a list of other specifications: \code{nb}, \code{pd}, SAR estimates.}
}
\references{
Hennig, C. and Hausdorf, B. 2004. Distance-based parametric bootstrap tests for clustering of
species ranges. \emph{Computational Statistics and Data Analysis} 45, 875--896. \url{http://stat.ethz.
ch/Research-Reports/110.html}

Hausdorf, B. and Hennig, C. 2003. Biotic Element Analysis in Biogeography. \emph{Systematic Biology}
52, 717--723.

Hausdorf, B. and Hennig, C. 2003. Nestedness of north-west European land snail ranges as a
consequence of differential immigration from Pleistocene glacial refuges. \emph{Oecologia} 135, 102--109.

Hausdorf, B. and Hennig, C. 2007. Null model tests of clustering of species, negative co-occurrence
patterns and nestedness in meta-communities. \emph{Oikos} 116, 818--828.
}
\author{
\enc{Péter Sólymos}{Peter Solymos}, \email{solymos@ualberta.ca}
}
\seealso{
autoconst, prab.sarestimate, randpop.nb, regpop.sar
}
\examples{
m <- matrix(c(
   1,3,2,0,3,1,
   0,2,1,0,2,1,
   0,0,1,2,0,3,
   0,0,0,1,4,3
   ), 4, 6, byrow=TRUE)
m
nb <- nbmc(1:4,c(1,1,2,2),rep(1,4),"dl")
nb
x1 <- permatsar(m,nb,mtype="prab")
x1
x2 <- permatsar(m,nb,mtype="count")
x2
}
\keyword{multivariate}
\keyword{ spatial }
