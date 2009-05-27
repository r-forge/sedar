\name{STImodels}
\alias{STImodels}
\alias{quickSTI}
\alias{manovRDa}
\alias{restrictedPerm}
\title{ Space-time interaction in ANOVA without replicates}
\description{
Function \code{\link{STImodels}} performs two-way ANOVA to test space-time interaction without replicates using one among a set of possible models. Function \code{\link{quickSTI}} allows performing space-time ANOVA in a simplified way. In many models, degrees of freedom are saved by coding space and/or time parsimoniously using principal coordinates of neighbour matrices (\code{\link{PCNM}}). 
}
\usage{
STImodels(Y, S, T,  model="5", nperm=999, nS=-1,nT=-1, Sfixed=TRUE, Tfixed=TRUE, COD.S=NULL, COD.T=NULL, print.res=TRUE)
quickSTI(Y, S, T, nperm=999, alpha = 0.05, COD.S=NULL, COD.T=NULL,print.res=TRUE)
}

\arguments{
  \item{Y}{ Site-by-species response data table (assumes row blocks corresponding to times, i.e. within each block all sites occur)}
  \item{S}{ The number of spatial points or a matrix of spatial coordinates. }
  \item{T}{ The number of time campaigns or a matrix (a vector) of temporal coordinates. }
  \item{model}{ Linear space-time model to be used (can be either "2", "3a", "3b", "4","5","6a", "6b", or "7"). }
  \item{nperm}{ The number of permutations to be done. }
  \item{alpha}{ Confidence level for the interaction test. Depending on the decision for the interaction test, the main factors are tested differently. }
  \item{nS}{ Number of space PCNMs to use (by default -1 means \code{s/2} rounding down). }
  \item{nT}{ Number of time PCNMs to use (by default -1 means \code{t/2} rounding down). }
  \item{Sfixed}{ A logical value: is Factor Space fixed, or not (if FALSE, it is considered a random factor). }
  \item{Tfixed}{ A logical value: is Factor Time fixed, or not (if FALSE, it is considered a random factor). }
  \item{COD.S}{ The spatial coding functions to be used instead of PCNM. Their number must be smaller than \code{S}. }
  \item{COD.T}{ The temporal coding functions to be used instead of PCNM. Their number must be smaller than \code{T}. }
  \item{print.res}{ If TRUE prints additional information. }
}

\details{
In \code{STImodels} tests for space-time interaction and space or time main effects are conducted using one of the different models. Depending on the model choosen, the interaction test is not available (concretely, with Models 2 and 6). 

Model 2  - Space and Time are coded using Helmert contrasts for the main effects. No interaction is tested.
Model 3a - Space is coded using PCNM variables whereas Time is coded using Helmert contrasts.
Model 3b - Space is coded using Helmert contrasts whereas Time is coded using PCNM variables.
Model 4  - Both Space and Time are coded using PCNM variables for all tests.
Model 5  - Space and Time are coded using Helmert contrasts for the main factor effects, but they are coded using PCNM variables for the interaction term.
Model 6a - Nested model. Testing for the existence of temporal structure (common or separated) using PCNM variables to code for Time.
Model 6b - Nested model. Testing for the existence of spatial structure (common or separated) using PCNM variables to code for Space.
Model 7  - Space and Time are coded using PCNM variables for the main factor effects, but they are coded using Helmert contrasts for the interaction term (not recommended).

When using \code{quickSTI}, space-time interaction is first tested using Model 5. Depending on the outcome of this test, the main factors are tested using different strategies. If the interaction is not significant then the test of main factors is also done following Model 5. If the interaction is significant, then a nested model is used to know whether separate spatial structures exist and another to know whether separate temporal structures exist. In \code{quickSTI} function space and time are always considered fixed factors (F ratios are constructed using residual MS in the denominator).

For the interaction the permutations are unrestricted, whereas for the main factors the permutations are restricted within time or space blocks. By default, the function computes \code{\link{PCNM}} for space and time coding, but the coding can be also provided by the user, through \code{COD.S} and \code{COD.T}. 
}

\value{

  \item{testS}{ An object with the result of the space effect test, including the mean squares for the F numerator (\code{MS.num}), the mean squares for the F denominator (\code{MS.den}), the proportion of explained variance (\code{R2}), the adjusted proportion of explained variance (\code{R2.adj}), the F statistics (\code{F}) and its p-value computed from permutation test (\code{Prob}). }
  \item{testT}{ An object with the result of the time effect test, like \code{testS}.}
  \item{testSTI}{ An object with the result of the space-time interaction effect test, like \code{testS}.}
}

\author{ Pierre Legendre and Miquel De Caceres, Universite de Montreal }

\references{ 
Borcard, D. and P. Legendre. 2002. All-scale spatial analysis of ecological data by means of principal coordinates of neighbour matrices. Ecological Modelling 153: 51-68.

Dray, S., P. Legendre and P. R. Peres-Neto. 2006. Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). Ecological Modelling 196: 483-493.



Legendre, P., M. De Caceres and  D. Borcard 2009. Community surveys through space and time to assess environmental changes: testing the space-time interaction in the absence of replication. Ecology (in press).
}

\seealso{\code{\link{PCNM}}, \code{\link{trichoptera}} }

\examples{
data(trichoptera)

# Log-transform the species data (excluding site and time colums)
trich.log = log(trichoptera[,-c(1,2)]+1)

# Run space-time interaction test using model "5"
STImodels(trich.log, S=22, T=10, nperm=999, model="5")

# Run space-time analysis with tests for main effects after testing interaction 
# (which is significant)
quickSTI(trich.log, S=22, T=10, nperm=999)

# Run space-time analysis for time blocks 6 and 7.
# Interaction is then not significant 
# and the tests of the main effects are done using model 5
quickSTI(trich.log[111:154,], S=22, T=2, nperm=999)
}

\keyword{ multivariate }
\keyword{ spatial }
