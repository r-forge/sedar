\name{quickPCNM}
\alias{quickPCNM}

\title{ Quick PCNM spatial eigenfunction analysis }
\description{
Function \code{\link{quickPCNM}} computes principal coordinates of neighbour matrices (PCNM) and carries out the spatial eigenfunction analysis. }
\usage{
quickPCNM(Y, space, thresh=0, method="fwd", myPCNM=NULL, alpha=0.05, rangexy=FALSE, detrend=TRUE, perm.max=NULL)
}

\arguments{
  \item{Y}{ A rectangular data table of class \code{data.frame} or \code{matrix} containing response data, ready for the analysis, and pretransformed if necessary. }
  \item{space}{ A rectangular data table of class \code{data.frame} or \code{matrix} containing the geographic (x or x-y) coordinates of the sites. }
  \item{thresh}{ A user-provided threshold for truncation of the geographic distance matrix (optional). }

  \item{method}{ specifies the method of selection of PCNM variables: }
  \item{}{ "none": no selection. All spatial eigenfunctions are used. May be useful when the user provides precomputed and preselected spatial eigenfunctions. }
  \item{}{ "fwd": regression-based forward selection. The fastest method, default, double stopping criterion [alpha level and adjusted R2]. This is the strongly recommended choice. }

  \item{myPCNM}{ An optional matrix containing precomputed spatial eigenfunctions like PCNM, MEM or AEM, or a polynomial of the geographic coordinates. It may also contain almost all PCNMs with positive and negative eigenvalues computed by the function \code{\link{PCNM}} with options 'all=TRUE, exclude.zero=TRUE'. The user will have to keep the number of explanatory variables below n-1, however.}
  \item{alpha}{ Level of significance for the tests (default value: 0.05). }
  \item{rangexy}{ If TRUE, the spatial coordinates are rescaled to a minimum of 0 on both axes and a maximum of 10 on the axis wih the largest range, without altering the ratio between the x and y coordinates. If only one coordinate is provided, it is rescaled into a range [0,10]. Useful when the x-y coordinates are provided in extremely large or small values. Default: FALSE. }
  \item{detrend}{ If TRUE (default), prior ot the PCNM analysis, the function tests if a linear trend is present and significant at the level given by 'alpha'. If that is the case, a linear detrending of the response data is done. One should not waste sine-shaped PCNM variables to model a linear trend. }
  \item{perm.max}{ Maximum number of permutations in tests of siignificance, e.g. perm.max=1000. If perm.max=NULL, vegan's default method will be used. }
}

\details{
If a truncation value is not provided, the largest distance in a minimum spanning tree linking all objects on the map is computed. That value is used as the truncation value.

BEWARE - In the forward selection of explanatory eigenfunctions, if the rda test involving all the PCNM eigenfunctions is nonsignificant at the selected alpha level, the procedure is stopped, because of the high risk of type I error associated with forward selection of explanatory variables (Blanchet et al. 2008). Users of the PCNM method MUST run an overall test prior to ANY type of variable selection. Borcard and Legendre (2002) have shown that the overall test involving all PCNMs has a correct level of type I error. Furthermore, a selection of variables should never produce a model with an adjusted R-square higher than the value obtained for the global model containing all the variables (Blanchet et al. 2008). If that happens, a warning is issued. If the difference is important, the user should try a more conservative selection method.

The function provides spatial plots of the several first canonical axes (our thanks to Guillaume Blanchet for improving these plots) and a list of results containing the PCNM variables, RDA results and a global test of the RDA.

Several summary diagnostics are displayed on-screen at the end of the computation.

If you want to run the function again to compare results, don't forget to ask for a new graphics window. Otherwise the active window created during the first run will be overwritten.

When the run is completed, type \code{summary(name_of_object)} to find the RDA results in the output object.
}

\value{

  \item{PCNM }{The PCNM eigenfunctions (principal coordinates of the truncated distance matrix). They are either the classical PCNMs with positive eigenvalues computed by this function, or the table of PCNMs that have been entered by the user. }
  \item{PCNM_eigenvalues }{The PCNM eigenvalues. }
  \item{fwd.sel }{The forward selection results produced by function \code{forward.sel} of the \code{packfor} library. }
  \item{PCNM_reduced_model }{A table containing only the PCNMs selected by forwartd selection. }
  \item{RDA }{The RDA results computed by function \code{\link[vegan]{rda}} of the \pkg{vegan} package. }
  \item{RDA_test }{The anova table reporting the permutation test of significance of the F-statistic of the RDA, computed by function \code{\link[vegan]{anova.cca}} of the \pkg{vegan} package. }
  \item{RDA_axes_tests }{The anova table reporting the permutation tests of significance of the F-statistic of the RDA canonical axes, computed by function \code{\link[vegan]{anova.cca}} of the \pkg{vegan} package. }
}

\references{ 
Blanchet F. G., P. Legendre, and D. Borcard. 2008. Forward selection of explanatory variables. Ecology 89: 2623-2632.

Borcard, D. and P. Legendre. 1994. Environmental control and spatial structure in ecological communities: an example using Oribatid mites (Acari, Oribatei). Environmental and Ecological Statistics 1: 37-53.
     
Borcard, D. and P. Legendre. 2002. All-scale spatial analysis of ecological data by means of principal coordinates of neighbour matrices. Ecological Modelling 153: 51-68.

Borcard, D., P. Legendre, C. Avois-Jacquet and H. Tuomisto. 2004. Dissecting the spatial structure of ecological data at multiple scales. Ecology 85: 1826-1832 [Statistical report].

Dray, S., P. Legendre and P. R. Peres-Neto. 2006. Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). Ecological Modelling 196: 483-493.
}

\author{ Daniel Borcard, Universite de Montreal }

\seealso{\code{\link{pcoa.all}}, \code{\link{PCNM}} }

\examples{
# Oribatid mite data from Borcard and Legendre (1994)
library(vegan)
data(mite)      # 70 peat cores, 35 species
data(mite.xy)   # Core geographic coordinates

mite.hel <- decostand(mite, "hellinger")

# A. Usual case with defaut settings; works well in most situations
quick.1 <- quickPCNM(mite.hel, mite.xy)
summary(quick.1)

# B. Same with user-provided truncation level
quick.2 <- quickPCNM(mite.hel, mite.xy, thresh=1.500)

# C. Same as A but with precomputed PCNM variables, a modified  
# alpha level and rescaling of x-y coordinates
quick.3 <- quickPCNM(mite.hel, mite.xy, myPCNM=quick.1$PCNM, alpha=0.01, range=TRUE)
}

\keyword{ multivariate }
\keyword{ spatial }
