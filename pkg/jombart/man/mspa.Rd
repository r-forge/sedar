\name{mspa}
\alias{mspa}
\alias{scatter.mspa}
\alias{print.mspa}
\title{Multi-Scale Pattern Analysis}
\description{
  The multi-scale pattern analysis (MSPA, Jombart et al 2009) investigates the main scales
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
  \item{\dots}{other parameters passed to \code{\link[pkg:ade4]{s.arrow}}.}
}
\value{
  An object having the classes \code{mspa} and \code{\link[pkg:ade4]{dudi}}: \code{mspa}
  objects are \code{\link[pkg:ade4]{dudi}} objects with the following
  extra slots:\cr
  - ls: principal components of the MSPA. These are the coordinates of
             variables onto principal axes, to be used for
             plotting. Correspond to matrix \bold{B} in Appendix A of
             Jombart et al (2009). \cr
  - R2: matrix of R2 between variables and MEMs. Corresponds to \bold{S}
             in Jombart et al (2009).\cr
  - meanPoint: coordinates of the 'mean variable' onto principal
             axes. The 'mean variable' is an hypothetic variable whose
             scale profile is the average of those of all variables of
             the analysed dataset.\cr
  - varweights: the weights of variables. Corresponds to \bold{d} in
             Jombart et al. (2009).\cr
}
\references{
  Jombart T, Dray S, and Dufour, A-B. (2009) Finding essential scales of
  spatial variation in ecological data: a multivariate
  approach. \emph{Ecography} \bold{32}: 161-168.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\seealso{\code{\link{chooseCN}} to obtain a list of spatial weights.}
\examples{
############################
### using simulated data ###
############################

## load data
data(mspaIllus)
X <- mspaIllus
dim(X)
lw <- attr(mspaIllus,"lw")

## perform a PCA
pcaX <- dudi.pca(X,scannf=FALSE,nf=4)

## perform the MSPA
mspa1 <- mspa(pcaX,lw,scannf=FALSE,nf=3)

## plot the results
scatter(mspa1)
scatter(mspa1,xax=1,yax=3)



####################################
### using oribatib mites dataset ###
####################################

## load data
data(oribatid)

## get the list of spatial weights
cn <- chooseCN(oribatid$xy,res="listw",ask=FALSE,type=1)

## Hellinger transformation
hellTrans <- function(X){
  if (!( is.matrix(X) | is.data.frame(X) )) stop("Object is not a matrix.")  
  if (any(is.na(X))) stop("na entries in table.")
  
  sumRow <- apply(X,1,sum)
  Y <- X/sumRow
  Y <- sqrt(Y)
  
  return(Y)
}


## ENVIRONMENTAL VARIABLES ##
## Hill and Smith analysis for environmental variables
## (for a mixture of quantitative / qualitative variables)
hsEnv <- dudi.hillsmith(oribatid$envir,scannf=FALSE)

## detrending of the analysis (residuals of regression onto xy coordinates)
hsEnv.detr <- pcaivortho(hsEnv,oribatid$xy,scannf=FALSE)

## MSPA of the detrended analysis
mspaEnv <- mspa(hsEnv.detr,cn,scannf=FALSE,nf=2)
scatter(mspaEnv)



## SPECIES DATA ##
## PCA of species abundances, after Hellinger transformation
pcaFau <- dudi.pca(hellTrans(oribatid$fau),scale=FALSE,scannf=FALSE)

## detrending of this PCA
pcaFau.detr <- pcaivortho(pcaFau,oribatid$xy,scannf=FALSE)

# MSPA of the detrended analysis
mspaFau <- mspa(pcaFau.detr,cn,scannf=FALSE,nf=2)
scatter(mspaFau)



## CANONICAL MSPA ##
## RDA species ~ envir
## (species abundances predicted by environment)
## note: RDA = 'PCAIV' (PCA with Instrumental Variables)
rda1 <- pcaiv(dudi=pcaFau.detr, df=oribatid$envir,scannf=FALSE,nf=2)

## canonical MSPA (species predicted by environment)
mspaCan1 <- mspa(dudi=rda1, lw=cn, scannf=FALSE, nf=2)
scatter(mspaCan1)

## same analysis, using a non-parametric centring
mspaCan1NP <- mspa(dudi=rda1, lw=cn, scannf=FALSE, nf=2,cent="sim",nper=999)
scatter(mspaCan1NP) # basically no change



## PARTIAL CANONICAL MSPA ##
## partial RDA species ~ envir
## (species abundances not predicted by environment)
rda2 <- pcaivortho(dudi=pcaFau.detr,df=oribatid$envir,scannf=FALSE,nf=2)

## partial canonical MSPA
mspaCan2 <- mspa(dudi=rda2, lw=cn, scannf=FALSE, nf=2)
scatter(mspaCan2) # nothing left

}
\keyword{multivariate}
\keyword{spatial}
