'PCNM' <- 
function(matdist, thresh=NULL, all=FALSE, include.zero=FALSE)
#
# Compute the PCNM eigenfunctions corresponding to all eigenvalues (+, 0, -).
# Input file: distance matrix produced by the function "dist".
# Computation of the threshold requires a function of the library "ape".
# 
# Original PCNM function: Stephane Dray, November 11, 2004
# The present version: Pierre Legendre, August 2007, January 2009
{
	require(vegan)
	a <- system.time({
	cat("Truncation level =",thresh,'\n')
	matdist <- as.matrix(matdist)

	## Truncation of distance matrix
	if(is.null(thresh)) {
		spanning <- vegan::spantree(as.dist(matdist))
		thresh <- max(spanning$dist)
		}
	matdist[matdist > thresh] <- 4*thresh

	mypcnm.all <- pcoa.all(matdist, all=all, include.zero=include.zero, rn=rownames(matdist))
	})
	a[3] <- sprintf("%2f",a[3])
	cat("Time to compute PCNMs =",a[3]," sec",'\n')
	return(mypcnm.all)
}
