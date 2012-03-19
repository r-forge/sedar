aem.time <- function(n, w=NULL, moran=FALSE,plot.moran=FALSE){
#
# Construct AEM eigenfunctions for a regular time series
# n = number of points
# w = vector of weights. The weights can be the inverse of the interval lengths if the observations are not regularly spaced. Default: equal weights.
# moran: compute Moran's I for each AEM eigenfunction
#
# Authors: Pierre Legendre and F. Guillaume Blanchet, March 2012
	
	# Normalize a vector (to length 1)
	normalize <- function(vec)  vec/sqrt(sum(vec^2)) 
	###  End internal functions
	#
	epsilon <- sqrt(.Machine$double.eps) 
	#
	# Construct matrix E
	E <- matrix(0,n,(n-1))
	rownames(E) <- paste("site",1:n,sep=".")
	colnames(E) <- paste("E",1:(n-1),sep="")
	for(i in 2:n) E[i,1:(i-1)] <- 1
	#
	# Apply weights if provided
	if(!is.null(w)) {
		if(length(w) != (i-1)) stop("Length of vector w not equal to (n-1)")
		E <- E %*% diag(w)
	} else {
		w <- rep(1,n-1)
	}
	#
	# Compute AEM eigenfunctions
	E.c <- scale(E, center=TRUE, scale=FALSE)
	E.svd <- svd(E.c)
	k <- length(which(E.svd$d > epsilon))
	# Normalize the AEM eigenfunctions
	E.svd$u[,1:k] <- apply(E.svd$u[,1:k], 2, normalize)
	
	xy <- 1:n
	if(moran) {
		nb <- cell2nb(n,1)
		fr.to.aem <- rm.double.link(listw2sn(nb2listw(nb))[,1:2])
		res <- moran.I.multi(E.svd$u[,1:k], link=fr.to.aem, weight=w,plot.res=plot.moran)
		Moran <- res$res.mat[,1:2]
		positive <- rep(FALSE,k)
		positive[which(Moran[,1] > res$expected)] <- TRUE
		Moran <- cbind(as.data.frame(Moran), positive)
		colnames(Moran) <- c("Moran","p.value","Positive")
		out <- list(E=E, values=E.svd$d[1:k]^2/(n-1), aem=E.svd$u[,1:k],  
			Moran=Moran, expected_Moran=res$expected)
		} else {
			out <- list(E=E, values=E.svd$d[1:k]^2/(n-1), aem=E.svd$u[,1:k])
		}
	out
}
