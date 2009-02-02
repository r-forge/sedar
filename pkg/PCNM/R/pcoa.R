pcoa <- function(D, correction="none", rn=NULL)
#
# Principal coordinate analysis (PCoA) of a square distance matrix D
# with correction for negative eigenvalues.
# 
# References:
# Gower, J. C. 1966. Some distance properties of latent root and vector methods
#    used in multivariate analysis. Biometrika. 53: 325-338.
# Gower, J. C. and P. Legendre. 1986. Metric and Euclidean properties of 
#    dissimilarity coefficients. J. Classif. 3: 5-48.
# Legendre, P. and L. Legendre. 1998. Numerical ecology, 2nd English edition. 
#    Elsevier Science BV, Amsterdam. [PCoA: Section 9.2]
#
#      Pierre Legendre, October 2007
{
centre <- function(D,n)
# Centre a square matrix D by matrix algebra
# mat.cen = (I - 11'/n) D (I - 11'/n)
{	One <- matrix(1,n,n)
	mat <- diag(n) - One/n
	mat.cen <- mat %*% D %*% mat
}

# ===== The PCoA function begins here =====

# Preliminary actions
	D <- as.matrix(D)
	n <- nrow(D)
	epsilon <- sqrt(.Machine$double.eps)
	if(length(rn)!=0) {
		names <- rn
		} else {
		names <- rownames(D)
		}
	CORRECTIONS <- c("none","lingoes","cailliez")
	correction <- pmatch(correction, CORRECTIONS)
	# cat("correction =",correction,'\n')

# Gower centring of matrix D
# delta1 = (I - 11'/n) [-0.5 d^2] (I - 11'/n)
	delta1 <- centre((-0.5*D^2),n)
	trace <- sum(diag(delta1))

# Eigenvalue decomposition
	D.eig <- eigen(delta1)
	
# Negative eigenvalues?
	min.eig <- min(D.eig$values)
	zero.eig <- which((D.eig$values > -epsilon) & (D.eig$values < epsilon))
	D.eig$values[zero.eig] <- 0

# No negative eigenvalue
	if(min.eig > -epsilon) { 
		correction <- 3
		eig <- D.eig$values
		k <- length(which(eig > epsilon))
		rel.eig <- eig[1:k]/trace
		cum.eig <- cumsum(rel.eig) 
		vectors <- sweep(D.eig$vectors[,1:k], 2, sqrt(eig[1:k]), FUN="*")
		bs <- broken.stick(k)[,2]
		cum.bs <- cumsum(bs)

res <- data.frame(eig[1:k], rel.eig, bs, cum.eig, cum.bs)
colnames(res) <- c("Eigenvalues","Relative_eig","Broken_stick","Cumul_eig","Cumul_br_stick")

rownames(vectors) <- names
colnames(vectors) <- colnames(vectors, do.NULL = FALSE, prefix = "Axis.")
out <- (list(values=res, vectors=vectors, trace=trace))

# Negative eigenvalues present
	} else {
		k <- n 
		eig <- D.eig$values
		rel.eig <- cumsum(eig)
		rel.eig.cor <- (eig - min.eig)/(trace - (n-1)*min.eig) # Eq. 9.27 for a single dimension
		rel.eig.cor = c(rel.eig.cor[1:(zero.eig[1]-1)], rel.eig.cor[(zero.eig[1]+1):n], 0)
		cum.eig.cor <- cumsum(rel.eig.cor) 
		k2 <- length(which(eig > epsilon))
		k3 <- length(which(rel.eig.cor > epsilon))
		vectors <- sweep(D.eig$vectors[,1:k2], 2, sqrt(eig[1:k2]), FUN="*")
		# Only the eigenvectors with positive eigenvalues are shown

# Negative eigenvalues: three ways of handling the situation
	if((correction==1) | (correction==2)) {

		# Lingoes correction: compute c1, then the corrected D
		if(correction == 1) {
			c1 <- -min.eig
			note <- paste("Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 -",c1,", except diagonal elements")
			D <- -0.5*(D^2 + 2*c1)

		# Cailliez correction: compute c2, then the corrected D
		} else if(correction == 2) {
			delta2 <- centre((-0.5*D),n)
			upper <- cbind(matrix(0,n,n), 2*delta1)
			lower <- cbind(-diag(n), -4*delta2)
			sp.matrix <- rbind(upper, lower)
			c2 <- max(Re(eigen(sp.matrix, symmetric=FALSE, only.values=TRUE)$values))
			note <- paste("Cailliez correction applied to negative eigenvalues: D' = -0.5*(D +",c2,")^2, except diagonal elements")
			D <- -0.5*(D + c2)^2
		}

	diag(D) <- 0
	mat.cor <- centre(D,n)
	toto.cor <- eigen(mat.cor)
	trace.cor <- sum(diag(mat.cor))

	# Negative eigenvalues present?
	min.eig.cor <- min(toto.cor$values)
	zero.eig.cor <- which((toto.cor$values < epsilon) & (toto.cor$values > -epsilon))
	toto.cor$values[zero.eig.cor] <- 0

	# No negative eigenvalue after correction: result OK
	if(min.eig.cor > -epsilon) { 
		eig.cor <- toto.cor$values
		rel.eig.cor <- eig.cor[1:k]/trace.cor
		cum.eig.cor <- cumsum(rel.eig.cor) 
		k2 <- length(which(eig.cor > epsilon))
		vectors.cor <- sweep(toto.cor$vectors[,1:k2], 2, sqrt(eig.cor[1:k2]), FUN="*")
		bs <- broken.stick(k2)[,2]
		bs <- c(bs, rep(0,(k-k2)))
		cum.bs <- cumsum(bs)

	# Negative eigenvalues still present after correction: incorrect result
	} else { 
	if(correction == 1) cat("Probleme! Negative eigenvalues are still present after Lingoes",'\n') 
	if(correction == 2) cat("Probleme! Negative eigenvalues are still present after Cailliez",'\n') 
	rel.eig.cor <- cum.eig.cor <- bs <- cum.bs <- rep(NA,n)
	vectors.cor <- matrix(NA,n,2)
	}

	res <- data.frame(eig[1:k], eig.cor[1:k], rel.eig.cor, bs, cum.eig.cor, cum.bs)
	colnames(res) <- c("Eigenvalues", "Corr_eig", "Rel_corr_eig", "Broken_stick", "Cum_corr_eig", "Cum_br_stick")

	rownames(vectors) <- names
	colnames(vectors) <- colnames(vectors, do.NULL = FALSE, prefix = "Axis.")
	out <- (list(correction=note, values=res, vectors=vectors, trace=trace, vectors.cor=vectors.cor, trace.cor=trace.cor))

	} else {
			
	note <- "No correction was applied to the negative eigenvalues"
	bs <- broken.stick(k3)[,2]
	bs <- c(bs, rep(0,(k-k3)))
	cum.bs <- cumsum(bs)

	res <- data.frame(eig[1:k], rel.eig, rel.eig.cor, bs, cum.eig.cor, cum.bs)
	colnames(res) <- c("Eigenvalues","Relative_eig","Rel_corr_eig","Broken_stick","Cum_corr_eig","Cumul_br_stick")

	rownames(vectors) <- names
	colnames(vectors) <- colnames(vectors, do.NULL = FALSE, prefix = "Axis.")
	out <- (list(correction=note, values=res, vectors=vectors, trace=trace))

			}
		}	# End negative engenvalues
	class(out) <- "pcoa"
	out
}	# End of PCoA
