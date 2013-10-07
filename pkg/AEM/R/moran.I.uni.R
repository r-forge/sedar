moran.I.uni<-function (x, mat.W, scaled = FALSE,normalize=FALSE, na.rm = FALSE, test.type="permutation",nperm=999,alternative = "greater") {
	
	#CC# Match arguments
	test.type<-match.arg(test.type,c("permutation","parametric"))
	
	n <- length(x)
	#CC# Test if mat.W is of a good format
	dim.mat.W<-dim(mat.W)
	if (dim.mat.W[1] != dim.mat.W[2]) 
		stop("'mat.W' must be a square matrix")
	if (dim(mat.W)[1] != n) 
		stop("'mat.W' must have as many rows as observations in 'x'")
	
	#CC# Deal with NA's
	nas <- is.na(x)
	if (any(nas)) {
		if (na.rm) {
			x <- x[!nas]
			n <- length(x)
			mat.W <- mat.W[!nas, !nas]
		}else{
			warning("'x' has missing values: maybe you wanted to set na.rm=TRUE?")
			return(list(observed = NA, expected = ei, sd = NA, 
				p.value = NA))
		}
	}
	
	#CC# Row normalization (Gittleman & Kot, eq. 3 and paragraph before)
	if(normalize){
		ROWSUM <- rowSums(mat.W)
		ROWSUM[ROWSUM == 0] <- 1
		mat.W <- mat.W/ROWSUM
	}
	
	if(test.type=="parametric"){
		moran.res<-moran.I.basic(x,mat.W,scaled=scaled)
		#=============================================
		#### Calculate Standard deviation of Moran's I
		#=============================================
		#CC# (Gittleman & Kot, eq. 7)
		S1 <- 0.5 * sum((mat.W + t(mat.W))^2)
	
		#CC# (Gittleman & Kot, eq. 8)
		S2 <- sum((apply(mat.W, 1, sum) + apply(mat.W, 2, sum))^2)
	
		#CC# (Gittleman & Kot, eq. 6)
		s.sq <- moran.res$s^2
		k <- (sum(moran.res$y^4)/n)/(moran.res$v/n)^2
		sdi <- sqrt((n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * s.sq) - 
			k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * s.sq))/((n - 
			1) * (n - 2) * (n - 3) * s.sq) - 1/((n - 1)^2))
	}
	#CC# Mean of Moran's I (Gittleman & Kot, eq. 5)
	ei <- -1/(n - 1)
	
	#____________________
	### Test of Moran's I
	#____________________
	
	alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
	
	if(test.type=="parametric"){
		#==============
		#### Parametric
		#==============
		pv <- pnorm(moran.res$obs, mean = ei, sd = sdi)
		if (alternative == "two.sided"){
			pv <- if (moran.res$obs <= ei){ 
				2 * pv
			}else{
				2 * (1 - pv)
			}
		}
		if (alternative == "greater") pv <- 1 - pv
	}
	if(test.type=="permutation"){
		#================
		#### Permutation
		#================
		moran.I.boot<-function(varia,i,...){
			return(moran.I.basic(varia[i],...)$obs)
		}
		boot.res <- boot(x, statistic = moran.I.boot, R = nperm, sim = "permutation", mat.W = mat.W,scaled=scaled)
		
		if(alternative=="less"){
			pval<-(length(which(boot.res$t<=boot.res$t0))+1)/(nperm+1)
		}
		if(alternative=="greater"){
			pval<-(length(which(boot.res$t>=boot.res$t0))+1)/(nperm+1)
		}
		if(alternative=="two.sided"){
			pval<-"not calculated"
		}
	}
	if(test.type=="parametric"){
		res<-list(observed = moran.res$observed, expected = moran.res$expected, sd = sdi, p.value = pv)
	}
	if(test.type=="permutation"){
		res<-list(observed = boot.res$t0, expected = -1/(length(x) - 1), p.value = pval)
	}
	
	return(res)
}
