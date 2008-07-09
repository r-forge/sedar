moran.I.uni<-function (x, mat.W, scaled = FALSE,normalize=FALSE, na.rm = FALSE, alternative = "two.sided") {
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
	
	#======================
	#### Claculate Moran's I
	#======================
	#CC# S0 (Gittleman & Kot, eq. 2)
	s <- sum(mat.W)
	
	#CC# Calculate centred y
	m <- mean(x)
	y <- x - m
	
	#CC# Numerator
	cv <- sum(mat.W * y %o% y)
	
	#CC# Denominator
	v <- sum(y^2)
	
	#CC# Observed Moran's I (Gittleman & Kot, eq. 1)
	obs <- (n/s) * (cv/v)
	
	#CC# Scaling Moran's I (Gittleman & Kot, eq. 4)
	if (scaled) {
		i.max <- (n/s) * (sd(rowSums(mat.W) * y)/sqrt(v/(n - 1)))
		obs <- obs/i.max
	}
	
	#============================================
	#### Calculate Standard deviation of Moran's I
	#============================================
	#CC# (Gittleman & Kot, eq. 7)
	S1 <- 0.5 * sum((mat.W + t(mat.W))^2)
	
	#CC# (Gittleman & Kot, eq. 8)
	S2 <- sum((apply(mat.W, 1, sum) + apply(mat.W, 2, sum))^2)
	
	#CC# (Gittleman & Kot, eq. 6)
	s.sq <- s^2
	k <- (sum(y^4)/n)/(v/n)^2
	sdi <- sqrt((n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * s.sq) - 
		k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * s.sq))/((n - 
		1) * (n - 2) * (n - 3) * s.sq) - 1/((n - 1)^2))
	
	#CC# Mean of Moran's I (Gittleman & Kot, eq. 5)
	ei <- -1/(n - 1)

	#===============================
	#### Parametric test of Moran's I
	#===============================
	alternative <- match.arg(alternative, c("two.sided", "less", 
		"greater"))
	pv <- pnorm(obs, mean = ei, sd = sdi)
	if (alternative == "two.sided") 
		pv <- if (obs <= ei) 
			2 * pv
		else 2 * (1 - pv)
	if (alternative == "greater") 
		pv <- 1 - pv
	list(observed = obs, expected = ei, sd = sdi, p.value = pv)
}
