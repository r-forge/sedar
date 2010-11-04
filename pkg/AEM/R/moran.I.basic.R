moran.I.basic<-function(x,mat.W,scaled){
	#CC# number of sites
	n<-length(x)
	
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
	
	#CC# Mean of Moran's I (Gittleman & Kot, eq. 5)
	ei <- -1/(n - 1)

	return(list(observed=obs,expected=ei,s=s,y=y,v=v))
}
