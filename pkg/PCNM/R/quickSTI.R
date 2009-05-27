quickSTI <- function(Y, S, T, nperm=999, alpha = 0.05, COD.S=NULL, COD.T=NULL,print.res=TRUE)
#
# Space-time tests using alternative coding for space-time interaction (Model 5). Depending on the outcome of this test the 
# main factors are tested using different strategies. If the interaction is not significant then the test of main 
# factors is done using model 5. If the interaction is significant then a nested
# model is used to know whether separate spatial structures exist and another to know whether separate temporal structures exist.
# In this function space and time are considered fixed factors (F ratios are constructed using residual MS in the denominator).
#
# Pierre Legendre & Miquel De Caceres 
# February 2009
#
# PARAMETERS:
#
# Y     	is the site-by-species data table (assumes row blocks corresponding to times, i.e. within each block 
#			all sites occur)
# S     	contains the number of spatial points or a matrix of spatial coordinates
# T     	contains the number of time campaigns or a matrix (a vector) of temporal coordinates
# COD.S	contains the spatial coding functions to be used instead of s/2 PCNM. Their number must be lower than s.
# COD.T	contains the temporal coding functions to be used instead of t/2 PCNM. Their number must be lower than t.
# nperm 	the number of permutations to be done
# alpha	confidence level for the interaction test. Depending on the decision for this test the main factors are tested 
# 			differently.
#
# print.res 	if TRUE prints additional information
#
#
{
	if(!is.logical(print.res)) {
		stop("Wrong operator; 'print.res' should be either 'FALSE' or 'TRUE'")
	}

	aa <- system.time( {

	    # Sets the number and location of spatial and temporal points
		S = as.matrix(S)
		if(dim(S)[1]==1 && dim(S)[2]==1) {
			s=S[1,1]
			sitesX = c(1:s)
		} else {
			s = dim(S)[1]
			sitesX = S
		}		
		T = as.matrix(T)
		if(dim(T)[1]==1 && dim(T)[2]==1) {
			t=T[1,1]
			timesX = c(1:t) 
		} else {
			t = dim(T)[1]
			timesX = T
		}

		
		#Total number of rows in matrix
		n = s*t			
		
		# Check response data file containing species data
		Y = as.matrix(Y)
		p = dim(Y)[2]
		if(dim(Y)[1] != n) stop("The number of rows in species file is not (s x t)!") 

		# Center response data
		Y = apply(Y,2,scale,center=TRUE,scale=FALSE)
		
				
		if(print.res) {
			cat("=========================================================\n")
			cat("        Space-time ANOVA without replicates\n")
			cat("                                                  \n")
			cat("  Pierre Legendre, Miquel De Caceres, Daniel Borcard\n")
			cat("---------------------------------------------------------\n\n")

			cat(" No. space points (s) =", s,'\n')
			cat(" No. time points (t) =", t,'\n')
			cat(" No. observations (n = s*t) =", n,'\n')
			cat(" No. response variables (p) =", p,'\n')
			cat(" Significance level for the interaction test (alpha) =", alpha,'\n','\n')
			
		}

		# Generates spatial PCNMs if not given by user
		if(is.null(COD.S)) {			
			if(print.res) cat(" Computing PCNMs to code for space\n")
			nS = trunc(s/2)
			sitesX.D=dist(sitesX)
			pA = PCNM(sitesX.D, moran=FALSE, all=TRUE, include.zero=TRUE)
			npos = 0		# Compute how many positive eigenvalues
			for(i in 1:length(pA$values)) if(pA$values[i]>0) npos=npos+1
			# Do not get the eigenvector corresponding to the zero eigenvalue
			SS=as.matrix(pA$vectors[,-(npos+1)][,1:nS])
			PCNM.S=SS
			for(j in 2:t) PCNM.S = rbind(PCNM.S,SS)
		} else {
			PCNM.S=apply(as.matrix(COD.S),2,scale,center=TRUE,scale=TRUE)
			nS = dim(PCNM.S)[2]
			if(nS>=s) stop("The number of spatial coding functions must be lower than s!") 
		}
			
		# Generate PCNM variables for time if not given by user
		if(is.null(COD.T)) {
			nT = trunc(t/2)
			if(print.res) cat(" Computing PCNMs to code for time\n")
			timesX.D=dist(timesX)
			pB = PCNM(timesX.D, all=TRUE, moran=FALSE, include.zero=TRUE)
			npos = 0		# Compute how many positive eigenvalues
			for(i in 1:length(pB$values)) if(pB$values[i]>0) npos=npos+1
			# Do not get the eigenvector corresponding to the zero eigenvalue
			TT=as.matrix(as.matrix(pB$vectors[,-(npos+1)])[,1:nT])
			T.temp = TT[1,]
			for(i in 2:s) T.temp = rbind(T.temp,TT[1,])
			PCNM.T = as.matrix(T.temp)
			for(j in 2:t) {
				T.temp = TT[j,]
				for(i in 2:s) T.temp = rbind(T.temp,TT[j,])
				PCNM.T = as.matrix(rbind(PCNM.T,T.temp))
			}
		} else {
			PCNM.T=apply(as.matrix(COD.T),2,scale,center=TRUE,scale=TRUE)
			nT = dim(PCNM.T)[2]
			if(nT>=t) stop("The number of temporal coding functions must be lower than t!")
		}
		
		if(s*t-s-t-nT*nS-1<=0) stop("Not enough degrees of freedom for testing interaction!")


		if(print.res) {
			cat(" No. space coding functions =", nS,'\n')
			cat(" No. time coding functions =", nT,'\n\n')
		}

		# Generates space and time helmert contrasts
		A = as.factor(rep(1:s,t))
		B = rep(1,s)
		for(i in 2:t) B= c(B,rep(i,s))
		B = as.factor(B)
		HM = model.matrix(~ A + B, contrasts = list(A="contr.helmert", B="contr.helmert"))
		HM.S = as.matrix(HM[,2:s])
		HM.T = as.matrix(HM[,(s+1):(s+t-1)])
		
		
		#
		# Test significance for interaction effect 
		#
		# Defines X (variables for the factor of interest) and W (covariables) for the space-time test
		#
		XSTI = PCNM.S*PCNM.T[,1]
		if(dim(PCNM.T)[2]>1) for(j in 2:dim(PCNM.T)[2]) XSTI = cbind(XSTI,PCNM.S*PCNM.T[,j])
		if(print.res) {
			cat("------------------------------------------\n")
			cat(" Testing space-time interaction (model 5)\n")
			cat("------------------------------------------\n\n")
			cat("   No. space variables =", dim(HM.S)[2],'\n')
			cat("   No. time variables =", dim(HM.T)[2],'\n')
			cat("   No. interaction variables =", dim(XSTI)[2],'\n')
			cat("   No. residual degrees of freedom =", (s*t-dim(XSTI)[2]-dim(HM.S)[2]-dim(HM.T)[2]-1),"\n\n")
		}
				
		res = manovRDa(Y=Y,s=s,t=t,S.mat=HM.S,T.mat=HM.T,STI.mat=XSTI, S.test=TRUE, T.test=TRUE, STI.test=TRUE, model = "5", nperm=nperm)
		
		if(print.res) {
			cat(' Interaction test:  R2 =', res$testSTI$R2,'  F =',res$testSTI$F,'  P(',nperm,'perm) =',res$testSTI$Prob,"\n\n")
		}
						
		if(res$testSTI$Prob<=alpha) {
			XS = PCNM.S
			for(j in 1:(t-1)) XS = cbind(XS,PCNM.S*HM.T[,j])
			XT = HM.T
			if(print.res) {
				cat("----------------------------------------------------------------------\n")
				cat(" Testing for the existence of separate spatial structures (model 6a)\n")
				cat("----------------------------------------------------------------------\n\n")
				cat("   No. space variables =", dim(XS)[2],'\n')
				cat("   No. time variables =", dim(XT)[2],'\n')
				cat("   No. residual degrees of freedom =", (s*t-dim(XS)[2]-dim(XT)[2]-1),"\n\n")
			}			
			res2 = manovRDa(Y=Y,s=s,t=t,S.mat=XS,T.mat=XT,STI.mat=NULL, S.test=TRUE, T.test=FALSE, STI.test=FALSE, model="6a", nperm=nperm)
			res$testS = res2$testS
			cat(' Space test:   R2 =', res$testS$R2,'  F =',res$testS$F,'  P(',nperm,'perm) =',res$testS$Prob,"\n\n")
			XT = PCNM.T
			for(j in 1:(s-1)) XT = cbind(XT,PCNM.T*HM.S[,j])
			XS = HM.S
			if(print.res) {
				cat("-----------------------------------------------------\n")
				cat(" Testing for separate temporal structures (model 6b)\n")
				cat("-----------------------------------------------------\n\n")
				cat("   No. space variables =", dim(XS)[2],'\n')
				cat("   No. time variables =", dim(XT)[2],'\n')
				cat("   No. residual degrees of freedom =", (s*t-dim(XS)[2]-dim(XT)[2]-1),"\n\n")
			}			
			res2 = manovRDa(Y=Y,s=s,t=t,S.mat=XS,T.mat=XT,STI.mat=NULL, S.test=FALSE, T.test=TRUE, STI.test=FALSE, model="6b", nperm=nperm)
			res$testT = res2$testT
			cat(' Time test:   R2 =', res$testT$R2,'  F =',res$testT$F,'  P(',nperm,'perm) =',res$testT$Prob,"\n\n")					
		} else {
			if(print.res) {
				cat("---------------------------------------------------------------------\n")
				cat(" Testing for common spatial and common temporal structures (model 5)\n")
				cat("---------------------------------------------------------------------\n\n")
				cat(' Space test:   R2 =', res$testS$R2,'  F =',res$testS$F,'  P(',nperm,'perm) =',res$testS$Prob,'\n')
				cat(' Time test:   R2 =', res$testT$R2,'  F =',res$testT$F,'  P(',nperm,'perm) =',res$testT$Prob,'\n')
			}	
		} 		
		
	})
	
	aa[3] <- sprintf("%2f",aa[3])
	if(print.res) {
		cat("---------------------------------------------------------\n")
		cat("      Time for this calculus =",aa[3]," sec",'\n')
		cat("=========================================================\n\n")
	}
	invisible(res)
}