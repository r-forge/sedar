STImodels <- function(Y, S, T,  model="5", nperm=999, nS=-1,nT=-1, Sfixed=TRUE, Tfixed=TRUE, COD.S=NULL, COD.T=NULL, print.res=TRUE)
#
# Two way ANOVA models to test space-time interaction without replicates. Some degrees of freedom are saved by coding
# space and/or time parsimoniously. By default, the function computes PCNM for space and time coding, but the coding can be
# also provided by the user, through COD.S and COD.T.
#
# Tests for space-time interaction and space or time main effects are conducted using one of the different models.
# Depending on the model choosen, the interaction test is not available (concretely, with Models 2 and 6). 
# For the interaction the permutations are unrestricted, whereas for the main factors the permutations are restricted 
# within time or space blocks.
#
# Pierre Legendre & Miquel De Caceres 
# September 2007
#
# PARAMETERS:
#
# Y     		is the site-by-species data table (assumes row blocks corresponding to times, i.e. within each block 
#           	all sites occur)
# S     		the number of spatial points or a matrix of spatial coordinates
# T     		the number of time campaigns or a matrix (a vector) of temporal coordinates
# nperm 		the number of permutations to be done
#
# model 		linear space-time model to be used (now can be either "2", "3a", "3b", "4","5","6a", "6b", or "7")
#
# nS        	number of S PCNMs to use (by default -1 means s/2 rounding down)
# nT    		number of S PCNMs to use (by default -1 means t/2 rounding down)
# Sfixed		a logical value: is Factor Space fixed, or not (if FALSE, it is considered a random factor)
# Tfixed		a logical value: is Factor Time fixed, or not (if FALSE, it is considered a random factor)
# COD.S		the spatial coding functions to be used instead of PCNM. Their number must be lower than s.
# COD.T		the temporal coding functions to be used instead of PCNM. Their number must be lower than t.
#
# print.res 	if TRUE prints additional information
#
#
{
	if(!is.logical(print.res)) {
		stop("Wrong operator; 'print.res' should be either 'FALSE' or 'TRUE'")
	}
	if(model!="2" && model!="3a"&& model!="3b"&& model!="4"&& model!="4"&& model!="5"&& model!="6a"&& model!="6b"&& model!="7") {
		stop(paste("Unrecognized model ",model,"; 'model' should be '2', '3a', '3b', '4', '5', '6a','6b' or '7'.", sep=""))
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
		p = dim(Y)[2]

		# Check response data file containing species data
		Y=as.matrix(Y)
		p = dim(Y)[2]
		
		if(dim(Y)[1] != n) stop("The number of rows in species file is not (s x t)!") 

		# Center response data
		Y = apply(Y,2,scale,center=TRUE,scale=FALSE)

		
		if(print.res) {
			cat("=======================================================\n")
			cat("        Space-time ANOVA without replicates\n")
			cat("                                                  \n")
			cat("  Pierre Legendre, Miquel De Caceres, Daniel Borcard\n")
			cat("=======================================================\n\n")
			cat(" No. space points (s) =", s,'\n')
			cat(" No. time points (t) =", t,'\n')
			cat(" No. observations (n = s*t) =", n,'\n')
			cat(" No. response variables (p) =", p,'\n','\n')
		}

		# Generates space PCNM variables (if necessary)
		if(model=="3a"|| model=="6a" || model=="7"|| model=="4"|| model=="5") {
			if(is.null(COD.S)) {			# Generates spatial PCNMs if not given by user
				if(print.res) cat(" Computing PCNMs to code for space\n")
				if(nS==-1) nS = trunc(s/2)
				sitesX.D=dist(sitesX)
				pA = PCNM(sitesX.D, all=TRUE, moran=FALSE, include.zero=TRUE)
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
		}
		
		# Generate PCNM variables for time (if necessary)
		if(model=="3b" || model=="6b"|| model=="7"||model=="4" || model=="5") {
			if(is.null(COD.T)) {# Generate PCNM variables for time if not given by user
				if(print.res) cat(" Computing PCNMs to code for time\n")
				if(nT==-1) nT = trunc(t/2)
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
		}
		
		if(s*t-s-t-nT*nS-1<=0 && model!=2 && model!=6) stop("Not enough degrees of freedom for testing interaction!")


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
		
		test.STI = TRUE
		test.S = TRUE
		test.T = TRUE
		
		# Defines X (variables for the factor of interest) and W (covariables) for each test to be done
		if(model=="5") {
			XSTI = PCNM.S*PCNM.T[,1]
			if(dim(PCNM.T)[2]>1) for(j in 2:dim(PCNM.T)[2]) XSTI = cbind(XSTI,PCNM.S*PCNM.T[,j])
			XS = HM.S
			XT = HM.T
			if(print.res) {
				cat(" MODEL V: HELMERT CONTRAST FOR TESTING MAIN FACTORS. \n")
				cat("          SPACE AND TIME PCNMs FOR TESTING INTERACTION.",'\n')
			}					
		} else if(model=="4") {
			XSTI = PCNM.S*PCNM.T[,1]
			if(dim(PCNM.T)[2]>1) for(j in 2:dim(PCNM.T)[2]) XSTI = cbind(XSTI,PCNM.S*PCNM.T[,j])
			XS = PCNM.S
			XT = PCNM.T
			if(print.res) {
				cat(" MODEL IV: PCNMs FOR BOTH SPACE AND TIME.",'\n')
			}		
		} else if(model=="3a") {
			XSTI = PCNM.S*HM.T[,1]
			if(t>1) for(j in 2:(t-1)) XSTI = cbind(XSTI,PCNM.S*HM.T[,j])
			XS = PCNM.S
			XT = HM.T
			if(print.res) {
				cat(" MODEL IIIa: PCNMs FOR SPACE AND HELMERT CONTRASTS FOR TIME.",'\n')
			}	
		} else if(model=="3b") {
			XSTI = PCNM.T*HM.S[,1]
			for(j in 2:(s-1)) XSTI = cbind(XSTI,PCNM.T*HM.S[,j])
			XS = HM.S
			XT = PCNM.T
			if(print.res) {
				cat(" MODEL IIIb: HELMERT CONTRASTS FOR SPACE AND PCNMs FOR TIME.",'\n')
			}		
		} else if(model=="7") {
			XSTI = HM.S*HM.T[,1]
			if(dim(HM.T)[2]>1) for(j in 2:dim(HM.T)[2]) XSTI = cbind(XSTI,HM.S*HM.T[,j])
			XS = PCNM.S
			XT = PCNM.T
			if(print.res) {
				cat(" MODEL VII: PCNMs FOR BOTH SPACE AND TIME BUT HELMERT CONTRAST FOR INTERACTION.",'\n')
			}		
		} else if(model=="6a") {
			XS = PCNM.S
			for(j in 1:(t-1)) XS = cbind(XS,PCNM.S*HM.T[,j])
			XT = HM.T
			XSTI=NULL
			if(print.res) {
				cat(" MODEL VIa: NESTED MODEL.\n")
				cat("           TESTING FOR THE EXISTENCE OF SPATIAL STRUCTURE (COMMON OR SEPARATED)",'\n')
			}	
			test.STI = FALSE
		} else if(model=="6b") {
			XT = PCNM.T
			for(j in 1:(s-1)) XT = cbind(XT,PCNM.T*HM.S[,j])
			XS = HM.S						
			XSTI=NULL
			if(print.res) {
				cat(" MODEL VIb: NESTED MODEL.\n")
				cat("           TESTING FOR THE EXISTENCE OF TEMPORAL STRUCTURE (COMMON OR SEPARATED).",'\n')
			}	
			test.STI = FALSE
		} else if(model=="2") {
			XS = HM.S
			XT = HM.T
			XSTI=NULL
			if(print.res) {
				cat(" MODEL II: HELMERT CONTRAST FOR SPACE AND TIME. NO INTERACTION TERM.",'\n')
			}	
			test.STI = FALSE
		} 		
			if(print.res) {
				cat("   No. space variables =", dim(XS)[2],'\n')
				cat("   No. time variables =", dim(XT)[2],'\n')
				if(test.STI) {
					cat("   No. interaction variables =", dim(XSTI)[2],'\n')
					cat("   No. residual degrees of freedom =", (s*t-dim(XS)[2]-dim(XT)[2]-dim(XSTI)[2]-1),"\n\n")
				} else {
					cat("   No. residual degrees of freedom =", (s*t-dim(XS)[2]-dim(XT)[2]-1),"\n\n")
				}
			}		
		
		
		res = manovRDa(Y=Y,s=s,t=t,S.mat=XS,T.mat=XT,STI.mat=XSTI, Sfixed= Sfixed, Tfixed=Tfixed, S.test=test.S, T.test=test.T, STI.test=test.STI, model = model, nperm=nperm)
		
		if(test.STI==TRUE) {
			cat(' Interaction test:   R2 =', res$testSTI$R2,'  F =',res$testSTI$F,'  P(',nperm,'perm) =',res$testSTI$Prob,'\n')
		}
		if(test.S==TRUE) {
			cat(' Space test:   R2 =', res$testS$R2,'  F =',res$testS$F,'  P(',nperm,'perm) =',res$testS$Prob,'\n')
		}
		if(test.T==TRUE) {
			cat(' Time test:   R2 =', res$testT$R2,'  F =',res$testT$F,'  P(',nperm,'perm) =',res$testT$Prob,'\n')
		}
		
	})
	aa[3] <- sprintf("%2f",aa[3])
	if(print.res) {
		cat("-------------------------------------------------------\n")
		cat("      Time for this calculation =",aa[3]," sec",'\n')
		cat("=======================================================\n\n")
	}
	invisible(res)
}