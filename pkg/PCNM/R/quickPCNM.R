'quickPCNM' <- 
	function(Y, space, thresh=NULL, method="fwd", myPCNM=NULL, alpha=0.05, rangexy=FALSE, detrend=TRUE, perm.max=NULL, original = FALSE)
{
#                                        Daniel Borcard
#                                        Universite de Montreal
#                                        May 2007 - April 2010

require(ade4)   # Needed for function 's.value'
require(vegan)

a <- system.time({

Y <- as.matrix(Y)
space <- as.matrix(space)
n <- nrow(Y)
epsilon <- sqrt(.Machine$double.eps)

# ----------------------------------------------------------------------------

if(rangexy==TRUE) {

if(ncol(space)==1){
space <- (space-min(space))/((max(space)-min(space))*0.1)
                  }
else{
mini <- apply(space,2,"min")
maxi <- apply(space,2,"max")
xy.trans <- sweep(space,2,mini)
range.max <- 0.1*(max((maxi[1]-mini[1]),(maxi[2]-mini[2])))
space <- as.matrix(xy.trans/range.max)
    }
                }

# ----------------------------------------------------------------------------

if (is.null(myPCNM)) {

### Building the PCNM variables

prePCNM <- PCNM(dist(space),thresh=thresh, dbMEM=FALSE, moran=TRUE, all=FALSE, include.zero=FALSE)

dmin <- prePCNM$thresh

positMoran <- which(prePCNM$Moran_I[,3]==TRUE)

if(original==TRUE) {
ev <- prePCNM$values
nb.ev <- length(which(ev > epsilon))
PCNMbase <- as.data.frame(prePCNM$vectors[1:nrow(space),])
                    }

else {
ev <- prePCNM$values[positMoran]
nb.ev <- length(which(ev > epsilon))
PCNMbase <- as.data.frame(prePCNM$vectors[1:nrow(space),positMoran])
      }
assign("PCNMbase", PCNMbase, envir=.GlobalEnv)# compensation for internal loss of object by R
                     }
                     
else {
PCNMbase <<- as.data.frame(myPCNM)
dmin <- "implicit in PCNM file"
meanPCNM <- sum(apply(PCNMbase,2,mean))
if(abs(meanPCNM) > epsilon ) {
cat("\n------------------------------------------------------------------")
cat("\nWARNING: the user-provided spatial variables are not centred")
cat("\nto zero mean. Are you sure that they are correct?")
cat("\n------------------------------------------------------------------")
                                  }
sumcorPCNM <- sum(cor(PCNMbase))
if(abs(sumcorPCNM) > ncol(PCNMbase)+epsilon ) {
cat("\n------------------------------------------------------------------")
cat("\nWARNING: the user-provided spatial variables are not orthogonal")
cat("\nto one another. Are you sure that they are correct?")
cat("\n------------------------------------------------------------------")
cat("\n")
                                    }
nb.ev <- ncol(PCNMbase)
ev <- apply(PCNMbase^2,2,sum)
     }
if(nb.ev >= n-1)  stop("Too many explanatory variables. They should be less than n-1",call.=FALSE)

# ----------------------------------------------------------------------------

### RDA "response matrix x PCNM"

## Preliminary step: linear detrending of response data if trend is significant
if(detrend==TRUE){
   trace.correl <-sum(diag(cor(Y,space)))
   if(abs(trace.correl)<1e-12){
     Y.det <- Y
     assign("Y.det", Y.det, envir=.GlobalEnv)  # compensation for internal loss of object by R
     temp2.test<-matrix(rep(1,5),1)
                              }
   else{  
   temp2 <- vegan::rda(Y,space)
   temp2.test <- vegan::anova.cca(temp2,alpha=alpha,step=100,perm.max=perm.max)
   if(temp2.test[1,5] <= alpha) {
      Y.det <- resid(lm(Y~space))
     assign("Y.det", Y.det, envir=.GlobalEnv)   # compensation for internal loss of object by R
       }                 
   else {
      Y.det <- Y
      assign("Y.det", Y.det, envir=.GlobalEnv)   # compensation for internal loss of object by R
      temp2.test<-matrix(rep(1,5),1)
        }
                 }
                              }
else {
   Y.det <- Y
   assign("Y.det", Y.det, envir=.GlobalEnv)  # compensation for internal loss of object by R
     }

## RDA with all PCNM variables(complete model)

mod1 <- vegan::rda(Y.det~.,data=PCNMbase)
if(is.null(perm.max)) {
	global.test <- vegan::anova.cca(mod1,alpha=alpha)
	} else {
	global.test <- vegan::anova.cca(mod1,alpha=alpha,step=100,perm.max=perm.max)
	}
mod1.sum <- summary(mod1,scaling=1)
R2glob <- mod1.sum$constr.chi/mod1.sum$tot.chi
R2glob.a <- 1-((n-1)/(n-global.test[1,1]-1))*(1-R2glob)

if(global.test[1,5] >= alpha) {
   cat("\n------------------------------------------------------------------")
   cat ("\n*** Procedure stopped ***")
   cat ("\np-value of global test: ",global.test[1,5])
   cat ("\nNo significant spatial structure detected by global PCNM analysis")
   cat ("\nSelection of PCNM variables would lead to spurious model")
   cat("\n------------------------------------------------------------------","\n")
     }

# ----------------------------------------------------------------------------

else {

METHODS <- c("none", "fwd")
    method <- match.arg(method, METHODS)

if(method == "none"){ 
mod <- mod1
mod.test <- global.test
mod.sum <- mod1.sum
R2 <- R2glob
R2adj <- R2glob.a
vars.sign <- c(1:ncol(PCNMbase))
nb.sig.ev <- length(vars.sign)
                    }

else{

  if(method == "fwd"){   # 1.2.1 open  fwd

  ## Regression-based forward selection of PCNM variables and RDA on significant
  ##     PCNM variables

  require(packfor)  #  packfor version 0.0-7 or later (to allow adjR2thresh)
   
  # If there is only one response variable that is normally distributed, save
  # time by replacing permutation tests by parametric tests. Otherwise and for
  # a multivariate response matrix, forward selection with the adjusted R2 as 
  # aditional stopping criterion.
    if(ncol(Y)==1){
    norm.test <- shapiro.test(Y.det)
       if(norm.test[2] > 0.05) {
       cat("\n------------------------------------------------------------------")
       cat("\nOnly one, normally distributed response variable found. Parametric")
       cat("\nforward selection on standardized response variable has been run.")
       cat("\n------------------------------------------------------------------","\n")
       Y.det <-scale(Y.det)
       fwd.sel <- packfor::forward.sel.par(Y.det,PCNMbase,alpha=alpha,adjR2thresh=R2glob.a)
                               } else {
       cat("\n------------------------------------------------------------------")
       cat("\nThe only response variable is not normally distributed.")
       cat("\nPermutational forward selection has been run.")
       cat("\n------------------------------------------------------------------","\n")
       fwd.sel <- packfor::forward.sel(Y.det,PCNMbase,alpha=alpha,adjR2thresh = R2glob.a)
                                      }
                  } else {
  fwd.sel <- packfor::forward.sel(Y.det,PCNMbase,alpha=alpha,adjR2thresh = R2glob.a)
                         }

# Compensation for a bug in packfor: if the stopping criterion is adjR2thresh, 
# one nonsignificant variable is retained.

          if(fwd.sel[nrow(fwd.sel),5] > R2glob.a + epsilon) {
             fwd.sel <- fwd.sel[-nrow(fwd.sel),]       }

  nb.sig.ev <- nrow(fwd.sel)
  vars.sign <- sort(fwd.sel[,2])}   # 1.2.1 close  fwd


PCNMred <- as.data.frame(PCNMbase[,c(vars.sign)])
assign("PCNMred",PCNMred,envir=.GlobalEnv)  # compensation for internal loss of object by R
mod <- vegan::rda(Y.det~.,data=PCNMred)
mod.sum <- summary(mod,scaling=1)}

if(is.null(perm.max)) {
	mod.test <- vegan::anova.cca(mod, alpha=alpha)
	} else {
	mod.test <- vegan::anova.cca(mod, alpha=alpha, step=100, perm.max=perm.max)
	}

R2 <- mod.sum$constr.chi/mod.sum$tot.chi
R2adj <- 1-((n-1)/(n-mod.test[1,1]-1))*(1-R2)

   }                        # close all choices (incl. no selection)

# Warning if the adjusted R-square of minimal model is greater than the
# adjusted R-square of the global model with all PCNM variables.

# if(R2adj > (1.05*R2glob.a)){
# cat("\n------------------------------------------------------------------")
# cat("\nWARNING: the adjusted R-square of the reduced model,",round(R2adj,4))
# cat("\nexceeds the adjusted R-square of the global model,",round(R2glob.a,4))
# cat("\nby more that 5%.This means that the selection has been overly liberal.")
# cat("\nChoose another, mode conservative method (see explanations).")
# cat("\n------------------------------------------------------------------")
#                                    }
# else if(R2adj > (R2glob.a) & R2adj <= (1.05*R2glob.a)){
# if(R2adj > (R2glob.a) & R2adj <= (1.05*R2glob.a)){
# cat("\n------------------------------------------------------------------")
# cat("\nWARNING: the adjusted R-square of the reduced model,",round(R2adj,4))
# cat("\nexceeds the adjusted R-square of the global model,",round(R2glob.a,4))
# cat("\nby 5% or less. This means that the selection has been a little bit")
# cat("\ntoo liberal. This small amount should not harm, however.")
# cat("\n------------------------------------------------------------------")
#                                    }

# At this point, anova.cca for all axes doesn't seem to work properly: it looses
# the track of one or the other object defined higher (generally PCNMred).
# The workaround is to export PCNMred from the function during the run
# and assign it to the global R environment (Line 247). PCNMred is thus present
# in the main R workspace and should be deleted prior to another quickPCNM run.

if(ncol(Y)==1 || nb.sig.ev == 1) {
   nb.ax <- 1                    }
 else {
 if(is.null(perm.max)) {
	 mod.axes.test <- vegan::anova.cca(mod, by="axis", cutoff=0.10)
	 } else {
	 mod.axes.test <- vegan::anova.cca(mod,by="axis",step=100,perm.max=perm.max,cutoff=0.10)
	 }
 ## Count the significant axes 
nb.ax <- length(which(mod.axes.test[,5]<=alpha))

        }         

# ----------------------------------------------------------------------------

## Number of axes to draw (arbitrary rule, an alterative to the tests above)
#  if(mod.test[1,1]<=2){
#     nb.ax=mod.test[1,1]} else {
#     nb.ax=ceiling(sqrt(mod.test[1,1]))}

## Plot of significant axes
fitted.scores=vegan::scores(mod,display="lc",choices=1:nb.ax)
par(mfrow=c(round(sqrt(nb.ax)),ceiling(sqrt(nb.ax))))

if(ncol(space)==2){
   for(i in 1:nb.ax){
   ade4::s.value(space,fitted.scores[,i],addaxes=FALSE,include.origin=FALSE,sub=paste("Axis ",i),csub=1.5)
                    }
} else {
   for(i in 1:nb.ax){
   plot(space,fitted.scores[,i],type="l",ylab="Fitted site scores")
                    }
       }

# ----------------------------------------------------------------------------

## Screen output
if(detrend==TRUE){
   if(temp2.test[1,5] <= alpha) {
      cat("\n-------------------------------------------------------")
      cat("\nA significant linear trend has been found in the response data.")
      cat("\nThe data have been detrended prior to PCNM analysis.")
                                }
    else{
      cat("\n-------------------------------------------------------")
      cat("\nNo significant linear trend has been found in the response data.")
      cat("\nThe data have NOT been detrended prior to PCNM analysis.")     
        }
                 }


cat("\n-------------------------------------------------------")
cat("\nThe truncation value used for PCNM building is",dmin,"\n")
cat(nb.ev," PCNM eigenvectors have been produced","\n")
cat("Adjusted R2 of global model = ",round(R2glob.a,4),"\n")
if(method != "none") {
if(nb.sig.ev==1){
cat(nb.sig.ev," PCNM variable has been selected","\n")}
else{
cat(nb.sig.ev," PCNM variables have been selected","\n")}
cat("R2 of minimum model = ",round(R2,4),"                    ","\n")
cat("Adjusted R2 of minimum model = ",round(R2adj,4),"        ","\n")
                     }
# cat("The minimum model has ",nb.ax,"significant canonical axes","\n")
cat("---------------------------------------------------------")
cat("\n")

})
a[3] <- sprintf("%2f",a[3])
cat("Time to compute quickPCNM =",a[3]," sec",'\n')

## Extraction of results to be returned

if(method == "none") {
   if(ncol(Y)>1 && nb.sig.ev > 1) {
#      table <- list(PCNMbase, ev[1:nb.ev], mod.sum, mod.test, mod.axes.test)
      table <- list(PCNMbase, ev[1:nb.ev], mod, mod.test, mod.axes.test)
      names(table) <- c("PCNM","PCNM_eigenvalues","RDA","RDA_test","RDA_axes_tests")
      } else {
#      table <- list(PCNMbase,ev[1:nb.ev],mod,mod.test)
      table <- list(PCNMbase,ev[1:nb.ev],mod,mod.test)
      names(table) <- c("PCNM","PCNM_eigenvalues","RDA","RDA_test")   
      }
} else {

   if(method == "fwd") {
      if(ncol(Y)>1 && nb.sig.ev > 1) {
#   table <- list(PCNMbase, ev[1:nb.ev], fwd.sel, PCNMred, mod.sum, mod.test, mod.axes.test)
   table <- list(PCNMbase, ev[1:nb.ev], fwd.sel, PCNMred, mod, mod.test, mod.axes.test)
   names(table) <- c("PCNM","PCNM_eigenvalues","fwd.sel","PCNM_reduced_model","RDA","RDA_test","RDA_axes_tests")
       } else {
#   table <- list(PCNMbase,ev[1:nb.ev], fwd.sel, PCNMred, mod.sum, mod.test)
   table <- list(PCNMbase,ev[1:nb.ev], fwd.sel, PCNMred,mod, mod.test)
   names(table) <- c("PCNM","PCNM_eigenvalues","fwd.sel","PCNM_reduced_model","RDA","RDA_test")
       }

                      }
    }
return(table)
}
