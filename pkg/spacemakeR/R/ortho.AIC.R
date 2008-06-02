"ortho.AIC" <-
function(Y,X,ord.var=FALSE){
# Fast Forward Selection AIC if X is orthonormal (XtX=I)
# return a vector of AICc
# if ord.var=TRUE, a list containing also order of variables is returned
if (sum(apply(as.matrix(apply(X,2,mean)),1,function(x) identical(all.equal(x, 0), TRUE)))!=ncol(X)) stop("X variables are not centered")
if (!(sum(identical(all.equal(sum(crossprod(as.matrix(X))-diag(ncol(X))), 0), TRUE)))) stop("X variables are not orthonormalized")

    f1 <- function(resp,X){
        R2 <- t(as.matrix(X))%*%as.matrix(Y)
        R2 <- t(R2)%*%R2
        return(sum(diag(R2)))
    }
    
    Y <- scale(Y,scale=FALSE)
    Y <- as.matrix(Y)
    R2 <- apply(X,2,f1,resp=Y)
    SSTot <- sum(diag(t(Y)%*%Y))
    RSS <- SSTot-R2
    ordre <- order(RSS)
    RSS <- sort(RSS)
    R2 <- R2[ordre]
    RSScum <- cumsum(c(RSS[1],-R2[-1]))
    RSScum <- c(SSTot,RSScum)
    # By default, Y is centered
    # K is the number of othogonal vectors + 1 (for intercept)
    
    K <- (1+(0:ncol(X)))
    AICtri <- nrow(X)*log(ifelse(RSScum<=0,NA,RSScum/nrow(X)))+2*K
    correct <- 2*(K*(K+1))/(nrow(X)-K-1)
    correct <- ifelse(is.finite(correct)&(correct>0),correct,NA)
    AICc <- AICtri+correct
    if(ord.var) {AICc <- list(AICc=AICc[-1],AICc0=AICc[1],ord=ordre,R2=cumsum(R2/SSTot))}
    if(!ord.var) AICc <- AICc[-1]
    return(AICc)
}
