#################################################
#
# Functions to perform a Multi Scale
# Ordination of Variables (MSPA).
# 
# T. Jombart (jombart@biomserv.univ-lyon1.fr)
# 2007
#################################################





#####################
# Function mspa
#####################
mspa <- function(dudi, lw, scannf = TRUE, nf = 2, centring=c("param","sim"), nperm=1000){
  
  # arguments checks
  if(!require(ade4)) stop("package ade4 is required")
  if(!require(spdep)) stop("package spdep is required")
  if (!inherits(dudi, "dudi")) stop("object of class 'dudi' expected")
  if(!inherits(lw,"listw")) stop("lw must be a listw object (package spdep)")
  if(lw$style != "W") stop("lw must be row-nowmalized (style W)")

  df <- dudi$tab
  varweights <- dudi$cw/sum(dudi$cw)
  n <- nrow(df)
  p <- ncol(df)
  centring <- match.arg(centring)
  
  ## auxiliary variables and functions
  
  # retrieve var.idx with correct variable names
  findname <- function(vec){
    if(length(vec) == 1) return(vec)
    res <- sub("[.][^.]*$","",vec[1])
    return(res)
  }

  # var.idx
  var.idx <- dudi$assign
  if(!is.null(var.idx)){
    temp <- split(colnames(df), var.idx)
    newlev <- sapply(temp, findname)
    levels(var.idx) <- newlev
  }
  
  genlab <- function(base,n){
    f1 <- function(cha, n) {
      if (nchar(cha) < n) {
        cha <- paste("0", cha, sep = "")
        return(f1(cha, n))
      }      
      return(cha)
    }
    w <- as.character(1:n)
    max0 <- max(nchar(w))
    w <- sapply(w, function(cha) f1(cha, max0))
    return(paste(base, w, sep = ""))
  } # end genlab

  
  # matrix centring and scaling
  X <- scalewt(df, center=TRUE, scale=TRUE)
  
  # computation of centred and scaled eigenvectors of the connection network
  # denoted U
  U <- as.matrix(orthobasis.listw(lw))
  r <- ncol(U)
  colnames(U) <- genlab("u_",r)  

  ## projection of X onto U
  ## only R-squared of each vector of U are kept
 
  ## model computations
  R <- t(X) %*% U /n
  R2 <- R*R

  ## handle centring of R2
  if(centring=="param"){
      newdf <- as.data.frame(R2-(1/(n-1)))
  } else{ # i.e. if centring is non-parametric
      tempX <- t(X)
      fPerm <- function(X){
          permX <- X[, sample(1:n)] # X has to be transposed, that is, variables in rows, obs in columns
          res <- permX %*% U /n
          res <- res * res
          return(res)
      } # end fPerm
      listR2sim <- lapply(1:nperm, function(i) fPerm(tempX))
      meanR2sim <- listR2sim[[1]]
      for(i in 2:nperm){
          meanR2sim <- meanR2sim + listR2sim[[i]]
      } # end for

      meanR2sim <- meanR2sim / nperm

      newdf <- as.data.frame(R2 - meanR2sim)
  } # end non-parametric centring

  ## keep only positive deviation in R2
  newdf[newdf < 0] <- 0

  ## we proceed to the analysis of this matrix
  res <- as.dudi(newdf, scannf=scannf, nf=nf, row.w = varweights,
                 col.w = rep(1, r), call=match.call(), type="mspa")

  res$ls <- as.data.frame(as.matrix(R2) %*% as.matrix(res$c1))  
  colnames(res$ls) <- colnames(res$li)
  row.names(res$ls) <- row.names(res$li)

  xmoy <- apply(R2, 2, function(c) weighted.mean(c,varweights))
  bary <- as.vector(t(xmoy) %*% as.matrix(res$c1))
  names(bary) <- colnames(res$c1)

  res$R2 <- R2
  res$meanPoint <- bary
  res$varweights <- varweights
  names(res$varweights) <- colnames(X)
  if(!is.null(var.idx)) res$assign <- var.idx
  if(centring=="sim") {
      res$centring <- meanR2sim
  }
  
  # return result
  return(res)
}





#############################
# function scatter.mspa
#############################
scatter.mspa <- function(x, xax = 1, yax = 2, clab.var = 0.75, clab.sca = 1, 
    posieig = "top", sub = NULL, ratio = 1/4, bary=TRUE, circle=TRUE, ...){

  if(!inherits(x,"mspa")) stop("to be used with mspa objects only")

  opar <- par(mar = par("mar"))
  on.exit(par(opar))
  
  s.arrow(x$c1[,c(xax,yax)], clab=clab.sca, sub=sub, ...)
  
  extrem <- chull(x$c1[,c(xax,yax)])

  par(xpd=TRUE)
  polygon(x$c1[extrem,c(xax,yax)],col="lightgrey")
  if(circle) {symbols(x=0,y=0,circles=1,add=TRUE,inches=FALSE)}

  s.arrow(x$c1[,c(xax,yax)], clab=clab.sca, add.p=TRUE)  
  s.label(x$ls[,c(xax,yax)], clab=clab.var, add.p=TRUE)
  
  # compute coordinates of factors
  # = weighted mean of its levels
  if(!is.null(x$assign)) {
    temp <- apply(x$li[,c(xax,yax)],2,function(c) tapply(c, x$assign,mean))
    rownames(temp) <- unique(x$assign)
    s.label(temp, add.plot=TRUE, clab=clab.var)
  }

  if(bary) points(x$meanPoint[xax],x$meanPoint[yax],pch=20,cex=2)
  
  add.scatter.eig(x$eig, x$nf, xax, yax, posi = posieig, ratio = ratio)

  par(mar=rep(.1,4))
  box(which="figure")

  return(invisible(match.call))
}
