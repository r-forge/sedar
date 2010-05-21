"scores.listw" <-
  function(listw,echo=FALSE, MEM.autocor = c("all","positive", "negative")) {
    if (!inherits(listw, "listw")) 
      stop("not a listw object")
    MEM.autocor <- match.arg(MEM.autocor)
    w <- listw2mat(listw)
    sumW <- sum(w)
    n <- nrow(w)
    symmetric <- isSymmetric.matrix(w,check.attributes=FALSE)
    if (symmetric==FALSE){
      cat (paste("listw not symmetric, (w+t(w)) used in the place of w","\n\n"))
      w=(w+t(w))/2
    }
    
    row.mean <- apply(w, 1, mean)
    col.mean <- apply(w, 2, mean)
    tot.mean <- mean(w)
    w <- sweep(w, 1, row.mean)
    w <- sweep(w, 2, col.mean)
    w <- w + tot.mean
    ## a voir w=w/max(w)
    
    res <- eigen(w,EISPACK=TRUE,symmetric=TRUE)
    eq0 <- apply(as.matrix(res$values/max(abs(res$values))),1,function(x) identical(all.equal(x, 0), TRUE))
    if(sum(eq0)==0) {stop(" Illegal matrix: no null eigenvalue")}
    if(echo) {cat(paste("vector number",which(eq0),"corresponding to null eigenvalue is removed","\n"))}
    res$values <- res$values[-which(eq0)]
    res$vectors <- res$vectors[,-which(eq0)]
    if(MEM.autocor == "positive"){
      
      posi <- which(res$values > -sumW/((n-1)*n))
      res$values <- res$values[posi]
      res$vectors <- res$vectors[,posi]
    }
    if( MEM.autocor == "negative"){
      neg <- sort(which(res$values < -sumW/((n-1)*n)), decreasing = TRUE)
      res$values <- res$values[neg]
      res$vectors <- res$vectors[,neg]
    }  
    res$call=match.call()
    res
  }
