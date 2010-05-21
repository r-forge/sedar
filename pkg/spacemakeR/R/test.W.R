"test.W" <-
  function(Y,nb,xy,MEM.autocor = c("all","positive", "negative"),f=NULL, ...){
    mycall <- pairlist(...)
    res <- list()
    MEM.autocor <- match.arg(MEM.autocor)
    if(!(is.null(f))){
      nbdist <- nbdists(nb,as.matrix(xy))
      if(!(is.null(mycall))){
        param <- expand.grid(as.list(mycall))
        m1 <- match(names(param),names(formals(f)))
        for (i in 1:nrow(param)){
          formals(f)[m1] <- unclass(param[i,])
          res[[i]] <- scores.listw(nb2listw(nb,style="B",glist=lapply(nbdist,f)), MEM.autocor = MEM.autocor)
        }
      }
      else {res[[1]] <- scores.listw(nb2listw(nb,style="B",glist=lapply(nbdist,f)), MEM.autocor = MEM.autocor)}
    }
    else { res[[1]] <- scores.listw(nb2listw(nb,style="B"), MEM.autocor = MEM.autocor) }
    res2 <- lapply(res,function(x) ortho.AIC(Y=Y,X=x$vec, ord.var=TRUE))
    if(!(is.null(mycall))){
      res3 <- data.frame(AICc=unlist(lapply(res2,function(x) min(x[[1]],na.rm=TRUE))),NbVar= unlist(lapply(res2,function(x) which.min(x[[1]]))))
      res3 <- cbind(param,res3)
    }
    else{
      res3 <- data.frame(AICc=unlist(lapply(res2,function(x) min(x[[1]],na.rm=TRUE))),NbVar= unlist(lapply(res2,function(x) which.min(x[[1]]))))
    }

    thebest <- which.min(res3$AICc)
    cat (paste("\n\nAICc for the null model:",res2[[thebest]]$AICc0,"\n"))
    cat ("\nBest spatial model:\n")
    print(res3[thebest,])
    
    return(list(all=res3,best=c(res[[thebest]],res2[[thebest]])))
    
}
