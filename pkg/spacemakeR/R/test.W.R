"test.W" <-
function(Y,nb,xy,f=NULL,...){
    mycall <- pairlist(...)
    res <- list()
    if(!(is.null(f))){
        nbdist <- nbdists(nb,as.matrix(xy))
        if(!(is.null(mycall))){
            param <- expand.grid(as.list(mycall))
            m1 <- match(names(param),names(formals(f)))
            for (i in 1:nrow(param)){
                formals(f)[m1] <- unclass(param[i,])
                res[[i]] <- scores.listw(nb2listw(nb,style="B",glist=lapply(nbdist,f)))
            }
        }
        else {res[[1]] <- scores.listw(nb2listw(nb,style="B",glist=lapply(nbdist,f)))}
    }
    else { res[[1]] <- scores.listw(nb2listw(nb,style="B")) }
    res2 <- lapply(res,function(x) ortho.AIC(Y=Y,X=x$vec, ord.var=TRUE))
    if(!(is.null(mycall))){
        res3 <- data.frame(AICc=unlist(lapply(res2,function(x) min(x[[1]],na.rm=TRUE))),NbVar= unlist(lapply(res2,function(x) which.min(x[[1]]))))
        res3 <- cbind(param,res3)

        }
    else{ res3 <- data.frame(AICc=unlist(lapply(res2,function(x) min(x[[1]],na.rm=TRUE))),NbVar= unlist(lapply(res2,function(x) which.min(x[[1]]))))}
   
    cat ("\n\nBest model:\n\n\n")
    thebest <- which.min(res3$AICc)
    print(res3[thebest,])

    return(list(all=res3,best=c(res[[thebest]],res2[[thebest]])))
    
}
