"variogmultiv" <-function(Y,xy, dmin=0,dmax=max(dist(xy)),nclass=20){

    dxy <- seq(dmin,dmax,le=nclass+1)
    distgeo <- dist(xy)
    distfac <- cut(distgeo,bre=dxy)
    n.w <- table(distfac)
    colfac <- col(as.matrix(distgeo),as.factor=TRUE)[lower.tri(col(as.matrix(distgeo)))]
    rowfac <- row(as.matrix(distgeo),as.factor=TRUE)[lower.tri(row(as.matrix(distgeo)))]
    n.c <- apply(ifelse(table(distfac,colfac)+table(distfac,rowfac)>0,1,0),1,sum)

    Y2 <- dist(Y)^2
    dxy2 <- dxy[-1]
    dxy <- dxy[-(nclass+1)]
    
    
    res <- sapply(split(Y2,distfac),sum)/n.w/2
    
    results=list(d=(dxy+dxy2)/2,var=as.vector(res),n.c=as.vector(n.c),n.w=as.vector(n.w),dclass=levels(distfac))
    return(results)
}
