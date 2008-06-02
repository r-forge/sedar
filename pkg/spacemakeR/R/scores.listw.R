"scores.listw" <-
function(listw,echo=FALSE) {
    if (!inherits(listw, "listw")) 
          stop("not a listw object")
    w <- listw2mat(listw)
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
    # a voir w=w/max(w)
    
    res <- eigen(w,EISPACK=TRUE,symmetric=TRUE)
    eq0 <- apply(as.matrix(res$values/max(abs(res$values))),1,function(x) identical(all.equal(x, 0), TRUE))
    if(sum(eq0)==0) {stop(" Illegal matrix: no null eigenvalue")}
    if(echo) {cat(paste("vector number",which(eq0),"corresponding to null eigenvalue is removed","\n"))}
    res$values <- res$values[-which(eq0)]
    res$vectors <- res$vectors[,-which(eq0)]
    res$call=match.call()
    res
    
}
