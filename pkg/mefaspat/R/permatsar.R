`permatsar` <-
function(m, nb, pd=NULL, mtype="count", times=100, round.sar=TRUE)
{
    require(vegan)
    vspec <- apply(m>0, 2, sum)
    if (any(vspec==0)) stop("all species must have > 0 column totals")
    require(prabclus)
    mtype <- match.arg(mtype, c("prab", "count"))
    ini <- prabinit(prabmatrix=t(m), neighborhood=nb)
    if (is.null(pd)) pd <- autoconst(ini)$pd
    sarestimate <- prab.sarestimate(ini)
    perm <- list()
    n.row <- nrow(m)
    n.col <- ncol(m)
    if (mtype == "prab")
        for (i in 1:times)
            perm[[i]] <- randpop.nb(nb, p.nb=pd, n.species=n.col, n.regions=n.row, 
            vector.species=vspec, species.fixed=TRUE, count=FALSE)
    if (mtype == "count") {
        for (i in 1:times) {
            perm[[i]] <- regpop.sar(ini, p.nb=pd, sarestimate=sarestimate)
            if (round.sar) {
                perm[[i]] <- round(perm[[i]], 0)
                d <- sum(m) - sum(perm[[i]])
                if (d != 0) {
                while(d != 0) {
                    if (d < 0) {
                        ch <- -1
                        th <- 1
                        } else {
                        ch <- 1
                        th <- 0}
                    rr <- sample(n.row, 1)
                    rc <- sample(n.col, 1)
                    if (perm[[i]][rr,rc] > th)
                        perm[[i]][rr,rc] <- perm[[i]][rr,rc] + ch
                    d <- sum(m) - sum(perm[[i]])
                    }}
                }
            }}
    specs <- list(nb=nb, pd=pd, sar=sarestimate)
    out <- list(call=match.call(), orig=m, perm=perm, specs=specs)
    attr(out, "mtype") <- mtype
    attr(out, "ptype") <- "sar"
    attr(out, "method") <- NA
    attr(out, "fixedmar") <- "none"
    attr(out, "times") <- times
    attr(out, "shuffle") <- NA
    class(out) <- c("permat", "list")
    return(out)
}

