"metacycle" <-
    function(steps=1, d, A, p=rep(1, length(A)), y=1, x=1, e = min(A), alpha=1,
             locations=NULL, ...)
{
    edis <- as.matrix(exp(-alpha*d))
    diag(edis) <- 0
    edis <- sweep(edis, 2, A, "*")
    E <- e/A^x
    E <- ifelse(E > 1, 1, E)
    if(is.null(locations))
        locations <- cmdscale(d)
    pmat <- matrix(0, nrow=length(p), ncol=steps + 1)
    pmat[,1] <- p
    for (i in 1:steps)
        pmat[,i+1] <- metastep(pmat[,i], edis, E, y)
    out <- list(p = pmat, d=edis, A=A, y=y, x=x, e=e, alpha=alpha, locations=locations)
    out$J.obs <- rowSums(pmat[,-1, drop=FALSE])/steps
    out$P.obs <- colSums(pmat)
    S <- rowSums(edis)
    C <- S^2/(S^2 + y)
    out$J.pot <- C/(C+E-C*E)
    out$S.pot <- S
    out$C.pot <- C
    class(out) <- "metacycle"
    out
}
