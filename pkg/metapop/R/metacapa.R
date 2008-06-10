`metacapa` <-
function(d, A, alpha=1, ...)
{
    M <- as.matrix(exp(-alpha*d))
    diag(M) <- 0
    M <- M * outer(A, A)
    tmp <- eigen(M)
    ev <- tmp$values[1]
    cap <- list(capacity = ev, x = tmp$vectors[,1]^2, d=d, A=A)
    class(cap) <- "metacapa"
    cap
}

