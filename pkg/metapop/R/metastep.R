`metastep` <-
function (p, edis, E, y) 
{
    S <- rowSums(edis[, p>0, drop=FALSE])
    C <- S^2/(S^2 + y^2)
    for (i in 1:length(p)) {
       if(p[i]==0 && runif(1) < C[i]) p[i] <- 1
       else if(p[i]==1 && runif(1) < (1-C[i])*E[i] ) p[i] <- 0
    }
    p
}

