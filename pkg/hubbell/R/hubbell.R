"hubbell" <-
function(comm, D, m=0, P=NULL)
{
    J <- length(comm)
    wipe <- sample(J, D, replace=FALSE)
    comm[wipe] <- 0
    tmp.table <- table(comm)[-1]
    tmp.ind <- as.numeric(names(tmp.table))
    probs <- vector("numeric", max(tmp.ind))
    left <- sum(tmp.table)
    probs[tmp.ind] <- tmp.table/left
    invade <- runif(D) <= m
    aliens <- sum(invade)
    if (aliens > 0)
        comm[wipe[invade]] <- sample(length(P), aliens, replace=TRUE, prob=P)
    comm[wipe[!invade]] <- sample(max(tmp.ind), D-aliens, replace=TRUE, prob=probs)
    class(comm) <- "hubbell"
  comm
}
