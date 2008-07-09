`find.angles` <-
function (coords) 
{
    n.angles <- (nrow(coords) * (nrow(coords) - 1))/2
    angles <- vector(length = n.angles)
    opp <- matrix(nrow = nrow(coords), ncol = nrow(coords))
    adj <- matrix(nrow = nrow(coords), ncol = nrow(coords))
    for (i in 1:nrow(coords)) {
        for (j in 1:nrow(coords)) {
            opp[i, j] <- coords[j, 2] - coords[i, 2]
            if (opp[i, j] != abs(opp[i, j])) {
                opp[i, j] <- NA
            }
            adj[i, j] <- abs(coords[j, 1] - coords[i, 1])
        }
    }
    opp.adj <- opp/adj
    angles <- pi/2 - atan(opp.adj)
    return(angles)
}

