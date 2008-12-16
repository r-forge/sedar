"give.thresh" <-
function(distxy){
    spanning <- ade4::mstree(distxy)
    return(max(neig2mat(spanning)*as.matrix(distxy)))
     }
