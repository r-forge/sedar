"give.thresh" <-
function(distxy){
    spanning=mstree(distxy)
    return(max(neig2mat(spanning)*as.matrix(distxy)))
     }
