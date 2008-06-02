"mst.nb" <-
function(dxy){
    mymst=mstree(dxy)
    return(neig2nb(mymst))
}
