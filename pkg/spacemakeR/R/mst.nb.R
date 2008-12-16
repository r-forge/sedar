"mst.nb" <-
function(dxy){
    mymst <- ade4::mstree(dxy)
    return(neig2nb(mymst))
}
