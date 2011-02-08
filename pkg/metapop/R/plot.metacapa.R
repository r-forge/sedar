"plot.metacapa" <-
    function(x, locations = NULL, cex = 4, ...)
{
    if (is.null(locations))
        locations <- cmdscale(x$d)
    opal <- palette()
    on.exit(opal)
    palette(heat.colors(255))
    cl <-  (1 - x$x/max(x$x))*254 + 1
    ## Make largest dot cex=5
    dia <- sqrt(max(x$A))/cex
    plot(locations, asp=1, cex=sqrt(x$A)/dia, xlab="", ylab="", pch=21,
         col="blue", bg=cl, ...)
    invisible(locations)
}
