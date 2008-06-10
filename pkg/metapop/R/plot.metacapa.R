"plot.metacapa" <-
    function(x, locations = NULL,  ...)
{
    if (is.null(locations))
        locations <- cmdscale(x$d)
    opal <- palette()
    on.exit(opal)
    palette(heat.colors(255))
    cl <-  (1 - x$x)*254 + 1
    plot(locations, asp=1, cex=sqrt(x$A*10), xlab="", ylab="", pch=21,
         col="blue", bg=cl, ...)
    invisible(locations)
}
