"plot.metacapa" <-
    function(x, locations = NULL,  ...)
{
    if (is.null(locations))
        locations <- cmdscale(x$d)
    opal <- palette()
    on.exit(opal)
    palette(heat.colors(255))
    cl <-  (1 - x$x/max(x$x))*254 + 1
    ## Find adjustment for size: largest circle is 1/25 of plot range
    ran <- max(diff(apply(locations, 2, range)))
    dia <- sqrt(max(x$A))
    adj <- ran/dia/25
    plot(locations, asp=1, cex=sqrt(x$A*adj), xlab="", ylab="", pch=21,
         col="blue", bg=cl, ...)
    invisible(locations)
}
