"plot.hubbell" <-
function (x, sp.max = 64, ...) 
{
    op <- par(no.readonly=TRUE)
    on.exit(par(op))
    opal <- palette()
    on.exit(palette(opal), add = TRUE)
    palette(rainbow(sp.max))
    layout(matrix(c(1, 1, 2, 1, 1, 2), 2, 3, byrow = TRUE))
    dim <- ceiling(sqrt(length(x)))
    xy <- expand.grid(1:dim, 1:dim)[1:length(x), ]
    plot(xy, xlab = "", ylab = "", axes = FALSE, pch = 16, col = x, 
         ...)
    box()
    rad <- rev(sort(table(x)))
    plot(1:length(rad), rad, log = "y", type = "l", xlab = "Rank", 
         ylab = "Abundance", xlim = c(1, sp.max),
         ylim = c(1, length(x)), ...)
    points(1:length(rad), rad, col = names(rad), pch=16, ...)
    par(ylog=op$ylog)
    invisible()
}
