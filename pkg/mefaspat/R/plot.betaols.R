`plot.betaols` <-
function(x, leg=TRUE, leg.pos="topright", leg.text=c("Group 1", "Group 2"), crit=0.05,
xlab, ylab, xlim, ylim, ...)
{
kmch <- ""
geo <- x$geogr
if (geo == "great circle") {
    kmch <- "(km)"
    geo <- "Great circle"}
diso <- if (attr(x, "dis.only") == TRUE) "Dissimilarity" else "Beta diversity"
if (missing(xlim)) xlim <- range(c(range(x$mat$gcd1), range(x$mat$gcd2)))
if (missing(ylim)) ylim <- range(c(range(x$mat$beta1), range(x$mat$beta2)))
if (missing(xlab)) xlab <- paste(geo, " distance ", kmch, sep="")
if (missing(ylab)) ylab <- paste(diso, ", index=\"", x$betadiv, "\"", sep="")

plot(0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
abline(lm(x$mat$beta1 ~ x$mat$gcd1), lty=2)
abline(lm(x$mat$beta2 ~ x$mat$gcd2), lty=3)
if (!is.null(crit)) for (i in 1:length(x$x.loc)) {
    xx <- x$x.loc[i]
    y1 <- x$coef.gr1[1] + xx * x$coef.gr1[2]
    y2 <- x$coef.gr2[1] + xx * x$coef.gr2[2]
    p.ty <- if (x$test[(i+1),2] < crit) 1 else 2
    arrows(xx,y1,xx,y2,lty=p.ty, code=0)
    }
points(x$mat$beta1 ~ x$mat$gcd1, pch=19)
points(x$mat$beta2 ~ x$mat$gcd2, pch=21)

if (leg) legend(leg.pos, pch=c(19,21), lty=c(2,3), legend=leg.text)

invisible(NULL)
}
