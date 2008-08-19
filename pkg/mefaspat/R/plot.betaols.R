`plot.betaols` <-
function(x, leg=TRUE, leg.pos="topright", leg.text=c("Group 1", "Group 2"), crit=0.05, ...)
{
kmch <- ""
geo <- x$geogr
if (geo == "great circle") {
    kmch <- "(km)"
    geo <- "Great circle"}
diso <- if (attr(x, "dis.only") == TRUE) "Dissimilarity" else "Beta diversity"
xlim <- range(c(range(x$mat$gcd1), range(x$mat$gcd2)))
xlab <- paste(geo, " distances ", kmch, sep="")
ylab <- paste(diso, ", index=\"", x$betadiv, "\"", sep="")

plot(0, type="n", xlim=xlim, ylim=c(0,1), xlab=xlab, ylab=ylab, ...)
abline(lm(x$mat$beta1 ~ x$mat$gcd1), lty=2)
abline(lm(x$mat$beta2 ~ x$mat$gcd2), lty=3)
if (!is.null(crit)) for (i in 1:length(x$x.loc)) {
    p.ty <- if (x$test[(i+1),2] < crit) 1 else 2
    lines(c(x$x.loc[i],x$x.loc[i]),c(-1,2),lty=p.ty)
    }
points(x$mat$beta1 ~ x$mat$gcd1, pch=19)
points(x$mat$beta2 ~ x$mat$gcd2, pch=21)

if (leg) legend(leg.pos, pch=c(19,21), lty=c(2,3), legend=leg.text)

invisible(NULL)
}
