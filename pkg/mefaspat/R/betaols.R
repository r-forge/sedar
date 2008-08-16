`betaols` <-
function(data, latlong, gr1, gr2 = NULL, x.loc = NULL, method = "sim", gcd = TRUE, times = 100)
{
require(vegan)
x <- data
fields.ok <- if (gcd) require(fields, quietly=TRUE) else FALSE
if (!fields.ok & gcd) warning("euclidean distances were used instead of great circle distances")

if (is.null(gr2)) gr2 <- setdiff(c(1:ncol(x)), gr1)
if (length(gr1) >= ncol(x)) stop("length(gr1) should be less than ncol(x)")
if (length(c(gr1, gr2)) > ncol(x)) stop("length(c(gr1, gr2)) should be less than or equal to ncol(x)")
if (sum(intersect(gr1, gr2)) > 0) stop("gr1 and gr2 overlap")

if (is.null(x.loc)) x.loc <- 0 else x.loc <- c(0, x.loc)

## internal function START
subset.internal <-
function(x, latlong, gr1, gr2, method)
{
part1 <- x[,gr1]
part2 <- x[,gr2]
xy1 <- subset(latlong, apply(part1 > 0, 1, sum) != 0)
xy2 <- subset(latlong, apply(part2 > 0, 1, sum) != 0)
part1 <- subset(part1, apply(part1 > 0, 1, sum) != 0)
part2 <- subset(part2, apply(part2 > 0, 1, sum) != 0)

if (fields.ok & gcd) {
    ## great circle distances, note: km and not miles!
    gcd1 <- as.dist(rdist.earth(xy1[,2:1], xy1[,2:1], miles = FALSE))
    gcd2 <- as.dist(rdist.earth(xy2[,2:1], xy2[,2:1], miles = FALSE))
        } else {
    gcd1 <- dist(xy1, method="euclidean")
    gcd1 <- dist(xy2, method="euclidean")
    }

gcd1 <- dist(xy1)
gcd2 <- dist(xy2)

## beta_ distances as in betadiver()
beta1 <- betadiver(part1, index=method)
beta2 <- betadiver(part2, index=method)
if (is.element(method, c("j", "sor", "rlb"))) {
    beta1 <- 1 - beta1
    beta2 <- 1 - beta2
    }

return(list(gcd1=gcd1, gcd2=gcd2, beta1=beta1, beta2=beta2))
} ## internal function END

m.orig <- subset.internal(x=x, latlong=latlong, gr1=gr1, gr2=gr2, method=method)

c1 <- coef(lm(m.orig$beta1 ~ m.orig$gcd1))
c2 <- coef(lm(m.orig$beta2 ~ m.orig$gcd2))
names(c1)[2] <- "slope.gr1"
names(c2)[2] <- "slope.gr2"

tmp <- x.loc
for (i in 1:length(x.loc)) {
    tmp[i] <- ((c1[1] + x.loc[i] * c1[2]) - (c2[1] + x.loc[i] * c2[2]))
    names(tmp)[i] <- paste("x.loc.", x.loc[i], sep="")
    }

difference <- abs(c((c1[2] - c2[2]), tmp))
 names(difference)[1] <- "slope"

d.perm <- array(NA, times * length(difference))
start <- 1

for (i in 1:times){ ## for times START

gr12 <- sample(c(1:ncol(x)), (length(gr1) + length(gr2)))
gr1p <- sample(gr12, length(gr1))
gr2p <- setdiff(gr12, gr1p)

m.perm <- subset.internal(x=x, latlong=latlong, gr1=gr1p, gr2=gr2p, method=method)

c1p <- coef(lm(m.perm$beta1 ~ m.perm$gcd1))
c2p <- coef(lm(m.perm$beta2 ~ m.perm$gcd2))

for (i in 1:length(x.loc)) {
    tmp[i] <- ((c1p[1] + x.loc[i] * c1p[2]) - (c2p[1] + x.loc[i] * c2p[2]))
    }

d.perm[start:(start + length(x.loc))] <- c((c1p[2] - c2p[2]), tmp)
start <- start + length(x.loc) + 1
} ## for times END

d.perm <- t(matrix(d.perm, times, byrow=TRUE))
rownames(d.perm) <- names(difference)

p.value <- apply(abs(d.perm) >= difference, 1, function(x) sum(x)/times)

final <- cbind(difference, p.value)

distance <- if (fields.ok & gcd) "great circle" else "Euclidean"

out <- list(
    call=match.call(),
    data=data,
    gr1=gr1,
    gr2=gr2,
    n.perm=times,
    m=m.orig,
    d.perm=data.frame(t(d.perm)),
    betadiv=method,
    geogr=distance,
    coef.gr1=c1,
    coef.gr2=c2,
    test=final
    )
class(out) <- "betaols"
return(out)
} ## function END

`plot.betaols` <-
function(x, leg=TRUE, leg.pos="topright", leg.text=c("Group 1", "Group 2"), ...)
{

kmch <- ""
if (x$geogr == "great circle") kmch <- "(km)"

xlim <- range(c(range(x$m$gcd1), range(x$m$gcd2)))
xlab <- paste(x$geogr, " distances ", kmch, sep="")
ylab <- paste("dissimilarity, index=\"", x$betadiv, "\"", sep="")


plot(x$m$beta1 ~ x$m$gcd1, pch=19, xlim=xlim, ylim=c(0,1), xlab=xlab, ylab=ylab, ...)
abline(lm(x$m$beta1 ~ x$m$gcd1), lty=2)
points(x$m$beta2 ~ x$m$gcd2, pch=21)
abline(lm(x$m$beta2 ~ x$m$gcd2), lty=3)

if (leg) legend(leg.pos, pch=c(19,21), lty=c(2,3), legend=leg.text)

invisible(NULL)
}

`print.betaols` <-
function(x, decimals=4, ...)
{

interc <- if (nrow(x$test) == 2) "intercept" else "intercepts"
cat("Permutation test for differences in slope and", interc, "of\n")
cat("OLS regression with",x$n.perm,"permutations, based on matrices of\n")
cat("\"", x$betadiv, "\" beta diversity index and ",x$geogr," distances.\n", sep="")
# cat("\nCall:")
# print(x$call)
# cat("\nCoefficients for species group 1:\n")
# print(round(x$coef.gr1,decimals))
# cat("\nCoefficients for species group 2:\n")
# print(round(x$coef.gr2,decimals))
# cat("\nPermutation test results:\n")
cat("\n")
print(round(x$test,decimals))
}

