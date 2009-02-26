`betaols2` <-
function(datamatrix, latlong, gr1, gr2 = NULL, x.loc = NULL, method = "sim", gcd = TRUE, times = 100, dis.only=TRUE)
{
    require(vegan)
    x <- datamatrix
    fields.ok <- if (gcd)
        require(fields, quietly=TRUE) else FALSE
    if (!fields.ok & gcd)
        warning("euclidean distances were used instead of great circle distances")
    if (is.null(gr2))
        gr2 <- setdiff(c(1:ncol(x)), gr1)
    if (length(gr1) >= ncol(x))
        stop("length(gr1) should be less than ncol(x)")
    if (length(c(gr1, gr2)) > ncol(x))
        stop("length(c(gr1, gr2)) should be less than or equal to ncol(x)")
    if (sum(intersect(gr1, gr2)) > 0)
        stop("gr1 and gr2 overlap")
    if (is.null(x.loc))
        x.loc <- 0 else x.loc <- unique(c(0, x.loc))

## internal function START
subset.internal <- function(x, latlong, gr1, gr2, method) {
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
    if (dis.only) {
        if (is.element(method, c("j", "sor", "rlb"))) {
            beta1 <- 1 - beta1
            beta2 <- 1 - beta2
        }
    }
    return(list(gcd1=gcd1, gcd2=gcd2, beta1=beta1, beta2=beta2))
} ## internal function END

    m.orig <- subset.internal(x=x, latlong=latlong, gr1=gr1, gr2=gr2, method=method)

#    c1 <- coef(lm(m.orig$beta1 ~ m.orig$gcd1))
#    c2 <- coef(lm(m.orig$beta2 ~ m.orig$gcd2))

    ddd <- data.frame(beta1 = array(m.orig$beta1), beta2 = array(m.orig$beta2),
        gcd1 = array(m.orig$gcd1), gcd2 = array(m.orig$gcd2))
    c1 <- lm.fit(x = model.matrix(attr(model.frame(~ddd$gcd1), "terms")), y = ddd$beta1)$coefficients
    c2 <- lm.fit(x = model.matrix(attr(model.frame(~ddd$gcd2), "terms")), y = ddd$beta2)$coefficients

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

#        c1p <- coef(lm(m.perm$beta1 ~ m.perm$gcd1))
#        c2p <- coef(lm(m.perm$beta2 ~ m.perm$gcd2))

        ddp <- data.frame(beta1 = array(m.perm$beta1), beta2 = array(m.perm$beta2),
            gcd1 = array(m.perm$gcd1), gcd2 = array(m.perm$gcd2))
        c1p <- lm.fit(x = model.matrix(attr(model.frame(~ddp$gcd1), "terms")), y = ddp$beta1)$coefficients
        c2p <- lm.fit(x = model.matrix(attr(model.frame(~ddp$gcd2), "terms")), y = ddp$beta2)$coefficients

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
        m=datamatrix,
        gr1=gr1,
        gr2=gr2,
        x.loc=x.loc,
        n.perm=times,
        mat=m.orig,
        d.perm=data.frame(t(d.perm)),
        betadiv=method,
        geogr=distance,
        coef.gr1=c1,
        coef.gr2=c2,
        test=final)
    attr(out, "dis.only") <- dis.only
    class(out) <- c("betaols", "list")
    return(out)
} ## function END

