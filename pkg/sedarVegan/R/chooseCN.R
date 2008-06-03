`chooseCN` <-
    function(xy, type = "delauney",  result.type="nb", d1 = NULL, d2 = NULL,
             k= NULL, a = NULL, dmin = NULL, support = c("vegan", "ade4"))
{
    support <- match.arg(support)
    if(is.data.frame(xy)) xy <- as.matrix(xy)
    if(ncol(xy) != 2) stop("xy does not have two columns.")
    if(any(is.na(xy))) stop("NA entries in xy.")
    result.type <- tolower(result.type)
    
    require(spdep) || stop("spdep library is required.")
    
    res <- list()
    d1.first <- d1
    d2.first <- d2
    k.first <- k
        
    ## check for uniqueness of coordinates
    if(any(duplicated(xy) && type != "inv.dist")){
        ## coords need not be unique if type==7 (inverse distances)
        xy <- jitter(xy)
        warning("Random noise was added to xy as duplicated coordinates existed.")
    }
        
    ## handle type argument
    type <- match.arg(type, c("delauney", "gabriel", "rel.neighbours", "mst",
                              "neighb.by.dist", "knearest", "inv.dist"))
    ## re-initialisation of some variables
    d1 <- d1.first
    d2 <- d2.first
    k <- k.first    
    ## graph types
    switch(type, delauney = {
        require(tripack, quiet=TRUE) || stop("tripack library is required.")
        cn <- tri2nb(xy)
    }, gabriel = {
        cn <- gabrielneigh(xy)
        cn <- graph2nb(cn, sym=TRUE)
    }, rel.neighbours = {
        cn <- relativeneigh(xy)
        cn <- graph2nb(cn, sym=TRUE)
    }, mst = {
        switch(support, ade4 = {
            require(ade4) ||  stop("ade4 library is required")
            cn <- mstree(dist(xy))
        }, vegan = {
            require(vegan) || stop("vegan library is required")
            cn <- spantree2neig(spantree(dist(xy)))
        })
        cn <- neig2nb(cn)
    }, neighb.by.dist = {
        if(is.null(d1) || is.null(d2)) {
            stop("args d1 and d2 must be given with type =", type)
        }
        ## avoid that a point is its neighbour
        dmin <- mean(dist(xy))/100000
        if(d1<dmin) d1 <- dmin
        if(d2<d1) stop("d2 < d1")
        cn <- dnearneigh(x=xy, d1=d1, d2=d2)
    }, knearest = {
        if(is.null(k)) {
            stop("arg k must be given with type = ", type)
        }
        cn <- knearneigh(x=xy, k=k)
        cn <- knn2nb(cn, sym=TRUE)
    }, inv.dist = {
        if(is.null(a)) {
            stop("arg a must be given with type = ", type)
        }
        cn <- as.matrix(dist(xy))
        if(is.null(dmin)) {
            stop("arg dmin must be given with type =", type)
        }
        if(a < 1)  a <- 1 
        thres <- mean(cn)/1e8
        if(dmin > thres) dmin <- thres
        cn[cn < dmin] <- dmin
        cn <- 1/(cn^a)
        diag(cn) <- 0
        cn <- prop.table(cn,1)
        plot.nb <- FALSE
        edit.nb <- FALSE
        result.type <- "listw"
    })
    ## end graph types
        
    if(result.type == "listw") {
        if(type != "inv.dist") {
            cn <- nb2listw(cn, style="W", zero.policy = TRUE)
        } else {
            cn <- mat2listw(cn)
            cn$style <- "W"
        }
    }
    res <- cn
    attr(res,"xy") <- xy
    class(res) <- c("chooseCN", class(res))
    res
}
