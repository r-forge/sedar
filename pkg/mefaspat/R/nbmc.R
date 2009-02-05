`nbmc` <-
function(loc, reg, hab=NULL, process=c("dl","ss","ls","me", "zd", "ud"), cont=rep(1, length(loc)))
{
    process <- match.arg(process)
    if (is.null(hab))
        hab <- rep(1,length(loc))
    ## internal function
    metac <- function(i, hab, reg, process, cont)
    {
        out <- switch(process,
            ## dispersal limitation
            "dl" = setdiff(loc[reg == reg[i]], loc[i]),
            ## species sorting
            "ss" = setdiff(loc[hab == hab[i]], loc[i]),
            ## limitation & sorting
            "ls" = setdiff(intersect(loc[reg == reg[i]], loc[hab == hab[i]]), loc[i]),
            ## mass effect
            "me" = setdiff(union(loc[reg == reg[i]], loc[hab == hab[i]]), loc[i]),
            ## zero dispersal
            "zd" = numeric(0),
            ## unlimited dispersal
            "ud" = setdiff(loc, loc[i])
        )
        intersect(out, loc[cont == cont[i]])
    }
    res <- lapply(1:length(loc), function(x) metac(x, reg=reg, hab=hab, process=process, cont=cont))
    names(res) <- loc
    return(res)
    }

