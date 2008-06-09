`nbmc` <-
function(loc, reg, hab, process=c("dl","ss","ls","me", "zd", "ud"), cont=rep(1,length(loc))){
    # internal
    metac <- function(i, hab, reg, process, cont)
        {
        ## dispersal limitation
        if(process=="dl") out <- setdiff(loc[reg == reg[i]], loc[i])
        ## species sorting
        if(process=="ss") out <- setdiff(loc[hab == hab[i]], loc[i])
        ## limitation & sorting
        if(process=="ls") out <- setdiff(intersect(loc[reg == reg[i]], loc[hab == hab[i]]), loc[i])
        #mass effect
        if(process=="me") out <- setdiff(union(loc[reg == reg[i]], loc[hab == hab[i]]), loc[i])
        #zero dispersal
        if(process=="zd") out <- numeric(0)
        #unlimited dispersal
        if(process=="ud") out <-setdiff(loc, loc[i])
        out <- intersect(out, loc[cont == cont[i]])
        return(out)}
    res <- lapply(1:length(loc), function(x) metac(x, reg=reg, hab=hab, process=process, cont=cont))
    names(res) <- loc
    return(res)
    }

