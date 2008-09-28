`heatload` <-
function(aspect, adjust=45, flat=0)
{
    if (!is.numeric(aspect)) {
        aspect <- tolower(substr(aspect,1,1))
        codes <-  c("n", "e", "s", "w", "f")
        if (sum(unique(aspect) %in% codes) != 5)
            stop("aspect as character misspecified")
        values <- c(0, 90, 180, 270, flat) + adjust
        tmp <- numeric(length(aspect))
        for (i in 1:5) tmp[aspect==codes[i]] <- values[i]
        aspect <- tmp
    } else {
        aspect <- aspect - adjust}
    hli <- (1 - cos(aspect * pi / 180)) / 2
return(hli)}

