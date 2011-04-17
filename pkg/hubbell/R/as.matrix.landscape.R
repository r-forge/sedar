### Changes 'landgame' result to a community data frame where each
### site is a row, and each distinct species is a column
as.matrix.landscape <-
    function(x, ...)
{
    spec <- unique(x)
    nr <- prod(dim(x)[1:2])
    dim(x) <- c(nr, dim(x)[3])
    df <- matrix(0, nr, length(spec))
    colnames(df) <- spec
    for(i in 1:nrow(x)) {
        tbl <- table(x[i,])
        df[i, names(tbl)] <- tbl
    }
    df
}
