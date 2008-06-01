`spantree2neig` <-
    function(tree)
{
    n <- length(tree$kid) + 1
    neig(edges = cbind(2:n, tree$kid))
}

