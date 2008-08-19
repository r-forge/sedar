`print.betaols` <-
function(x, decimals=4, ...)
{
locs <- if (nrow(x$test) == 2) "Location" else "Locations"
cat("Object of class 'betaols'\n")
cat("Number of permutations:",x$n.perm,"\n")
cat("Index of beta diversities: \"", x$betadiv, "\"\nGeographic distances: ",x$geogr,"\n", sep="")
cat(locs, "for comparison:",x$x.loc,"\n")
cat("Total number of species:",ncol(x$m),"\n")
cat("Number of species in group 1:",length(x$gr1),"\n")
cat("Number of species in group 2:",length(x$gr2),"\n")
}
