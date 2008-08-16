`summary.betaols` <-
function(object, decimals=4, ...)
{
x<-object
interc <- if (nrow(x$test) == 2) "intercept" else "intercepts"
cat("Permutation test for differences in slope and", interc, "of\n")
cat("OLS regression with",x$n.perm,"permutations, based on matrices of\n")
cat("\"", x$betadiv, "\" beta diversity index and ",x$geogr," distances.\n", sep="")
 cat("\nCall:")
 print(x$call)
 cat("\nCoefficients for species group 1:\n")
 print(round(x$coef.gr1,decimals))
 cat("\nCoefficients for species group 2:\n")
 print(round(x$coef.gr2,decimals))
 cat("\nPermutation test results:\n")
cat("\n")
print(round(x$test,decimals))
}

