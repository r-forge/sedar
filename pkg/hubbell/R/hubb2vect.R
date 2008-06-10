"hubb2vect" <-
function(vec)
{
  S <- length(vec)
  J <- sum(vec)
  lim <- c(0,cumsum(vec))
  hubb.x <- vector("numeric",J)
  for (i in 1:S) hubb.x[(lim[i]+1):lim[i+1]] <- i
  class(hubb.x) <- "hubbell"
  hubb.x
}
