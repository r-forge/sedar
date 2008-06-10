"rhubbell" <-
function(theta, J) 
{
  community <- NULL
  for (j in 0:(J-1)) {
    if (runif(1) < theta/(theta+j))
      community <- c(community, 1)
    else {
      species <- sample(length(community), 1, prob=community/j)
      community[species] <- community[species]+1
    }
  }
  return(community)   
}
