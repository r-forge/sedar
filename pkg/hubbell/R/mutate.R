"mutate" <-
function(genes)
{
  tmp <- col2rgb(genes)/255
  i <- sample(3,1)
  tmp[i,] <- runif(1)
  genes <- rgb(tmp[1,],tmp[2,],tmp[3,])
  genes
}
