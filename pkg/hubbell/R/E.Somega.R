"E.Somega" <-
function (theta, omega, J) 
{
  sol <- J
  iJ <- order(J)
  Jlist <- c(1,J[iJ])
  tmp <- 0
  for (i in 1:length(sol)) {
    for (j in Jlist[i]:Jlist[i+1]) {
      tmp <- tmp + theta/(theta+j-1)*j^(-omega)
    }
    sol[iJ[i]] <- tmp
  }
  sol
}
