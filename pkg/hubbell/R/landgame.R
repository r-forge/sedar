"landgame" <-
function(land, D, m, theta)
{
  tmp <- dim(land)
  J <- tmp[3]
  Nx <- tmp[1]
  Ny <- tmp[2]
  JM <- J*Nx*Ny
  nu <- theta/2/JM
                                        # Shuffle: prepare to kill first D
  for (i in 1:Nx)
    for (j in 1:Ny) 
      land[i,j,] <- land[i,j,sample(J)]
                                        # Recolonize first D slots
  for (i in 1:Nx)
    for (j in 1:Ny) {
                                        # Consider immigration
      imm <- rbinom(1, D, m)
      if (imm) {
        pool <- NULL
        for (step in c(-1,1)) {
          i.step <- i - step
          if (i.step < 1) i.step <- Nx
          if (i.step > Nx) i.step <- 1
          pool <- c(pool, land[i.step, j, (D+1):J])
          j.step <- j - step
          if (j.step < 1) j.step <- Ny
          if (j.step > Ny) j.step <- 1
          pool <- c(pool, land[i, j.step, (D+1):J])
        }
        land[i,j,1:imm] <- sample(pool, imm, replace=TRUE)
      }
                                        # Fill the rest from the local community
      if (imm < D)
        land[i,j, (imm+1):D] <- sample(land[i,j,(D+1):J], D-imm, replace=TRUE)
                                        # Mutate
      evolved <- runif(D) <= nu
      if (any(evolved))
        for (sp in which(evolved))
          land[i,j,sp] <- mutate(land[i,j,sp])
    }
  land
}
