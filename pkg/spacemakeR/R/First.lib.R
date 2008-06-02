".First.lib" <- function(lib, pkg) {
  cat("\nspacemakeR\n")
  cat("\nDray S., Legendre P. and Peres-Neto P. R. (2006)")
  cat("\nSpatial modeling: a comprehensive framework for principal")
  cat("\ncoordinate analysis of neighbor matrices (PCNM). Ecological Modelling, 3-4: 483-493\n\n")
  cat('Read the tutorial available (vignette in /inst/doc directory) using vignette("tutorial",package="spacemakeR").\n\n')
  require(ade4)
  require(spdep)
  require(tripack)
  cat("\n\n")
}
