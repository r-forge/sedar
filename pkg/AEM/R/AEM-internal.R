`.First.lib` <-
function(lib, pkg) {
  library.dynam("AEM", pkg, lib)
  require(spdep)
}

