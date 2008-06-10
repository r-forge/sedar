"makeland" <-
function(nrow, ncol=nrow, J)
{
  land <- array(rgb(1,1,1), dim=c(nrow,ncol,J))
  class(land) <- "landscape"
  land
}
