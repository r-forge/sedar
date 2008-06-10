"wtheta.fun" <-
function (S,J,theta) 
{
  pred <- E.S(theta,J)
  var.pred <- Var.S(theta,J)
  (S-pred)/sqrt(var.pred)
}
