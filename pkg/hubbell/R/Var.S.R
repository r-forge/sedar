"Var.S" <-
function(theta, J) { E.S(theta, J) + theta^2 * (trigamma(theta+J) - trigamma(theta)) }
