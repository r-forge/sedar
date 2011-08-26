\name{E.S}
\alias{E.S}
\alias{Var.S}
\alias{E.Somega}
\alias{wtheta.fun}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Functions for estimating Hubbell's theta from species --
  area curve }
\description{
  Functions for expected species number and its variance at certain
  sample size.  These functions can be used in non-linear regression to
  estimate Hubbell's ultimate diversity number \eqn{\theta}.
}
\usage{
E.S(theta, J)
Var.S(theta, J)
E.Somega(theta, omega, J)
wtheta.fun(S, J, theta)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{theta}{Hubbell's ultimate diversity number \eqn{\theta}. }
  \item{J}{Community size \eqn{J}. }
  \item{omega}{A descriptor for community isolation \eqn{\omega}.}
  \item{S}{Expected species richness \eqn{S}.}
}
\details{
  Hubbell predicts a specific form for species area curve, so that the
  expected number of species in a metacommunity is dependent only on
  number of individuals \eqn{J_M} and Hubbell's ultimate diversity
  number \eqn{\theta}.  Thus the species--area curve can be used to
  estimate that ultimate \eqn{\theta}. Further, Hubbell outlines analysis
  for local community.

  The expected species richness in Hubbell's model is approximately
  \eqn{S = 1 + \theta \cdot \ln (1 + \frac{J_M - 1}{\theta})}{S = 1 + \theta
    ln(1 + (J-1)/theta)}.  However, function \code{E.S} does not use this
  approximation but exact result from Hubbell's species generator (see
  \code{\link{rhubbell}}).  Function \code{Var.S} similarly gives
  the exact estimator of variance of expected species number from
  Hubbell's species generator.  Function \code{E.Somega} uses the
  species generator with isolation parameter \eqn{\omega \geq 0}{omega
    >= 0} for a local community; no variance estimator is available.
  
  Function \code{wtheta.fun} combines \code{E.S} and \code{Var.S} for
  variance weighted least squares fitting.  It returns the
  standardized residuals, or prediction errors divided by the standard
  error. 

  These functions together can be used for non-linear least squares
  regression with \code{\link{nls}}.

}
\value{
  Functions return the specified value each.
}
\references{
  Hubbell, S.P. (2001). The Unified Neutral Theory of Biodiversity and
  Biogeography.  Princeton Univ. Press.
}
\author{ Jari Oksanen }
\note{
  Hubbell gives only an approximate equation for expected species
  number.  The exact form and its variance were derived by J. Oksanen
  and may be unpublished.

  Hubbell's model for local community is based on wrong argumentation:
  Hubbell assumes falsely that his model is asymptotically similar to
  the classical Arrhenius model.  Even with that model, the
  argumentation in deriving \eqn{\omega} is invalid.  Moreover,
  \eqn{\omega} and \eqn{\theta} are so strongly correlated (positively)
  that their estimates are very unreliable.
}



\seealso{
  \code{\link{nls}} for non-linear regression.
}

\examples{
## Use Carabid data. First get the total count per pitfall (J) and species
## richness (S). 
##
data(carabid)
J <- apply(carabid, 1, sum)
S <- apply(carabid>0, 1, sum)
plot(J, S)
##
## Non-linear least squares
##
sparea <- nls(S ~ E.S(theta,J), start=list(theta=10))
sparea
i <- order(J)  ## Order on x axis for nice lines
lines(J[i], fitted(sparea)[i], col="blue", lwd=2)
##
## Then variance weighted least squares.
##
sparea.w <- nls(~ wtheta.fun(S, J, theta),start=list(theta=coef(sparea)))
sparea.w
## We used one sided model formula, and so the fitted values would be
## actually the weighted residuals.  Therefore we need to use E.S to
## plot fitted values.
lines(J[i], E.S(coef(sparea.w), J[i]), col="red", lwd=2)
##
## Now with dispersal limitation omega.
## This takes time: so increase iteration limit and ask to
## trace iteration, so that you believe someting is happening
##
sparea.om <- nls(S ~ E.Somega(theta,omega,J), start=list(theta=4, omega=0),
   control=nls.control(maxiter=500), trace=1)
sparea.om
lines(J[i], fitted(sparea.om)[i], lwd=2)
##
## Some traditional alternatives
##
arrhenius <- lm(log(S) ~ log(J))
lines(J[i], exp(fitted(arrhenius))[i])
# Doubling J brings along constant number of new species:
logadd <- lm(S ~ log(J))
lines(J[i], fitted(logadd)[i], lty=2)
}
\keyword{ models }% at least one, from doc/KEYWORDS

