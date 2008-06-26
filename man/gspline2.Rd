\name{gspline2}
\alias{gspline2}
%\alias{dgspline2}
\alias{rgspline2}
\title{
  Density and random number generation from a bivariate G-spline (penalized Gaussian mixture)
}
\description{
  Density has not been implemented yet.

  Bivariate G-spline (penalized Gaussian mixture) is distributed as
  \deqn{(\alpha_1,\,\alpha_2)' + \sum_{j_1=-K_1}^{K_1}\sum_{j_2=-K_2}^{K_2} w_{j_1,j_2}\mbox{N}_2((\tau_1\mu_{1,j_1},\,\tau_2\mu_{2,j_2})',\,\mbox{diag}(\tau_1^2\sigma_{1,j_1}^2,\;\tau_2^2\sigma_{2,j_2}^2))}{%
  (alpha[1], alpha[2])' +
  sum[j1=-K1][K1]sum[j2=-K2][K2]w[j1,j2]N((tau[1]*mu[1,j1],
  tau[2]*mu[2,j2])', diag(tau[1]^2*sigma[1,j1]^2, tau[2]^2*sigma[2,j2]^2)).}  
}
\usage{
rgspline2(n, mu1, mu2, sigma1, sigma2, weight, knotInd,
   intcpt=0, scale=1, logw=TRUE)
}
\arguments{
  \item{n}{number of observations to be generated}  
  \item{mu1}{a vector with G-spline knots
    \eqn{\mu_{1,j_1}\;(j_1=-K_1,\dots,K_1)}{mu[1,j1] (j1=-K1,...,K1)}
    (means of basis G-splines) in the 1st margin}
  \item{mu2}{a vector with G-spline knots
    \eqn{\mu_{2,j_2}\;(j_2=-K_2,\dots,K_2)}{mu[2,j2] (j2=-K2,...,K2)}
    (means of basis G-splines) in the 2nd margin}
  \item{sigma1}{basis standard deviation(s)
    \eqn{\sigma_{1,j_1}\;(j_1=-K_1,\dots,K_1)}{sigma[1,j1] (j1=-K1,...,K1)}
    in the 1st margin.

    If a single number is
    supplied then it is assumed that all basis G-splines have the same
    standard deviation. Alternatively a vector of the same length as
    \code{mu1} can be given in which case the basis G-splines do not
    necessarily have the same standard deviations}
  \item{sigma2}{basis standard deviation(s)
    \eqn{\sigma_{2,j_2}\;(j_2=-K_2,\dots,K_2)}{sigma[2,j2] (j2=-K2,...,K2)}
    in the 2nd margin.

    If a single number is
    supplied then it is assumed that all basis G-splines have the same
    standard deviation. Alternatively a vector of the same length as
    \code{mu2} can be given in which case the basis G-splines do not
    necessarily have the same standard deviations}  
  \item{weight}{a vector or matrix with G-spline (log-)weights.

    If \code{missing(knotInd)} then it should have the
    same length as \code{length(mu1)*length(mu2)}.

    If \code{knotInd} is supplied then \code{weight} must be a vector
    of the same length as \code{knotInd} which contains (log-)weights
    of G-spline components having non-zero weights.
  }
  \item{knotInd}{If supplied then it contains indeces of G-spline
    components which correspond to non-zero weights. Indeces should be
    on the scale from 0 to \code{length(mu1)*length(mu2)-1} (similarly
    like in the file knotInd.sim created, e.g.,  by the function \code{\link{cumlogitRE}})}
  \item{intcpt}{G-spline intercept value(s) \eqn{(\alpha_1,\,\alpha_2)'}{(alpha[1],alpha[2])'}

    If a single value is supplied then it is assumed that intercept
    values in both margins are the same.
  }
  \item{scale}{G-spline scale value(s) \eqn{(\tau_1,\,\tau_2)'}{(tau[1],tau[2])'}

    If a single value is supplied then it is assumed that scale
    values in both margins are the same.
  }
  \item{logw}{logical indicating whether logarithmic weights are
    supplied in \code{weight}}
}
\value{
  Values of the density or generated random numbers.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
knots1 <- c(-2, 0, 2)
knots2 <- c(-2, -1, 0, 1, 2)
sigma1 <- 0.5
sigma2 <- 0.5
intcpt <- c(3, -1)
scale <- c(0.2, 0.5)

weight <- matrix(c(1,2,1,  2,3,2,  3,4,3,  2,3,2,  1,2,1), ncol=5)
xA <- rgspline2(1000, mu1=knots1, mu2=knots2, sigma1=sigma1,
   sigma2=sigma2, weight=weight, intcpt=intcpt, scale=scale, logw=FALSE)
xB <- rgspline2(1000, mu1=knots1, mu2=knots2, sigma1=sigma1,
   sigma2=sigma2, weight=log(weight), intcpt=intcpt, scale=scale, logw=TRUE)

oldpar <- par(mfrow=c(1, 2), bty="n")
plot(xA[,2], xA[,1], pch=16, col="red")
abline(h=intcpt[1]+scale[1]*knots1, col="orange")
abline(v=intcpt[2]+scale[2]*knots2, col="orange")
plot(xB[,2], xB[,1], pch=16, col="blue")
abline(h=intcpt[1]+scale[1]*knots1, col="orange")
abline(v=intcpt[2]+scale[2]*knots2, col="orange")
par(oldpar)

### Only selected components with non-zero weights
logweight2 <- c(-1, -0.5, 1)
knotInd <- c(0, 5, 14)
yA <- rgspline2(1000, mu1=knots1, mu2=knots2, sigma1=sigma1,
   sigma2=sigma2, weight=logweight2, knotInd=knotInd,
   intcpt=intcpt, scale=scale, logw=TRUE)
yB <- rgspline2(1000, mu1=knots1, mu2=knots2, sigma1=sigma1,
   sigma2=sigma2, weight=exp(logweight2), knotInd=knotInd,
   intcpt=intcpt, scale=scale, logw=FALSE)
oldpar <- par(mfrow=c(1, 2), bty="n")
plot(yA[,2], yA[,1], pch=16, col="red")
abline(h=intcpt[1]+scale[1]*knots1, col="orange")
abline(v=intcpt[2]+scale[2]*knots2, col="orange")
plot(yB[,2], yB[,1], pch=16, col="blue")
abline(h=intcpt[1]+scale[1]*knots1, col="orange")
abline(v=intcpt[2]+scale[2]*knots2, col="orange")
par(oldpar)

logweight3 <- log(weight)[-c(5,8,11)]
knotInd3 <- c(0,1,2, 3,5, 6,8, 9,11, 12,13,14)
zA <- rgspline2(1000, mu1=knots1, mu2=knots2, sigma1=sigma1,
   sigma2=sigma2, weight=logweight3, knotInd=knotInd3,
   intcpt=intcpt, scale=scale, logw=TRUE)
zB <- rgspline2(1000, mu1=knots1, mu2=knots2, sigma1=sigma1,
   sigma2=sigma2, weight=exp(logweight3), knotInd=knotInd3,
   intcpt=intcpt, scale=scale, logw=FALSE)
oldpar <- par(mfrow=c(1, 2), bty="n")
plot(zA[,2], zA[,1], pch=16, col="red")
abline(h=intcpt[1]+scale[1]*knots1, col="orange")
abline(v=intcpt[2]+scale[2]*knots2, col="orange")
plot(zB[,2], zB[,1], pch=16, col="blue")
abline(h=intcpt[1]+scale[1]*knots1, col="orange")
abline(v=intcpt[2]+scale[2]*knots2, col="orange")
par(oldpar)
}
\keyword{smooth}
