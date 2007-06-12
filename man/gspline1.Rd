\name{gspline1}
\alias{gspline1}
\alias{dgspline1}
\alias{rgspline1}
\title{
  Density and random number generation from a univariate G-spline (penalized Gaussian mixture)
}
\description{
  Univariate G-spline (penalized Gaussian mixture) is distributed as
  \deqn{\alpha + \sum_{j=-K}^K w_j\mbox{N}(\tau\mu_j,\,\tau^2\sigma_j^2)}{%
  alpha + sum[j=-K][K]w[j]N(tau*mu[j], tau^2*sigma[j]^2).}  
}
\usage{
rgspline1(n, mu, sigma, weight, intcpt=0, scale=1, logw=TRUE)

dgspline1(x, mu, sigma, weight, intcpt=0, scale=1, logw=TRUE)
}
\arguments{
  \item{n}{number of observations to be generated}
  \item{x}{grid of values at which we evaluate the G-spline values}
  \item{mu}{a vector with G-spline knots
    \eqn{\mu_j\;(j=-K,\dots,K)}{mu[j] (j=-K,...,K)}(means of basis G-splines)}
  \item{sigma}{basis standard deviation(s)
    \eqn{\sigma_j\;(j=-K,\dots,K)}{sigma[j] (j=-K,...,K)}.
    If a single number is
    supplied then it is assumed that all basis G-splines have the same
    standard deviation. Alternatively a vector of the same length as
    \code{mu} can be given in which case the basis G-splines do not
    necessarily have the same standard deviations}
  \item{weight}{a vector with G-spline (log-)weights. It should have the
    same length as \code{mu}}  
  \item{intcpt}{G-spline intercept value \eqn{\alpha}{alpha})}
  \item{scale}{G-spline scale value \eqn{\tau}{tau})}
  \item{logw}{logical indicating whether logarithmic weights are
    supplied in \code{weight}}
}
\value{
  Values of the density or generated random numbers.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek@med.kuleuven.be}
}  
\examples{
knots <- c(-2, 0, 2)
sigma <- 1
weight <- c(0.3, 0.1, 0.6)
intcpt <- 3
scale <- 0.2
xgrid <- seq(1.8, 4.2, length=300)

dx <- dgspline1(xgrid, mu=knots, sigma=sigma, weight=weight,
   intcpt=intcpt, scale=scale, logw=FALSE)
x <- rgspline1(100, mu=knots, sigma=sigma, weight=weight,
   intcpt=intcpt, scale=scale, logw=FALSE)
hist(x, col="seagreen2", prob=TRUE, xlim=range(xgrid), xlab="x", ylab="g(x)")
lines(xgrid, dx, col="red", lwd=2)
}
\keyword{smooth}

