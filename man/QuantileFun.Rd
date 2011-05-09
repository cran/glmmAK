\name{QuantileFun}
\alias{QuantileFun}
\title{
  Sample quantiles
}
\description{
  This is (almost) the same as \code{quantile(x)}, or
  \code{apply(x, 1, quantile)} or \code{apply(x, 2, quantile)}.

  The motivation to write it was to validate my C++ function.
}
\usage{
QuantileFun(x, probs=seq(0, 1, 0.25), vals.in.cols=TRUE)
}
\arguments{
  \item{x}{values of the function}
  \item{probs}{numeric vector of probabilities with values in \eqn{[0,1]}{[0,1]}}
  \item{vals.in.cols}{if \code{TRUE} then it is assumed that function \eqn{f}{f} evaluated in a specific grid point
    over (MCMC) iterations is stored in a column of \code{x}. That is
    (MCMC) iterations correspond to rows.
  }  
}
\value{
  A \code{data.frame} with 1 row for each \code{probs} value.
}
\seealso{
  \code{\link[stats]{quantile}}.
}  
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
probs <- c(0, 0.25, 0.354, 0.5, 0.75, 1)

x <- rnorm(1001)
QuantileFun(x, probs=probs)
quantile(x, probs=probs)

n <- 1001
xx <- data.frame(x1=rnorm(n), x2=rgamma(n, shape=1, rate=1), x3=1:n)
QuantileFun(xx, probs=probs, vals.in.cols=TRUE)
apply(xx, 2, quantile, probs=probs)

xx2 <- t(xx)
QuantileFun(xx2, probs=probs, vals.in.cols=FALSE)
apply(xx2, 1, quantile, probs=probs)
}
\keyword{univar}
