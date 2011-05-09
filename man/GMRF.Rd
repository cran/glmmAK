\name{GMRF}
\alias{GMRF}
\alias{momentsGMRF}
\alias{rGMRF}
\alias{dGMRF}
\alias{dGMRF2}
\title{
  Gaussian Markov random fields
}
\description{
  Moments, density and random generation for the Gaussian Markov random field
  with mean equal to 'mean', precision matrix equal to 'Q' (or covariance
  matrix equal to 'Sigma') and possibly constrained by a linear
  constraint 'Ax=b'.

  Generation of random numbers is performed by Algorithm 2.6 in Rue and
  Held (2005, pp. 38). 
}
\usage{
momentsGMRF(mean=0, Q=1, Sigma, A, b=0)

rGMRF(n, mean=0, Q=1, Sigma, A, b=0)

dGMRF(x, mean=0, Q=1, Sigma, A, b=0, log=FALSE)

dGMRF2(x, mean=0, Q=1, Sigma, A, b=0, log=FALSE)
}
\arguments{
  \item{mean}{vector of mean. If \code{length(mean)} is equal to 1, it
    is recycled and all components have the same mean.}
  \item{Q}{precision matrix of the GMRF.}
  \item{Sigma}{covariance matrix of the GMRF. Only one of \code{Q} and
    \code{Sigma} must be given. If \code{Sigma} is supplied, precision
    is computed from \eqn{\Sigma}{Sigma} as \eqn{Q = \Sigma^{-1}}{Q =
      Sigma^{-1}}.}
  \item{A}{optional matrix defining the constraint \eqn{Ax=b}{Ax=b} for
    sampled vectors \eqn{x}{x}.

    If not supplied, the GMRF is assumed to be unconstrained.

    Currently at most 1 constraint is allowed, that is \eqn{A}{A} must
    be a vector and \eqn{b}{b} a number.
   }
  \item{b}{vector or the right-hand side of the constraint. If
    \code{length(b)} is equal to 1, it is recycled and all constraint
    right-hand sides are the same.}
  \item{n}{number of observations to be sampled.}
  \item{x}{vector or matrix of the points where the density should be
    evaluated.}
  \item{log}{logical; if \code{TRUE}, log-density is computed}
}
\value{
Some objects.
}
\section{Value for momentsGMRF}{
  A list with the components:
  \describe{
    \item{mean}{mean of the (constrained) GMRF}
    \item{Sigma}{covariance matrix of the (constrained) GMRF}
  }
  
  and the following attributes:
  \describe{
    \item{mean.unconstr}{mean of the GMRF before imposing the constraints}
    \item{Sigma.unconstr}{covariance matrix \eqn{\Sigma_u}{Sigma[u]} of the GMRF before imposing
      the constraints}
    \item{Q.unconstr.cholesky}{Cholesky decomposition of the matrix \eqn{Q_u=\Sigma_u^{-1}}{Q[u]=Sigma[u]^(-1)}}
    \item{nconstraint}{number of constraints}
    \item{A}{left-hand side of the constrains}
    \item{b}{right-hand side of the constraints}    
  }  
}
\section{Value for rGMRF}{
  A list with the components:
  \describe{
    \item{x}{vector or matrix with sampled values}
    \item{log.dens}{vector with the values of the log-density evaluated
      in the sampled values}
  }
}
\section{Value for dGMRF, dGMRF2}{
  A vector with evaluated values of the (log-)density
}  
\references{
  Rue, H. and Held, L. (2005).
  \emph{Gaussian Markov Random Fields: Theory and Applications}.
  Boca Raton: Chapman and Hall/CRC.  
}
\seealso{
  \code{\link{dnorm}}, \code{\link[mvtnorm]{Mvnorm}}.
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
set.seed(1977)

mu <- c(0, 6, 8)
L <- matrix(1:9, nrow=3)
L[upper.tri(L, diag=FALSE)] <- 0
Sigma <- L \%*\% t(L)
Q <- chol2inv(chol(Sigma))

A <- rep(1, nrow(Sigma))
b <- 0

##### Unconstrained GMRF
##### ==================
## Moments
momentsGMRF(mean=mu, Sigma=Sigma)
momentsGMRF(mean=mu, Q=Q)

## Random numbers
z <- rGMRF(1000, mean=mu, Sigma=Sigma)
apply(z$x, 2, mean)
var(z$x)

## Random numbers, again
z <- rGMRF(10, mean=mu, Sigma=Sigma)
print(z)

## Values of the log-density
dGMRF(z$x, mean=mu, Sigma=Sigma, log=TRUE)
dGMRF(z$x, mean=mu, Q=Q, log=TRUE)
dGMRF2(z$x, mean=mu, Sigma=Sigma, log=TRUE)
dGMRF2(z$x, mean=mu, Q=Q, log=TRUE)

## Values of the density
dGMRF(z$x, mean=mu, Sigma=Sigma)
dGMRF(z$x, mean=mu, Q=Q)
dGMRF2(z$x, mean=mu, Sigma=Sigma)
dGMRF2(z$x, mean=mu, Q=Q)


##### Constrained GMRF
##### ================
## Moments
momentsGMRF(mean=mu, Sigma=Sigma, A=A, b=b)
momentsGMRF(mean=mu, Q=Q, A=A, b=b)

## Random numbers
z <- rGMRF(1000, mean=mu, Sigma=Sigma, A=A, b=b)
apply(z$x, 2, mean)
var(z$x)

## Random numbers, again
z <- rGMRF(10, mean=mu, Sigma=Sigma, A=A, b=b)
print(z)
A \%*\% t(z$x)

## Values of the log-density
dGMRF(z$x, mean=mu, Sigma=Sigma, A=A, b=b, log=TRUE)
dGMRF(z$x, mean=mu, Q=Q, A=A, b=b, log=TRUE)
dGMRF2(z$x, mean=mu, Sigma=Sigma, A=A, b=b, log=TRUE)
dGMRF2(z$x, mean=mu, Q=Q, A=A, b=b, log=TRUE)

## Values of the log-density
dGMRF(z$x, mean=mu, Sigma=Sigma, A=A, b=b)
dGMRF(z$x, mean=mu, Q=Q, A=A, b=b)
dGMRF2(z$x, mean=mu, Sigma=Sigma, A=A, b=b)
dGMRF2(z$x, mean=mu, Q=Q, A=A, b=b)
}
\keyword{distribution}
\keyword{multivariate}
