\name{maxPosterProb}
\alias{maxPosterProb}
\title{
  G-spline utility
}
\description{
  For given G-spline basis and given data, it determines for each data
  point the G-spline component for which the value of the basis density
  is maximal.
}
\usage{
maxPosterProb(data, intercept, std.dev, K, delta, sigma)
}
\arguments{
  \item{data}{numeric vector or matrix with data. If given as a matrix
    then rows correspond to observations and columns to margins.
  }
  \item{intercept}{numeric vector of length 1 or \code{ncol(data)} with
    intercepts for each margin.

    If given as a number, it is recycled.
  }
  \item{std.dev}{numeric vector of length 1 or \code{ncol(data)} with
    standard deviations for each margin.

    If given as a number, it is recycled.
  }
  \item{K}{numeric vector of length 1 or \code{ncol(data)} which
    specifies, for each marginal G-spline, then number of knots
    on each side of the zero knot. That is, the
    \eqn{i}{i}-th marginal G-spline has \eqn{2K_i+1}{2K[i]+1} knots.

    If given as a number, it is recycled.  
  }
  \item{delta}{numeric vector of length 1 or \code{ncol(data)}
    which specifies the distance between two consecutive
    knots for each marginal G-spline. That is, the \eqn{i}{i}-th
    marginal G-spline has the following knots
    \deqn{\mu_{i,j} = j\,\delta_i,\; j=-K_i,\dots,K_i.}{%
      mu[i,j] = j*delta[i], j=-K[i],...,K[i].}

    If given as a number, it is recycled.
  }
  \item{sigma}{numeric vector of length 1 or \code{ncol(data)}
    with basis standard deviations for marginal G-splines.
    
    If given as a number, it is recycled.
  }  
}
\value{
  Matrix which specifies determined components (indeces are on scale \eqn{-K_i,\dots,K_i}{-K[i],...,K[i]}).
}
\seealso{
  \code{cumlogitRE}, \code{logpoissonRE}.
}  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
  N <- 100
  intcpt <- c(0, 5, 15)
  std.dev <- c(1, 0.5, 3)
  data <- data.frame(b1=rnorm(N, intcpt[1], std.dev[1]),
                     b2=rnorm(N, intcpt[2], std.dev[2]),
                     b3=rnorm(N, intcpt[3], std.dev[3]))
  alloc <- maxPosterProb(data=data, intercept=intcpt, std.dev=std.dev,
                         K=15, delta=0.3, sigma=0.2)

  par(mfrow=c(1, 3), bty="n")
  for (i in 1:3) hist(alloc[,i], prob=TRUE, col="seagreen3",
                      xlab="Allocation", breaks=(-15):15,
                      main=paste("Margin ", i, sep=""))
}
\keyword{models}

