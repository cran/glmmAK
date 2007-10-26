\name{BPvalue}
\alias{BPvalue}
\title{
  Univariate Bayesian P-value (pseudo-contour probability)
}
\description{
  For a sample (from the posterior distribution) of \eqn{\theta}{theta}
  this function computes
  \deqn{P=2\,\mbox{min}\bigl\{\mbox{P}(\theta<0),\,\mbox{P}(\theta>0)\bigr\},}{%
    P=2*min(P(theta<0, theta>0)),}
  which can be viewed as a counterpart of a classical two-sided
  P-value.

  Note that this is the same as a univariate pseudo-contour probability
  as described in Besag et al. (1995, p. 30) and in Held (2004).
}
\usage{
BPvalue(sample)
}
\arguments{
  \item{sample}{vector, matrix or data frame with sampled values. If it is
    a matrix or data frame then it is assumed that each column
    corresponds to a separate parameter and the P-value is computed
    separately for each column.    
  }
}
\value{
  Vector of computed P-values.
}
\references{
  Besag, J., Green, P., Higdon, D. and Mengersen, K. (1995).
  Bayesian computation and stochastic systems (with Discussion).
  \emph{Statistical Science,} \bold{10}, 3 - 66.

  Held, L. (2004).
  Simultaneous posterior probability statements from Monte Carlo output.
  \emph{Journal of Computational and Graphical Statistics,} \bold{13}, 20 - 35.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{komarek@karlin.mff.cuni.cz}
}  
\examples{
m <- 1000
sample <- rnorm(m, mean=1)
BPvalue(sample)
## compare with
2*pnorm(0, mean=1, lower.tail=TRUE)

sample <- data.frame(x1=rnorm(m), x2=rnorm(m, mean=-1), x3=rnorm(m, mean=2))
BPvalue(sample)
## compare with
2*pnorm(0, mean=0)
2*pnorm(0, mean=-1, lower.tail=FALSE)
2*pnorm(0, mean=2, lower.tail=TRUE)
}
\keyword{htest}
