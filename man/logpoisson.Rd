\name{logpoisson}
\alias{logpoisson}
\alias{print.logpoisson}
\alias{summary.logpoisson}
\title{
  Poisson log-linear regression model
}  
\description{
  Fits the poisson log-linear regression model using the maximum-likelihood. The
  log-likelihood is maximized using the Newton-Raphson algorithm (the
  same as Fisher scoring in this case). The
  function returns the inverse of the observed and expected information
  matrix. 
}
\usage{
logpoisson(y, x, offset=0, epsilon=1e-08, maxit=25, trace=FALSE)

\method{print}{logpoisson}(x, \dots)

\method{summary}{logpoisson}(object, \dots)
}
\arguments{
  \item{y}{response vector taking integer values or zero.}
  \item{x}{matrix or data.frame with covarites.

    Intercept is included by default in the model and should not be included
    in \code{x}.
  }
  \item{offset}{possible offset term. It is assumed to be equal to zero
    if not specified.}
  \item{epsilon}{positive convergence tolerance
    \eqn{\varepsilon}{epsilon}. The iterations converge when
    \deqn{\left\vert\frac{\ell_{new} -
	\ell_{old}}{\ell_{new}}\right\vert \leq \varepsilon,}{%
      abs((l[new] - l[old])/l[new]) <= epsilon,}
    where \eqn{\ell}{l} denotes the value of the log-likelihood.
  }
  \item{maxit}{integer giving the maximal number of iterations.}
  \item{trace}{logical indicating if output should be produced for each
    iteration.}
  \item{object}{an object of class "logpoisson".}
  \item{\dots}{other arguments passed to \code{print} or \code{summary}.}  
}
\value{
  An object of class "logpoisson". This has components
  \item{coefficients}{the coefficients of the linear predictor.}
  \item{loglik}{the value of the log-likelihood.}
  \item{score}{the score vector.}
  \item{vcov}{the inverse of the information matrix.}
  \item{linear.predictors}{the values of the linear predictor for each
    observation.}
  \item{fitted.values}{the values of fitted counts for
    each observation.}
  \item{converged}{logical indicating whether the optimization routine converged.}
  \item{iter}{number of iterations performed}
  \item{y}{}
  \item{x}{}  
}
\references{
Agresti, A. (2002). \emph{Categorical Data Analysis. Second
  edition}. Hoboken: John Wiley \& Sons. Section 7.2.
}  
\seealso{
  \code{\link{glm}}.
}  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek@med.kuleuven.be}
}  
\examples{
set.seed(1977)
n <- 100
x1 <- rbinom(n, 1, 0.4)
x2 <- runif(n, 0, 1)
eta <- 5 + 0.1*x1 -0.2*x2
mu <-  exp(eta)
y <- rpois(n, mu)

### Fit the model using poisson
Xmat <- data.frame(x1=x1, x2=x2)
fit <- logpoisson(y=y, x=Xmat)
summary(fit)

### Fit the model using standard glm
fit0 <- glm(y~x1+x2, family=poisson(link="log"))
summary(fit0)
}
\keyword{models}
