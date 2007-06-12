\name{cumlogit}
\alias{cumlogit}
\alias{print.cumlogit}
\alias{summary.cumlogit}
\title{
  Cumulative logit model for ordinal responses
}  
\description{
  Fits the cumulative logit model using the maximum-likelihood. The
  log-likelihood is maximized using the Newton-Raphson algorithm. The
  function returns the inverse of both observed and expected information
  matrix. Summary of the model produced by the \code{summary} function
  uses by default the inverse of the observed information matrix for the
  inference. This can be changed by the user such that the expected
  information is used instead.
}
\usage{
cumlogit(y, v, x, C=1, logit.order=c("decreasing", "increasing"),
   epsilon=1e-08, maxit=25, trace=FALSE)

\method{print}{cumlogit}(x, vcov=c("observed", "expected"), \dots)

\method{summary}{cumlogit}(object, vcov=c("observed", "expected"), \dots)
}
\arguments{
  \item{y}{response vector taking values \eqn{0, 1, \dots, C}{0,1,...,C}.}
  \item{v}{matrix or data.frame with covarites whose effect does not
    necessarily satisfy proportional odds assumption.

    Intercept is included by default in the model and should be included neither in
    \code{x}, nor in \code{v}.
  }
  \item{x}{vector, matrix or data.frame with covarites whose effect is
    assumed to satisfy proportional odds assumption.
  }  
  \item{C}{number of response categories minus 1.}
  \item{logit.order}{either "decreasing" or "increasing" indicating in
    which direction the logits are formed.

    For \code{logit.order="decreasing"}:   
        \tabular{rcl}{
          \eqn{\log\Bigl\{\frac{P(Y\geq 1)}{P(Y=0)}\Bigr\}}{log[P(Y>=1)/P(Y=0)]} \tab
	  \eqn{=}{=} \tab
	  \eqn{\beta'x + \gamma_1'v}{beta'x + gamma[1]'v} \cr
          \eqn{\log\Bigl\{\frac{P(Y\geq 2)}{P(Y\leq 1)}\Bigr\}}{log[P(Y>=2)/P(Y<=1)]} \tab
	  \eqn{=}{=} \tab
	  \eqn{\beta'x + \gamma_2'v}{beta'x + gamma[2]'v} \cr
	  \tab \eqn{\vdots}{...} \tab \cr
          \eqn{\log\Bigl\{\frac{P(Y=C)}{P(Y\leq C-1)}\Bigr\}}{log[P(Y=C)/P(Y<=C-1)]} \tab
	  \eqn{=}{=} \tab
	  \eqn{\beta'x + \gamma_C'v}{beta'x + gamma[C]'v}	  
	}

    For \code{logit.order="increasing"}:   	
        \tabular{rcl}{
          \eqn{\log\Bigl\{\frac{P(Y=0)}{P(Y\geq 1)}\Bigr\}}{log[P(Y=0)/P(Y>=1)]} \tab
	  \eqn{=}{=} \tab
	  \eqn{\beta'x + \gamma_1'v}{beta'x + gamma[1]'v} \cr
          \eqn{\log\Bigl\{\frac{P(Y\leq1)}{P(Y\geq 2)}\Bigr\}}{log[P(Y<=1)/P(Y>=2)]} \tab
	  \eqn{=}{=} \tab
	  \eqn{\beta'x + \gamma_2'v}{beta'x + gamma[2]'v} \cr
	  \tab \eqn{\vdots}{...} \tab \cr
          \eqn{\log\Bigl\{\frac{P(Y<=C-1)}{P(Y=C)}\Bigr\}}{log[P(Y<=C-1)/P(Y=C)]} \tab
	  \eqn{=}{=} \tab
	  \eqn{\beta'x + \gamma_C'v}{beta'x + gamma[C]'v}	  
	}      	
  }
  \item{vcov}{character indicating which type of the information matrix
    should be used for the inference}    
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
  \item{object}{an object of class "cumlogit".}
  \item{\dots}{other arguments passed to \code{print} or \code{summary}.}
}
\value{
  An object of class "cumlogit". This has components
  \item{coefficients}{the coefficients of the linear predictor.}
  \item{loglik}{the value of the log-likelihood.}
  \item{score}{the score vector.}
  \item{vcov}{the inverse of the observed information matrix.}
  \item{expect.vcov}{the inverse of the expected information matrix.}
  \item{logit.order}{character indicating the way in which the logits
    are formed.}
  \item{linear.predictors}{the values of the linear predictor for each
    observation and each logit.}
  \item{fitted.values}{the values of fitted category probabilities for
    each observation.}
  \item{converged}{logical indicating whether the optimization routine converged.}
  \item{iter}{number of iterations performed}
  \item{C}{see the function argument}
  \item{y}{see the function argument}
  \item{v}{see the function argument}
  \item{x}{see the function argument}  
}
\references{
Agresti, A. (2002). \emph{Categorical Data Analysis. Second
  edition}. Hoboken: John Wiley \& Sons. Section 7.2.
}  
\seealso{
  \code{\link{glm}}, \code{\link[MASS]{polr}}.
}  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek@med.kuleuven.be}
}  
\examples{
## Simulate the data for C = 3
## ============================
set.seed(775988621)
N3 <- 1000

  ## covariates:
x3 <- data.frame(x1=rnorm(N3), x2=runif(N3, 0, 1))
v3 <- data.frame(v1=rnorm(N3), v2=runif(N3, 0, 1))

  ## regression coefficients
  ## (theta = c(0.1, -0.2,   1, 0.1, -0.2,   0, 0.05, -0.25,   -0.5, 0.05, -0.25)):
alpha3 <- c(1, 0, -0.5)     
beta3 <- c(0.1, -0.2)
gamma3 <- list(c(0.1, -0.2), c(0.05, -0.25), c(0.05, -0.25))

 ## linear predictors and inverse logits:
eta3 <- data.frame(eta1=alpha3[1] + as.matrix(x3) \%*\% beta3 +
                        as.matrix(v3) \%*\% gamma3[[1]],
                   eta2=alpha3[2] + as.matrix(x3) \%*\% beta3 +
                        as.matrix(v3) \%*\% gamma3[[2]],
                   eta3=alpha3[3] + as.matrix(x3) \%*\% beta3 +
                        as.matrix(v3) \%*\% gamma3[[3]])                    
ilogit3 <- data.frame(ilogit1 = exp(eta3[,1])/(1 + exp(eta3[,1])),
                      ilogit2=exp(eta3[,2])/(1 + exp(eta3[,2])),
                      ilogit3=exp(eta3[,3])/(1 + exp(eta3[,3])))                         

 ## category probabilities:
pis3 <- data.frame(pi0=1-ilogit3[,1],
                   pi1=ilogit3[,1]-ilogit3[,2],
                   pi2=ilogit3[,2]-ilogit3[,3],
                   pi3=ilogit3[,3])

 ## response:
rtemp <- function(prob, C){
  return(sample(0:C, size=1, prob=prob))
}
y3 <- apply(pis3, 1, rtemp, C=3)

## Fit the model
## ===============
fit <- cumlogit(y=y3, x=x3, v=v3, C=3)
print(fit)

fit2 <- cumlogit(y=y3, x=x3, v=v3, C=3, logit.order="increasing")
print(fit2)
}
\keyword{models}
