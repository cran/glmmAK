\name{summaryGspline1}
\alias{summaryGspline1}
\title{
  Summary for a univariate G-spline (penalized Gaussian mixture)
}
\description{
  This function is primarily designed to work
  out the MCMC output from functions \code{\link{cumlogitRE}} and
  \code{\link{logpoissonRE}} in which a distribution of the univariate
  random effect was specified as a G-spline. It computes posterior
  pointwise mean and quantiles for a G-spline density based on the MCMC output.
}
\usage{
summaryGspline1(x, mu, sigma,
   standard=TRUE, intcpt, scale,
   probs, values=FALSE,
   dir=getwd(), wfile="logweight.sim", header=TRUE, logw=TRUE,
   skip=0, nwrite)
}
\arguments{
  \item{x}{grid of values at which we want to evaluate the G-spline density}
  \item{mu}{a vector with G-spline knots (means of basis G-splines)}
  \item{sigma}{basis standard deviation(s). If a single number is
    supplied then it is assumed that all basis G-splines have the same
    standard deviation. Alternatively a vector of the same length as
    \code{mu} can be given in which case the basis G-splines do not
    necessarily have the same standard deviations}
  \item{standard}{logical, if \code{TRUE} then the standardized (zero-mean,
    unit-variance) G-spline densities are computed and summarized}
  %\item{nsample}{length of the MCMC sample. It should be equal to the
  %  lengths of the arguments \code{intcpt} and \code{scale} if these are
  %  given.
  %}  
  \item{intcpt}{a vector with sampled intercept values. If not supplied it
    is assumed that all intercepts are equal to zero.    

    It does not have to be supplied if \code{standard=TRUE}.
  }
  \item{scale}{a vector with sampled values of the G-spline scale (in most
    of my papers denoted by \eqn{\tau}{tau}) If not supplied it is
    assumed that all scale values are equal to one.

    It does not have to be supplied if \code{standard=TRUE}.
  }
  \item{probs}{probabilities for which the (pointwise) sample quantiles
    of the G-spline density should be computed.

    If not given only average (and values) of the G-spline are computed}
  \item{values}{if \code{TRUE} also values of the G-spline at each
    (MCMC) iteration are returned.

    If \code{FALSE} only sample mean (and quantiles) of the G-spline are
    returned}
  \item{dir}{character giving the directory where the file with (sampled)
    G-spline (log-)weights is stored}
  \item{wfile}{character giving the name of the file with (sampled)
    G-spline (log-)weights}
  \item{header}{logical indicating whether the file \code{wfile}
    contains a header}
  \item{logw}{logical indicating whether the file \code{wfile} contains
    logarithms of the weights}
  \item{skip}{number of data rows that should be skipped at the beginning of
    the file \code{wfile}}
  %\item{by}{possible additional thinning of the sample stored in
  %  \code{wfile}}
  \item{nwrite}{frequency with which is the user informed about the
    progress of computation (every \code{nwrite}th iteration count of
    iterations change)}    
}
\value{
  A list with the following components (component \samp{values} is
  present only when the argument \code{values} was \code{TRUE}):
  \item{summary}{a \code{data.frame} with the following columns
    (the description below applies to the situation when \code{probs=0.5}):
    \tabular{ll}{
    \samp{x}    \tab a grid of values at which the G-spline density is
      evaluated \cr
    \samp{Mean} \tab pointwise posterior mean of the G-spline density \cr
    \samp{50\%} \tab pointwise posterior 50\% quantile of the G-spline
    density. \cr
       \tab There is one column of this type for each \code{probs} value\cr
    }
  }
  \item{values}{a matrix with one column for each \code{x} value and
    number of rows equal to the length of the MCMC. In each row, there
    is a G-spline density evaluated at one MCMC iteration.}
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
### See ex-Toenail.pdf, ex-Toenail.R
### and ex-Epileptic.pdf, ex-Epileptic.R
### available in the documentation
### to the package
}
\keyword{smooth}

