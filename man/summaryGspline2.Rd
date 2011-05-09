\name{summaryGspline2}
\alias{summaryGspline2}
\title{
  Summary for a bivariate G-spline (penalized Gaussian mixture)
}
\description{
  This function is primarily designed to work
  out the MCMC output from functions \code{\link{cumlogitRE}} and
  \code{\link{logpoissonRE}} in which a distribution of the bivariate
  random effect was specified as a G-spline. It computes posterior
  pointwise mean and quantiles for a G-spline density based on the MCMC output.
  
  Besides the summary for the joint bivariate G-spline density it also
  directly computes summaries for both marginal G-splines.
}
\usage{
summaryGspline2(x1, x2, mu1, mu2, sigma1, sigma2,
   standard=TRUE, intcpt, scale,
   probs, values=FALSE,
   dir=getwd(), wfile="weight.sim", indfile="knotInd.sim",
   header=TRUE, logw=FALSE, is.indfile=TRUE, 
   skip=0, nwrite)
}
\arguments{
  \item{x1}{grid of values for the first margin at which we evaluate the
    G-spline density}
  \item{x2}{grid of values for the second margin at which we evaluate the G-spline density}  
  \item{mu1}{a vector with G-spline knots (means of basis G-splines) for
    the first margin}
  \item{mu2}{a vector with G-spline knots (means of basis G-splines) for
    the second margin}
  \item{sigma1}{basis standard deviation(s) for the first margin. If a single number is
    supplied then it is assumed that all basis G-splines in the first
    margin have the same standard deviation. Alternatively a vector of the same length as
    \code{mu1} can be given in which case the basis G-splines in the
    first margin do not necessarily have the same standard deviations}
  \item{sigma2}{basis standard deviation(s) for the second margin. If a single number is
    supplied then it is assumed that all basis G-splines in the second
    margin have the same standard deviation. Alternatively a vector of the same length as
    \code{mu2} can be given in which case the basis G-splines in the
    second margin do not necessarily have the same standard deviations}
  \item{standard}{if \code{TRUE} then the standardized (zero-mean,
    unit-variance) G-splines are computed and summarized}
  %\item{nsample}{length of the MCMC sample. It should be equal to the
  %  lengths of the arguments \code{intcpt} and \code{scale} if these are
  %  given.
  %}  
  \item{intcpt}{a two-column matrix with sampled intercept values. If not supplied it
    is assumed that all intercepts are equal to zero.

    It does not have to be supplied if \code{standard=TRUE}.
  }
  \item{scale}{a two-column matrix with sampled intercept values of the G-spline scale (in most
    of my papers denoted by \eqn{\tau}{tau}) If not supplied it is
    assumed that all scale values are equal to one.

    It does not have to be supplied if \code{standard=TRUE}.
  }
  \item{probs}{probabilities for which the (pointwise) sample quantiles
    of the G-spline should be computed.

    If not given only average (and values) of the G-spline are computed}
  \item{values}{if \code{TRUE} also values of the G-spline at each
    (MCMC) iteration are returned.

    If \code{FALSE} only sample mean (and quantiles) of the G-spline are
    returned}
  \item{dir}{character giving the directory where the file with (sampled)
    G-spline (log-)weights is stored}
  \item{wfile}{character giving the name of the file with (sampled)
    G-spline (log-)weights}
  \item{indfile}{character giving the name of the file where it is
    indicated which G-spline components correspond to non-zero weights.

    It does not have to be supplied if \code{is.indfile=FALSE}.
  }  
  \item{header}{logical indicating whether the files \code{wfile}, \code{indfile}
    contain a header}
  \item{logw}{logical indicating whether the file \code{wfile} contains
    logarithms of the weights}
  \item{is.indfile}{logical.

    If \code{TRUE} then \code{wfile} contains
    only the non-zero weights and the G-spline is reconstructed using
    \code{indfile}.

    If \code{FALSE} then \code{wfile} must contain on
    each row weights of all components and \code{indfile} is ignored.}
  \item{skip}{number of data rows that should be skipped at the beginning of
    the files \code{wfile}, \code{indfile}}
  %\item{by}{possible additional thinning of the sample stored in
  %  \code{wfile}}
  \item{nwrite}{frequency with which is the user informed about the
    progress of computation (every \code{nwrite}th iteration count of
    iterations change)}    
}
\value{
  A list with the following components
  (component \samp{values}, \samp{values1}, \samp{values2} are
  present only when the argument \code{values} was \code{TRUE}).
  The description below applies to the situation when \code{probs=0.5}.
  \item{summary}{a list with the components
    \tabular{ll}{
      \samp{x1} \tab grid of values for the first margin at which we evaluate the G-spline
      density \cr
      \samp{x2} \tab grid of values for the second margin at which we evaluate the G-spline
      density \cr
      \samp{Mean} \tab a matrix with \code{length(x1)} rows and
      \code{length(x2)} columns giving \cr
      \tab the pointwise posterior mean of the joint G-spline density \cr
      \samp{50\%} \tab a matrix with \code{length(x1)} rows and
      \code{length(x2)} columns \cr
      \tab giving the pointwise posterior 50\% quantile of
      the joint G-spline density. \cr
      \tab There is a matrix of this type for each \code{probs} value. \cr      
    }  
  }
  \item{summary1}{a \code{data.frame} with the following columns
    \tabular{ll}{
    \samp{x}    \tab a grid of values at which the first marginal G-spline density is
      evaluated \cr
    \samp{Mean} \tab pointwise posterior mean of the first marginal G-spline density \cr
    \samp{50\%} \tab pointwise posterior 50\% quantile of the first
    marginal G-spline density. \cr
       \tab There is one column of this type for each \code{probs} value\cr
    }
  }
  \item{summary2}{a \code{data.frame} for the second marginal G-spline
    having the same structure as \samp{summary1}.}    
  \item{values}{a matrix with one column for each (\code{x1}, \code{x2}) value and
    number of rows equal to the length of the MCMC. In each row, there
    is a joint G-spline density evaluated at one MCMC iteration.}
  \item{values1}{a matrix with one column for each \code{x1} value and
    number of rows equal to the length of the MCMC. In each row, there
    is the first marginal G-spline density evaluated at one MCMC
    iteration.}
  \item{values1}{a matrix with one column for each \code{x1} value and
    number of rows equal to the length of the MCMC. In each row, there
    is the second marginal G-spline density evaluated at one MCMC iteration.}
}
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
### See ex-Toenail.pdf, ex-Toenail.R
### and ex-Epileptic.pdf, ex-Epileptic.R
### available in the documentation
### to the package
}
\keyword{smooth}

