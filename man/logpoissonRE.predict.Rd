\name{logpoissonRE.predict}
\alias{logpoissonRE.predict}
\title{
  Prediction for Poisson log-linear regression
  with random effects
}
\description{
  This function compute predictive (expected) counts 
  for specified combinations of covariates. It is based on
  the MCMC output obtained using \code{\link{logpoissonRE}}.
}
\usage{
logpoissonRE.predict(nobs, x, xb, offset, cluster,
      intcpt.random=FALSE, hierar.center=FALSE,     
      drandom=c("normal", "gspline"),
      betaF, betaR, varR, is.varR=TRUE,
      prior.gspline,
      probs, values=FALSE,
      dir=getwd(), wfile, indfile, header=TRUE, logw, is.indfile,
      skip=0, nwrite)                             
}
\arguments{
  \item{nobs}{number of covariate combinations for which we want to
    perform a prediction}
  \item{x}{covariate combinations for which we want to perform a
    prediction.

    It should have the same structure as in
    \code{\link{logpoissonRE}} used to obtain the MCMC output}
  \item{xb}{covariate combinations for which we want to perform a
    prediction.

    It should have the same structure as in
    \code{\link{logpoissonRE}} used to obtain the MCMC output}
  \item{offset}{optional offset vector.}
  \item{cluster}{vector defining pertinence of the single observations to
    clusters. It is useful when we want to predict longitudinal profiles.

    See also the same argument in \code{\link{logpoissonRE}}.}
  \item{intcpt.random}{see the same argument in \code{\link{logpoissonRE}}}
  \item{hierar.center}{see the same argument in \code{\link{logpoissonRE}}}
  \item{drandom}{see the same argument in \code{\link{logpoissonRE}}}
%  \item{nsample}{length of the MCMC sample. It should be equal to the
%    lengths or numbers of rows of arguments \code{betaF}, \code{betaR}, \code{varR}}
  \item{betaF}{sampled values of the fixed effects. This should be a
    (sub)sample from the MCMC output stored in the file \file{betaF.sim}}
  \item{betaR}{sampled values of the mean of random effects. This should be a
    (sub)sample from the MCMC output stored in the file \file{betaR.sim}

    It is only needed if \code{hierar.center} is \code{TRUE}.
  }  
  \item{varR}{sampled values of either (co)variance matrices or
    precision (matrices) for random effects if there are any. This
    should be a (sub)sample of either the first or second half of the
    columns stored in the file \file{varR.sim}}
  \item{is.varR}{logical indicating whether \code{varR} gives
    (co)variance (\code{is.varR} \code{TRUE}) or precisions (inverse
    variances) (\code{is.varR} \code{FALSE})}
  \item{prior.gspline}{if \code{drandom} is \bold{gspline} this is a list
    specifying the G-splines. It should have the same structure as the
    same argument in \code{\link{logpoissonRE}} used to obtain the MCMC
    output. However, it is satisfactory if the items
    \bold{K}, \bold{delta} and \bold{sigma} are given.}
  \item{probs}{probabilities for which the (pointwise) sample quantiles
    of the predictive counts should be computed.

    If not given only average (and values) of the predictive
    counts are computed}
  \item{values}{if \code{TRUE} also values of the predictive counts at each
    (MCMC) iteration are returned.

    If \code{FALSE} only sample mean (and quantiles) of the predictive
    probabilities  are returned}  
  \item{dir}{character giving the directory where the file with (sampled)
    G-spline (log-)weights is stored.

    Needed only if \code{drandom} is \bold{gspline}.
  }
  \item{wfile}{character giving the name of the file with (sampled)
    G-spline (log-)weights.

    Needed only if \code{drandom} is \bold{gspline}. In most cases, for
    \bold{univariate} G-spline this argument will be equal to
    \dQuote{\code{logweight.sim}}
    and for \bold{bivariate} G-spline equal to \dQuote{\code{weight.sim}}.
  }
  \item{indfile}{character giving the name of the file where we stored
    indeces of these G-spline components for which the weights are
    stored in the file given by \code{wfile}. The corresponding file
    should have the same structure as \file{knotInd.sim} created by
    \code{\link{logpoissonRE}}.

    Needed only if \code{is.indfile} is \code{TRUE}. In most cases, for
    \bold{univariate} G-spline it does not have to be specified and for
    \bold{bivariate} G-spline it will be equal to \dQuote{\code{knotInd.sim}}.
  }  
  \item{header}{logical indicating whether the files \code{wfile}, \code{indfile}
    contain a header.

    Needed only if \code{drandom} is \bold{gspline}. 
  }
  \item{logw}{logical indicating whether the file \code{wfile} contains
    logarithms of the weights.

    Needed only if \code{drandom} is \bold{gspline}. In most cases, for
    \bold{univariate} G-spline it will be \code{TRUE} and for
    \bold{bivariate} G-spline it will be \code{FALSE}.
  }
  \item{is.indfile}{logical.

    If \code{TRUE} then \code{wfile} contains
    only the non-zero weights and the G-spline is reconstructed using
    \code{indfile}.

    If \code{FALSE} then \code{wfile} must contain on
    each row weights of all components and \code{indfile} is ignored.

    Needed only if \code{drandom} is \bold{gspline} and random effects
    are bivariate.
  } 
  \item{skip}{number of data rows that should be skipped at the beginning of
    the files \code{wfile}, \code{indfile}.
  }
%  \item{by}{possible additional thinning of the sample stored in
%    \code{wfile}.
%  }
  \item{nwrite}{frequency with which is the user informed about the
    progress of computation (every \code{nwrite}th iteration count of
    iterations change)}    
}  %%% end of arguments
\value{
  A list with the following components (description below applies for
  the case with \code{prob=0.5})
  \item{Mean}{a matrix with \eqn{1} column giving in each row
    posterior predictive mean of the count \eqn{\mbox{E}(Y)}{E(Y)}
    for a given covariate combination.
  }
  \item{50\%}{a matrix with \eqn{1} column giving in each row
    posterior predictive quantile (here 50\% quantile) of the count
    for a given covariate combination.

    There is one component of this type in the resulting \code{list}
    for each value of \code{probs}.
  }
  \item{values}{a matrix with \eqn{n} columns, where
    \eqn{n} denotes the number of covariate combinations for which we
    perform the prediction, and number of rows equal to the length of
    the MCMC. Each column gives sampled counts a given covariate combination.

    It is returned only if \code{values} is \code{TRUE}.
  }  
}
\references{
  Komárek, A. and Lesaffre, E. (2008).
  Generalized linear mixed model with a penalized Gaussian mixture as a
  random-effects distribution.
  \emph{Computational Statistics and Data Analysis}, \bold{52}, 3441--3458.
}
\seealso{
  \code{\link{logpoissonRE}}, \code{\link{logpoisson}}, \code{\link{glm}}.
}  
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\keyword{models}


