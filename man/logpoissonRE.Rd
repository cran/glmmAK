\name{logpoissonRE}
\alias{logpoissonRE}
\title{
  Poisson log-linear regression
  with random effects
}
\description{
  This function implements MCMC sampling for the Poisson log-linear model.
  Details are given in \eqn{\mbox{Kom\'arek}}{Kom&#225rek} and
  Lesaffre (2007). On as many places as possible, the same notation as
  in this paper is used also in this manual page.

  In general, the following log-linear model for response \eqn{Y} is assumed:
  \deqn{\log\bigl\{\mbox{E}(Y)\bigr\} = \eta,}{\log(Y) = eta,}
  where the form of the linear predictor
  \eqn{\eta}{eta} depends on whether a
  hierarchical centering is used or not. In the following,
  \eqn{\beta}{beta} denotes fixed effects and \eqn{b}{b} random effects.

  \bold{No hierarchical centering (DEFAULT)}\cr
  The linear predictor has the following form
  \deqn{\eta = \beta'(x',\,x_b') + b'x_b,}{%
    eta= beta'(x', x(b)') + b'x(b),}
  where \eqn{b} is a vector of random effects with zero location. 
  
  \bold{Hierarchical centering}\cr
  The linear predictor has the following form
  \deqn{\eta = \beta'x + b'x_b,}{%
    eta= beta'x + b'x(b),}
  \eqn{b} is a vector of random effects with location \eqn{\alpha}{alpha}.

  For description of the rest of the model, see \code{\link{cumlogitRE}}.      
}
\usage{
logpoissonRE(y, x, xb, offset=0, cluster,                       
      intcpt.random=FALSE,
      hierar.center=FALSE,                       
      drandom=c("normal", "gspline"),
      prior.fixed,
      prior.random,
      prior.gspline,
      init.fixed,
      init.random,
      init.gspline,                 
      nsimul = list(niter=10, nthin=1, nburn=0, nwrite=10),
      store = list(ecount=FALSE, b=FALSE, alloc=FALSE, acoef=FALSE),
      dir=getwd(),
      precision=8)
}
\arguments{
  \item{y}{response vector taking integer values or zero.
  }
  \item{x}{vector, matrix or data.frame with covariates for \bold{fixed}
    effects.
  }
  \item{xb}{vector, matrix or data.frame with covariates for \bold{random}
    effects.

    If you want to include \bold{random intercept}, do it by setting the argument
    \code{intcpt.random} to \code{TRUE}. The intercept column should not
    be included in \code{xb}.    
  }
  \item{offset}{optional vector of the offset term. }
  \item{cluster}{see \code{\link{cumlogitRE}}.}
  \item{intcpt.random}{see \code{\link{cumlogitRE}}.}
  \item{hierar.center}{see \code{\link{cumlogitRE}}.}
  \item{drandom}{see \code{\link{cumlogitRE}}.}
  \item{prior.fixed}{see \code{\link{cumlogitRE}}.}
  \item{prior.random}{see \code{\link{cumlogitRE}}.}
  \item{prior.gspline}{see \code{\link{cumlogitRE}}.}  
  \item{init.fixed}{see \code{\link{cumlogitRE}}.}
  \item{init.random}{see \code{\link{cumlogitRE}}.}
  \item{init.gspline}{see \code{\link{cumlogitRE}}.}
  \item{nsimul}{see \code{\link{cumlogitRE}}.}
  \item{store}{list indicating which chains (out of these not stored by default) should be compulsory
    stored. The list has the logical components with the following
    names.
    \describe{
      \item{ecount}{if \code{TRUE} values of individual predictive
	(expected) counts are stored.
      }
      \item{b}{if \code{TRUE} values of cluster specific random effects
	are stored.
      }
      \item{alloc}{if \code{TRUE} values of allocation indicators are stored.
      }
      \item{acoef}{if \code{TRUE} and distribution of random effects is
	given as a \bold{bivariate} G-spline values of log-G-spline
	weights (\eqn{a}{a} coefficients) are stored for all components.
      }
    }      
  }
  \item{dir}{see \code{\link{cumlogitRE}}.}
  \item{precision}{see \code{\link{cumlogitRE}}.}
}
\value{
  See \code{\link{cumlogitRE}.}
}
\section{Files created}{
  \describe{
    \item{iteration.sim}{see \code{\link{cumlogitRE}}.}
    \item{betaF.sim}{
      sampled values of the fixed effects \eqn{\beta}{beta}.

      \bold{Note} that in models with \bold{G-spline} distributed random
      effects which are not hierarchically centered, the average effect
      of the covariates involved in the random effects (needed for
      inference) is obtained as a
      sum of the corresponding \eqn{\beta}{beta} coefficient and a
      scaled mean of the G-spline. \eqn{\beta}{beta} coefficients
      adjusted in this way are stored in the file \file{betaRadj.sim}
      (see below).
    }
    \item{betaR.sim}{sampled values of the location parameters \eqn{\alpha}{alpha}
      of the random effects when the \bold{hierarchical centering} was
      used.

      \bold{Note} that in models with \bold{G-spline} distributed random
      effects which are hierarchically centered, the average effect
      of the covariates involved in the random effects (needed for
      inference) is obtained as a
      sum of the corresponding \eqn{\alpha}{alpha} coefficient and a mean
      of the G-spline.  \eqn{\alpha}{alpha} coefficients adjusted in
      this way are stored in the file \file{betaRadj.sim} (see below).      
    }
    \item{varR.sim}{see \code{\link{cumlogitRE}}.}
    \item{loglik.sim}{see \code{\link{cumlogitRE}}.}
    \item{expectcount.sim}{
      sampled values of predictive (expected) counts for each observations.
      
      Created only if \code{store$ecount} is \code{TRUE}.
    }  
    \item{b.sim}{see \code{\link{cumlogitRE}}.}
  }
}  
\section{Files created for models with G-spline distributed
  random effects}{
  See \code{\link{cumlogitRE}}.
}
\references{
  Agresti, A. (2002).
  \emph{Categorical Data Analysis. Second edition}.
  Hoboken: John Wiley \& Sons.

  Gelfand, A. E., Sahu, S. K., and Carlin, B. P. (1995).
  Efficient parametrisations for normal linear mixed models.
  \emph{Biometrika}, \bold{82}, 479--488.

  Gilks, W. R. and Wild, P. (1992).
  Adaptive rejection sampling for Gibbs sampling.
  \emph{Applied Statistics,} \bold{41}, 337--348.

  Neal, R. M. (2003).
  Slice sampling (with Discussion).
  \emph{The Annals of Statistics,} \bold{31}, 705--767.
  
  \eqn{\mbox{Kom\'arek, A.}}{Kom&#225rek, A.} and Lesaffre, E. (2007).
  Generalized linear mixed model with a penalized Gaussian mixture as a
  random-effects distribution.
  \emph{SOME JOURNAL}, \bold{999}, 000--000.
  
  Molenberghs, G. and Verbeke, G. (2005).
  \emph{Models for Discrete Longitudinal Data}.
  New York: Springer Science+Business Media.
}  
\seealso{
  \code{\link{logpoisson}}, \code{\link{cumlogitRE}}, \code{\link{poisson}}, \code{\link{glm}}.
}  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek@med.kuleuven.be}
}  
\examples{
### See ex-Epileptic.pdf and ex-Epileptic.R
### available (soon) in the documentation
### to the package
}
\keyword{models}
