\name{cumlogit.help}
\alias{cumlogit.help}
\alias{design.cumlogit}
\alias{fit.cumlogit}
\alias{init00.cumlogit}
\alias{linear.predictors00.cumlogit}
\alias{prob00.cumlogit}
\title{
   Miscellaneous helping functions specific for the cumulative logit model
}
\description{
  These functions are not to be called by ordinary users.

  These are just sub-parts of other functions to make them more
  readable for the programmer.
}
\usage{
design.cumlogit(y, v, x, vb, xb, cluster, intcpt.random=FALSE,
   hierar.center=FALSE, C, predict=FALSE, nobs)

fit.cumlogit(y, v, x, C=1, logit.order=c("decreasing", "increasing"),
   epsilon=1e-08, maxit=25, trace=FALSE)

init00.cumlogit(y, v, x, C, only.intcpt=FALSE,
   logit.order=c("decreasing", "increasing"))

linear.predictors00.cumlogit(coef, v, x, C,
   logit.order=c("decreasing", "increasing"))

prob00.cumlogit(coef, v, x, C,
   logit.order=c("decreasing", "increasing"))
}
\arguments{
  \item{y}{}
  \item{v}{}
  \item{x}{}
  \item{vb}{}
  \item{xb}{}
  \item{cluster}{}
  \item{intcpt.random}{}
  \item{hierar.center}{}
  \item{C}{}
  \item{predict}{}
  \item{nobs}{}
  \item{logit.order}{}
  \item{epsilon}{}
  \item{maxit}{}
  \item{trace}{}
  \item{only.intcpt}{}
  \item{coef}{}
}
\value{
  Some lists.
}
\seealso{
  \code{\link{cumlogit}}, \code{\link{cumlogitRE}}.
}  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{komarek@karlin.mff.cuni.cz}
}
\keyword{internal}
