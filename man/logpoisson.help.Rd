\name{logpoisson.help}
\alias{logpoisson.help}
\alias{design.logpoisson}
\alias{fit.logpoisson}
\title{
   Miscellaneous helping functions specific for the cumulative logit model
}
\description{
  These functions are not to be called by ordinary users.

  These are just sub-parts of other functions to make them more
  readable for the programmer.
}
\usage{
design.logpoisson(y, x, xb, offset, cluster, intcpt.random=FALSE,
   hierar.center=FALSE, predict=FALSE, nobs)

fit.logpoisson(y, x, offset, epsilon=1e-08, maxit=25, trace=FALSE)
}
\arguments{
  \item{y}{}
  \item{x}{}
  \item{xb}{} 
  \item{offset}{}
  \item{cluster}{}
  \item{intcpt.random}{}
  \item{hierar.center}{}
  \item{predict}{}
  \item{nobs}{}
  \item{epsilon}{}
  \item{maxit}{}
  \item{trace}{}    
}
\seealso{
  \code{\link{logpoisson}}, \code{\link{logpoissonRE}}.
}  
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}
\keyword{internal}
