\name{glmmAK.files2coda}
\alias{glmmAK.files2coda}
\title{
  Conversion of sampled values into coda mcmc objects
}
\description{
  It takes the values sampled by \code{\link{cumlogitRE}} or
  \code{\link{logpoissonRE}} and stored in *.sim files and converts them
  into coda \code{\link[coda]{mcmc}} objects.
}
\usage{
glmmAK.files2coda(dir, drandom=c("none", "normal", "gspline"),
   quiet=FALSE, skip=0,
   params=list(prob=FALSE, ecount=FALSE, b=FALSE, alloc=FALSE))
}
\arguments{
  \item{dir}{character specifying a directory with sampled values}
  \item{drandom}{string specifying the distribution of random effects
    used in the original \code{\link{cumlogitRE}} or
    \code{\link{logpoissonRE}} call}
  \item{quiet}{logical, passed to \code{\link{scan}} function}
  \item{skip}{number of MCMC iterations that should be skipped at the
    beginning of the chain}
  \item{params}{a list of logical values specifying which of the
    optional sampled parameters should be read}
}
\value{
  A list of coda \code{\link[coda]{mcmc}} objects.
}
\seealso{
  \code{\link[coda]{mcmc}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek@med.kuleuven.be}
}  
%\examples{
%}
\keyword{manip}
