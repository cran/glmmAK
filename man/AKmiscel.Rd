\name{AKmiscel}
\alias{AKmiscel}
\alias{mfLand}
\alias{mfPort}
\title{
  Miscellaneous smaller functions
}
\description{
  Functions \code{mfLand} and \code{mfPort} try to split a graphical
  device in a ``nice'' way to produce several plots on 1 page in
  landscape or portrait format.
}
\usage{
mfLand(np)

mfPort(np)
}
\arguments{
\item{np}{number of plots that are to be produced on 1 page}
}
\value{
A 2-component vector giving the number of rows and columns into which
the graphical device could be splitted.
}
\seealso{
  \code{\link{par}}.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{komarek@karlin.mff.cuni.cz}
}  
\examples{
mfLand(6)
mfPort(6)
}
\keyword{dplot}
