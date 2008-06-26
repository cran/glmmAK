\name{densplotAK}
\alias{densplotAK}
\title{
  Probability density function estimate from MCMC output
}
\description{
  Displays a plot of the density estimate for each variable in x, calculated by the density function. 

  This is slightly modified version of \code{\link[coda]{densplot}} function
  of a \code{coda} package to conform to my personal preferences.
}
\usage{
densplotAK(x, plot=TRUE, show.obs=FALSE, bwf, bty="n", main ="",
   xlim, ylim, xlab, ylab, \dots) 
}
\arguments{
  \item{x}{an \code{\link[coda]{mcmc}} or \code{\link[coda]{mcmc.list}} object.}
  \item{plot}{if \code{TRUE} this function works more or less in the
    same way as \code{coda} function \code{\link[coda]{densplot}} function.
    If \code{FALSE} this function returns one data frame for each chain
    with computed density which can be used for future plotting.}
  \item{show.obs}{show observations along the x-axis?}
  \item{bwf}{function for calculating the bandwidth. If omitted, the
    bandwidth is calculate by 1.06 times the minimum of the standard
    deviation and the interquartile range divided by 1.34 times the
    sample size to the negative one fifth power.}
  \item{xlim, ylim, xlab, ylab}{further arguments passed to the
    \code{\link{plot.default}} function.}  
  \item{bty, main, \dots}{further arguments passed to the
    \code{\link{plot.default}} function.}
}
\seealso{
  \code{\link[coda]{densplot}}.
}  
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}
\keyword{hplot}
