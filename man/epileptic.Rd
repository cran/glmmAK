\name{epileptic}
\alias{epileptic}
\docType{data}
\title{Epileptic seizures}
\description{
  This data set considers information from a clinical trial of 59 epileptics,
  reported by Thall and Vail (1990). For each patient, the number of epileptic
  seizures was recorded during a baseline period of eight weeks. Patients were
  then randomized to treatment with the anti-epileptic drug progabide, or to
  placebo in addition to standard chemotherapy. The number of seizures was then
  recorded in four consecutive two-weeks intervals.
}
\usage{data(epileptic)}
\format{
  A data frame with 295 observations on the following 5 variables.
  \describe{
    \item{\code{id}}{a unique identifier for the subject in the study.}
    \item{\code{seizure}}{a numeric vector giving the number of epileptic seizures.}
    \item{\code{visit}}{a numeric vector giving the number of the visit, 0=baseline, and
                        1,2,3, and 4 for the four consecutive two-weeks intervals.}
    \item{\code{trt}}{a numeric vector giving the treatment group.}
    \item{\code{age}}{a numeric vector giving the age at the entry.}
  }
}
\source{
  Thall, P. F., and Vail, S. C. (1990)
  Some covariance models for longitudinal count data with
  overdispersion.
  \emph{Biometrics}, \bold{46}, 657--671.
}
\references{
  Breslow, N. E. and Clayton, D. G. (1993).
  Approximate inference in generalized linear mixed models.
  \emph{Journal of the American Statistical Association}, \bold{88},
  9--25.

  Kleinman, K. P. and Ibrahim, J. G. (1998).
  A semi-parametric Bayesian approach to generalized linear mixed
  models.
  \emph{Statistics in Medicine}, \bold{17}, 2579--2596.
}
\seealso{
  \code{\link{epilepticBC}}
}  
\examples{
data(epileptic)
## maybe str(epileptic); plot(epileptic) ...
}
\keyword{datasets}
