\name{epilepticBC}
\alias{epilepticBC}
\docType{data}
\title{Epileptic seizures, Breslow and Clayton transformed}
\description{
  This data set considers information from a clinical trial of 59 epileptics,
  reported by Thall and Vail (1990). For each patient, the number of epileptic
  seizures was recorded during a baseline period of eight weeks. Patients were
  then randomized to treatment with the anti-epileptic drug progabide, or to
  placebo in addition to standard chemotherapy. The number of seizures was then
  recorded in four consecutive two-weeks intervals.

  This is a transformed version of the original data
  (\code{\link{epileptic}}) which allows directly to fit models
  described in Breslow and Clayton (1993) or in Kleinman and Ibrahim
  (1998).

  There are 4 rows in the dataset for each of 59 patients.
}
\usage{data(epilepticBC)}
\format{
  A data frame with 236 observations on the following 10 variables.
  \describe{
    \item{\code{id}}{a unique identifier for the subject in the study.}
    \item{\code{visit}}{a numeric vector giving the number of the visit,
                        1,2,3, and 4 for the four two-weeks intervals
			following the baseline measurement.}
    \item{\code{seizure0}}{a numeric vector giving the number of
      epileptic seizures in the eight-week period prior the randomization.}
    \item{\code{age}}{a numeric vector giving the age (in years) at the entry.}    
    \item{\code{Seizure}}{a numeric vector giving the number of
      epileptic seizures in a given two-week interval.}
    \item{\code{Base}}{a numeric vector giving the transformed
      \code{seizure0}. \code{Base}=log(\code{seizure0}/4).}
    \item{\code{Trt}}{a numeric vector giving the treatment group.}    
    \item{\code{Base.Trt}}{a numeric vector giving the interaction
      covariate \code{Base}*\code{Trt}.}
    \item{\code{Age}}{a numeric vector giving the transformed
      age. \code{Age}=log(\code{age}).}
    \item{\code{Visit}}{a numeric vector giving the centered
      \code{visit}. \code{Visit}=(2*\code{visit}-5)/10.}    
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
  \code{\link{epileptic}}
}  
\examples{
data(epilepticBC)
## maybe str(epilepticBC); plot(epilepticBC) ...
}
\keyword{datasets}
