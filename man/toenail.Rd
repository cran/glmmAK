\name{toenail}
\alias{toenail}
\docType{data}
\title{Toenail infection}
\description{
  This data set considers information from
  a longitudinal clinical trial in dermatology which was set up to
  compare the efficacy of two oral treatments for toenail infection (De Backer at al., 1998). 
  One of the end points of the study was the degree of onycholysis which expresses the degree of separation of
  the nail plate from the nail-bed (0, absent; 1, mild; 2, moderate; 3,
  severe) and was evaluated at seven visits
  (approximately on weeks 0, 4, 8, 12, 24, 36 and 48). In total, 1 908 measurements on 294 patients are available.
  In this dataset, only a dichotomized onycholysis (0, absent or mild;
  1, moderate or severe) is given (variable \code{infect}).

  \bold{IMPORTANT NOTICE:} The data have kindly been made available by
  Novartis, Belgium. The source of the data must be acknowledged in any
  publication which uses them
  (see Lesaffre and Spiessens, 2001 for more details).
}
\usage{data(toenail)}
\format{
  A data frame with 1908 observations on the following 5 variables.
  \describe{
    \item{\code{idnr}}{identification number of the patient}
    \item{\code{infect}}{binary response}
    \item{\code{trt}}{treatment group}
    \item{\code{time}}{time of measurement (in months)}
    \item{\code{visit}}{visit number}    
  }  
}
\source{
  Lesaffre, E. and Spiessens, B. (2001).
  On the effect of the number of quadrature points in a logistic
  random-effects model: An example.
  \emph{Applied Statistics}, \bold{50}, 325--335.
}
\references{
  De Backer, M., De Vroey, C., Lesaffre, E., Scheys, I., and De Keyser,
  P. (1998).
  Twelve weeks of continuous onychomycosis caused by dermatophytes: A double blind comparative trial of terbafine
  250 mg/day versus itraconazole 200 mg/day.
  \emph{Journal of the American Academy of Dermatology}, \bold{38}, S57--S63.  
}
\examples{
data(toenail)
## maybe str(toenail); plot(toenail) ...
}
\keyword{datasets}
