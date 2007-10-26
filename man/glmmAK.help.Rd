\name{glmmAK.help}
\alias{glmmAK.help}
\alias{checknsimul.glmmAK}
\alias{prior.fixed.glmmAK}
\alias{init.fixed.glmmAK}
\alias{prior.random.glmmAK}
\alias{init.random.glmmAK}
\alias{prior.gspline.glmmAK}
\alias{init.gspline.glmmAK}
\alias{headers.glmmAK}
\title{
   Miscellaneous helping functions.
}
\description{
  These functions are not to be called by ordinary users.

  These are just sub-parts of other functions to make them more
  readable for the programmer.
}
\usage{
checknsimul.glmmAK(nsimul)

prior.fixed.glmmAK(prior.fixed, nFixed)

init.fixed.glmmAK(init.fixed, init.fixedML, nFixed)

prior.random.glmmAK(prior.random, nRandom,
   drandom=c("normal", "gspline"), hierar.center)

init.random.glmmAK(init.random, init.MrandomML, init.DrandomML,
   nRandom, drandom=c("normal", "gspline"), N, hierar.center)

prior.gspline.glmmAK(prior.gspline, nRandom,
   drandom=c("normal", "gspline"), simplified=FALSE)

init.gspline.glmmAK(init.gspline, prior.gspline, init.random,
   nRandom, drandom=c("normal", "gspline"))

headers.glmmAK(dir, store, DES, prob.init,
   init.fixedML, init.MrandomML, hierar.center,
   drandom=c("normal", "gspline"), model=c("cumlogit", "logpoisson"),
   prior.gspline)
}
\arguments{
  \item{nsimul}{a list}
  \item{prior.fixed}{}
  \item{init.fixed}{}
  \item{init.fixedML}{}  
  \item{nFixed}{}
  \item{prior.random}{}
  \item{init.random}{}
  \item{init.randomML}{}    
  \item{init.DrandomML}{}
  \item{init.MrandomML}{}    
  \item{nRandom}{}
  \item{drandom}{}
  \item{model}{}
  \item{hierar.center}{}
  \item{N}{}
  \item{prior.gspline}{}  
  \item{init.gspline}{}
  \item{simplified}{}
  \item{dir}{}
  \item{store}{}
  \item{DES}{}  
}
\value{
  Some lists.
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{komarek@karlin.mff.cuni.cz}
}
\keyword{internal}

