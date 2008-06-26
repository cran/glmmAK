\name{copula}
\alias{copula}
\alias{Cplackett}
\alias{Cgauss}
\alias{Cclayton}
\alias{cplackett}
\alias{cgauss}
\alias{cclayton}
\title{
  Copulas
}
\description{
  Functions to compute the cumulative distribution functions and
  densities for several bivariate copulas.

  These functions do not have anything to do with the GLMM's in this
  package. They are here simply because of an interest of the author to
  play with copulas a little bit.
}
\usage{
Cplackett(u, v, theta=1)
Cgauss(u, v, theta=0)
Cclayton(u, v, theta=0)
cplackett(u, v, theta=1)
cgauss(u, v, theta=0)
cclayton(u, v, theta=0)
}
\arguments{
  \item{u}{Unif(0, 1) quantiles for the first margin. It can be a vector
    or a matrix}
  \item{v}{Unif(0, 1) quantiles for the second margin. It can be a
    vector or a matrix}  
  \item{theta}{value of the association parameter}
}
\value{
  A vector or a matrix with the values of the corresponding cdf or density.
}
\references{
  Nelsen, R. B. (2006).
  \emph{An Introduction to Copulas,} Second Edition
  New York: Springer.  
}
\author{
  \eqn{\mbox{Arno\v{s}t Kom\'arek}}{Arno&#353t Kom&#225rek} \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
### Margin1 = normal mixture
### Margin2 = normal
intcpt <- c(1, 5)
sds <- c(1, 3)
x <- seq(intcpt[1]-2.5*sds[1], intcpt[1]+1+2.5*sds[1], length=40)
y <- seq(intcpt[2]-2.5*sds[2], intcpt[2]+2.5*sds[2], length=41)
Fx <- 0.6*pnorm(x, mean=intcpt[1], sd=sds[1]) +
   0.4*pnorm(x, mean=intcpt[1]+2, sd=0.5*sds[1])
Fy <- pnorm(y, mean=intcpt[2], sd=sds[2])
fx <- 0.6*dnorm(x, mean=intcpt[1], sd=sds[1]) +
   0.4*dnorm(x, mean=intcpt[1]+2, sd=0.5*sds[1])
fy <- dnorm(y, mean=intcpt[2], sd=sds[2])
u <- matrix(rep(Fx, length(y)), ncol=length(y))
v <- matrix(rep(Fy, length(x)), nrow=length(x), byrow=TRUE)
du <- matrix(rep(fx, length(y)), ncol=length(y))
dv <- matrix(rep(fy, length(x)), nrow=length(x), byrow=TRUE)


### Copula distribution functions
theta <- c(3, 0.3, 1)
CC <- list()
CC$plackett <- Cplackett(u, v, theta=theta[1])
CC$gauss <- Cgauss(u, v, theta=theta[2])
CC$clayton <- Cclayton(u, v, theta=theta[3])

### Copula densities
cc <- list()
cc$plackett <- cplackett(u, v, theta=theta[1]) * du * dv
cc$gauss <- cgauss(u, v, theta=theta[2]) * du * dv
cc$clayton <- cclayton(u, v, theta=theta[3]) * du * dv

### Figures
lcol <- "red"
pcol <- "seagreen2"
mains <- paste(c("Plackett ", "Gauss ", "Clayton "), "copula, theta=", theta, sep="")
zlab <- "F(x,y)"
zlab2 <- "f(x,y)"
tangle <- -25
ptangle <- 40

oldpar <- par(bty="n", mfcol=c(2, 3))
contour(x, y, cc$plackett, main=mains[1], col=lcol)
persp(x, y, cc$plackett, zlab=zlab2, main=mains[1], col=pcol, theta=tangle, phi=ptangle)
contour(x, y, cc$gauss, main=mains[2], col=lcol)
persp(x, y, cc$gauss, zlab=zlab2, main=mains[2], col=pcol, theta=tangle, phi=ptangle)
contour(x, y, cc$clayton, main=mains[3], col=lcol)
persp(x, y, cc$clayton, zlab=zlab2, main=mains[3], col=pcol, theta=tangle, phi=ptangle)

par(bty="n", mfcol=c(2, 3))
contour(x, y, CC$plackett, main=mains[1], col=lcol)
persp(x, y, CC$plackett, zlab=zlab, main=mains[1], col=pcol, theta=tangle, phi=ptangle)
contour(x, y, CC$gauss, main=mains[2], col=lcol)
persp(x, y, CC$gauss, zlab=zlab, main=mains[2], col=pcol, theta=tangle, phi=ptangle)
contour(x, y, CC$clayton, main=mains[3], col=lcol)
persp(x, y, CC$clayton, zlab=zlab, main=mains[3], col=pcol, theta=tangle, phi=ptangle)

par(bty="n")
layout(matrix(c(0,1,1,0,  2,2,3,3), nrow=2, byrow=TRUE))
contour(x, y, cc$plackett, main=mains[1], col=lcol)
contour(x, y, cc$gauss, main=mains[2], col=lcol)
contour(x, y, cc$clayton, main=mains[3], col=lcol)

par(bty="n")
layout(matrix(c(0,1,1,0,  2,2,3,3), nrow=2, byrow=TRUE))
persp(x, y, cc$plackett, zlab=zlab2, main=mains[1], col=pcol, theta=tangle, phi=ptangle)
persp(x, y, cc$gauss, zlab=zlab2, main=mains[2], col=pcol, theta=tangle, phi=ptangle)
persp(x, y, cc$clayton, zlab=zlab2, main=mains[3], col=pcol, theta=tangle, phi=ptangle)

par(bty="n")
layout(matrix(c(0,1,1,0,  2,2,3,3), nrow=2, byrow=TRUE))
contour(x, y, CC$plackett, main=mains[1], col=lcol)
contour(x, y, CC$gauss, main=mains[2], col=lcol)
contour(x, y, CC$clayton, main=mains[3], col=lcol)

par(bty="n")
layout(matrix(c(0,1,1,0,  2,2,3,3), nrow=2, byrow=TRUE))
persp(x, y, CC$plackett, zlab=zlab, main=mains[1], col=pcol, theta=tangle, phi=ptangle)
persp(x, y, CC$gauss, zlab=zlab, main=mains[2], col=pcol, theta=tangle, phi=ptangle)
persp(x, y, CC$clayton, zlab=zlab, main=mains[3], col=pcol, theta=tangle, phi=ptangle)

par(oldpar)
}
\keyword{distribution}
\keyword{multivariate}
