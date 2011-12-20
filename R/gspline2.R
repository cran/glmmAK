#*** gspline2.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              arnost.komarek[AT]mff.cuni.cz
##
##         CREATED:  12/04/2007
##
## PURPOSE: Values of the density based on the BIVARIATE G-spline
##          + random numbers generation
##
## FUNCTIONS:  ??/??/200?:  dgspline2 (not yet implemented)
##             12/04/2007:  rgspline2
##            
## 
#* ********************************************************************************* */


## Random numbers generation from the BIVARIATE G-spline
rgspline2 <- function(n, mu1, mu2, sigma1, sigma2, weight, knotInd, intcpt=0, scale=1, logw=TRUE)
{
  thispackage <- "glmmAK"
  #thispackage <- NULL
  
  if (n <= 0) stop("Incorrect n supplied")
  
  nknot1 <- length(mu1)
  if (!nknot1) stop("Incorrect mu1 supplied")

  nknot2 <- length(mu2)
  if (!nknot2) stop("Incorrect mu2 supplied")

  total.length <- nknot1 * nknot2
  
  if (length(sigma1) == 1) sigma1 <- rep(sigma1, nknot1)
  if (length(sigma1) != nknot1) stop("sigma1 and mu1 must be of the same length")

  if (length(sigma2) == 1) sigma2 <- rep(sigma2, nknot2)
  if (length(sigma2) != nknot2) stop("sigma2 and mu2 must be of the same length")

  if (missing(knotInd)){
    is.indfile <- 0
    k.effect <- total.length
    indWeight <- rep(1, total.length)    
  }else{
    is.indfile <- 1
    k.effect <- length(knotInd)
    if (!k.effect) stop("Incorrect indWeight supplied")
    if (any(knotInd < 0) | any(knotInd >= total.length)) stop("indWeight out of the range")
    if (length(weight) != k.effect) stop("Incorrect weight supplied")
    
    tempw <- weight
    weight <- rep(0, total.length)
    weight[knotInd+1] <- tempw

    indWeight <- rep(0, total.length)
    indWeight[knotInd+1] <- 1
  }  
  if (length(weight) != total.length) stop("Incorrect weight supplied")    

  if (length(intcpt) == 1) intcpt <- rep(intcpt, 2)
  if (length(intcpt) != 2) stop("Incorrect intcpt supplied")

  if (length(scale) == 1) scale <- rep(scale, 2)  
  if (length(scale) != 2) stop("Incorrect scale supplied")  
  if (any(scale <= 0)) stop("Value of scale must be positive")

  if (!logw) if (any(weight < 0)) stop("weights may not be negative")
  
  sample <- .C("rBiGsplineR",
                      x=double(2*n),             weight=as.double(weight),   ind.w.effect=as.integer(indWeight),
                      n=as.integer(n),
                      knots0=as.double(mu1),     knots1=as.double(mu2),
                      sigma0=as.double(sigma1),  sigma1=as.double(sigma2),   nknots=as.integer(c(nknot1, nknot2)),    total.length=as.integer(total.length),
                      intcpt=as.double(intcpt),  tau=as.double(scale),
                      logw=as.integer(logw),     is.indfile=as.integer(is.indfile),
     PACKAGE = thispackage)

  #cat("weight=       ", weight, "\n", sep="  ")  
  #cat("sample$weight=", sample$weight, "\n", sep="  ")
  #cat("indWeight=    ", indWeight, "\n", sep="  ")  
  #cat("ind.w.effect= ", sample$ind.w.effect, "\n", sep="  ")  

  return(matrix(sample$x, nrow=n, ncol=2, byrow=TRUE))  
}  
