#*** gspline1.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              arnost.komarek[AT]mff.cuni.cz
##
##         CREATED:  22/01/2007
##
## PURPOSE: Values of the density based on the UNIVARIATE G-spline
##          + random numbers generation
##
## FUNCTIONS:  22/01/2007:  dgspline1
##             29/01/2007:  rgspline1
##            
## 
#* ********************************************************************************* */


## Random numbers generation from the UNIVARIATE G-spline
rgspline1 <- function(n, mu, sigma, weight, intcpt=0, scale=1, logw=TRUE)
{
  thispackage <- "glmmAK"
  #thispackage <- NULL
  
  if (n <= 0) stop("Incorrect n supplied")
  
  nknot <- length(mu)
  if (!nknot) stop("Incorrect mu supplied")

  if (length(sigma) == 1) sigma <- rep(sigma, nknot)
  if (length(sigma) != nknot) stop("sigma and mu must be of the same length")
  
  if (length(weight) != nknot) stop("Incorrect weight supplied")

  if (length(intcpt) != 1) stop("Incorrect intcpt supplied")
  
  if (length(scale) != 1) stop("Incorrect scale supplied")  
  if (scale <= 0) stop("Value of scale must be positive")

  if (logw) weight <- exp(weight)
  weight <- weight/sum(weight)
  if (any(weight < 0)) stop("weights may not be negative")

  sample <- .C("rGspline1R",
                      x=double(n),               weight=as.double(weight),   n=as.integer(n),
                      knots=as.double(mu),       sigma=as.double(sigma),     nknots=as.integer(nknot),
                      intcpt=as.double(intcpt),  tau=as.double(scale),       logw=as.integer(0),
     PACKAGE = thispackage)
  
  return(sample$x)   
}  
  

## Values of the density based on the UNIVARIATE G-spline
dgspline1 <- function(x, mu, sigma, weight, intcpt=0, scale=1, logw=TRUE)
{
  thispackage <- "glmmAK"
  #thispackage <- NULL
  
  ngrid <- length(x)
  if (!ngrid) stop("Incorrect x supplied")
  
  nknot <- length(mu)
  if (!nknot) stop("Incorrect mu supplied")

  if (length(sigma) == 1) sigma <- rep(sigma, nknot)
  if (length(sigma) != nknot) stop("sigma and mu must be of the same length")
  
  if (length(weight) != nknot) stop("Incorrect weight supplied")

  if (length(intcpt) != 1) stop("Incorrect intcpt supplied")
  
  if (length(scale) != 1) stop("Incorrect scale supplied")  
  if (scale <= 0) stop("Value of scale must be positive")

  if (logw) weight <- exp(weight)
  weight <- weight/sum(weight)
  if (any(weight < 0)) stop("weights may not be negative")

#  gx <- x
#  for (i in 1:length(x)){
#    z <- (x[i] - intcpt)/scale
#    gkx <- dnorm(z, mean=mu, sd=sigma)
#    gx[i] <- sum(weight*gkx)/scale
#  }

  fit <- .C("eval_Gspline1",
                      average=double(ngrid),     value=double(ngrid),
                      weight=as.double(weight),  knots.tau=double(nknot),  sigma.tau=double(nknot),
                      grid=as.double(x),         ngrid=as.integer(ngrid),  standard=as.integer(0),
                      knots=as.double(mu),       sigma=as.double(sigma),   nknots=as.integer(nknot),
                      intcpt=as.double(intcpt),  tau=as.double(scale),     logw=as.integer(0),
     PACKAGE = thispackage)

  if (!is.null(dim(x))){
    attr(fit$value, "dim") <- attr(x, "dim")
    attr(fit$value, "dimnames") <- attr(x, "dimnames")
  }else{
    attr(fit$value, "names") <- attr(x, "names")
  }  

  return(fit$value)  
}  
