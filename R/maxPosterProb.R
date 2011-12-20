#*** maxPosterProb.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              arnost.komarek[AT]mff.cuni.cz
##
##         CREATED:  20/10/2006
##
## PURPOSE: Some G-spline utilities
##
## FUNCTIONS: maxPosterProb
## 
#* ********************************************************************************* */


##
## Wrapper to C++ max_poster_prob function from util_Gspline.cpp
##
maxPosterProb <- function(data, intercept, std.dev, K, delta, sigma)
{
  thispackage <- "glmmAK"
  #thispackage <- NULL
  
  if (is.null(dim(data))){
    Dim <- 1
    nData <- length(data)
  }
  else{
    Dim <- ncol(data)
    nData <- nrow(data)
  }
  if (!Dim | !nData) stop("No data supplied")

  if (length(intercept) == 1) intercept <- rep(intercept, Dim)
  if (length(intercept) != Dim) stop("Incorrect intercept supplied")

  if (length(std.dev) == 1) std.dev <- rep(std.dev, Dim)  
  if (length(std.dev) != Dim) stop("Incorrect std.dev supplied")
  if (any(std.dev <= 0)) stop("std.dev must be all positive")

  if (length(K) == 1) K <- rep(K, Dim)
  if (length(K) != Dim) stop("Incorrect K supplied")  
  if (any(K < 0)) stop("K must be all non-negative")

  if (length(delta) == 1) delta <- rep(delta, Dim)
  if (length(delta) != Dim) stop("Incorrect delta supplied")  
  if (any(delta <= 0)) stop("K must be all positive")
  
  if (length(sigma) == 1) sigma <- rep(sigma, Dim)  
  if (length(sigma) != Dim) stop("Incorrect sigma supplied")  
  if (any(sigma <= 0)) stop("sigma must be all positive")
  
  knots <- numeric()
  for (i in 1:Dim){
    Knots <- ((-K[i]):K[i])*delta[i]
    names(Knots) <- paste("mu", i, ".", (-K[i]):K[i], sep="")
    knots <- c(knots, Knots)
  }  

  fit <- .C("max_poster_prob", alloc=integer(Dim*nData),
                               sigma=as.double(sigma),
                               knots=as.double(knots),
                               K=as.integer(K),
                               data=as.double(t(data)),
                               Dim=as.integer(Dim),
                               nData=as.integer(nData),
                               intercept=as.double(intercept),
                               std.dev=as.double(std.dev),
            PACKAGE=thispackage
            )

  if (Dim == 1){
    Alloc <- fit$alloc
    names(Alloc) <- names(data)
  }
  else{
    Alloc <- matrix(fit$alloc, nrow=nData, ncol=Dim, byrow=TRUE)
    rownames(Alloc) <- rownames(data)
    colnames(Alloc) <- colnames(data)    
  }  

  return(Alloc)  
}  
