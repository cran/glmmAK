#*** QuantileFun.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              arnost.komarek[AT]mff.cuni.cz
##
##           CREATED:  19/01/2007
##
## PURPOSE: Sampled quantiles for sampled functional evaluated in a grid,
##          e.g. sampled quantiles for predictive density (in MCMC context)
##
## FUNCTIONS: QuantileFun
##
## ******************************************************************************************

QuantileFun <- function(x, probs=seq(0, 1, 0.25), vals.in.cols=TRUE)
{
  thispackage <- "glmmAK"
  #thispackage <- NULL
  
  if (is.null(dim(x))){
    namex <- names(x)
    if (is.null(namex)) namex <- paste(1:length(x))
    if (vals.in.cols){
      x <- matrix(x, ncol=1)
      colnames(x) <- "x1"
      rownames(x) <- namex    
    }else{
      x <- matrix(x, nrow=1)
      rownames(x) <- "x1"
      colnames(x) <- namex          
    }  
  }

  if (vals.in.cols) x <- t(x)  
  ngrid <- nrow(x)
  nsample <- ncol(x)  
  nprob <- length(probs)
  
  if (!ngrid | !nsample) stop("Incorrect x supplied")
  if (!nprob) stop("Incorrect prob supplied")
  if (any(is.na(x))) stop("Missing values not allowed in x")
  if (any(probs < 0)) stop("probs must be non-negative")
  if (any(probs > 1)) stop("probs must not be higher than 1")  

  result <- .C("Quantile",
               qs=double(ngrid*nprob),  x=as.double(x),      ngrid=as.integer(ngrid),  nsample=as.integer(nsample),
               prob=as.double(probs),   nprob=as.integer(nprob),
               PACKAGE = thispackage)
  
  result$qs <- as.data.frame(matrix(result$qs, nrow=nprob, ncol=ngrid, byrow=TRUE))
  colnames(result$qs) <- rownames(x)
  rownames(result$qs) <- paste(probs*100, "%", sep="")

  return(result$qs)
}  
