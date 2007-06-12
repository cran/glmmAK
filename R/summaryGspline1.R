#*** summaryGspline1.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              akom@email.cz
##
##           CREATED:  19/01/2007
##
## PURPOSE:  Summary for sampled UNIVARIATE G-spline
##          
##
## FUNCTIONS: summaryGspline1
##
## ******************************************************************************************

summaryGspline1 <-
  function(x, mu, sigma,
           standard=TRUE, intcpt, scale,
           probs, values=FALSE,
           dir=getwd(), wfile="logweight.sim", header=TRUE, logw=TRUE,
           skip=0, nwrite)
{
  by <- 1
  
  thispackage <- "glmmAK"
  #thispackage <- NULL

  ngrid <- length(x)
  if (!ngrid) stop("Incorrect x supplied")
  
  nknot <- length(mu)
  if (!nknot) stop("Incorrect mu supplied")

  if (length(sigma) == 1) sigma <- rep(sigma, nknot)
  if (length(sigma) != nknot) stop("sigma and mu must be of the same length")

  if (!("iteration.sim" %in% dir(dir))) stop("File iteration.sim not found")
  iters <- scanFH(paste(dir, "/iteration.sim", sep=""))[,1]
  nsample <- length(iters) - skip
  if (nsample <= 0) stop(paste("Too high value of skip (nsample=", length(iters), ", skip=", skip, ")", sep=""))

  if (missing(intcpt)){
    intcpt <- rep(0, nsample)
  }else{
    if (length(intcpt) == length(iters) & skip > 0) intcpt <- intcpt[-(1:skip)]
    if (length(intcpt) != nsample) stop("Incorrect intcpt supplied")    
  }
  
  if (missing(scale)){
    scale <- rep(1, nsample)
  }else{
    if (length(scale) == length(iters) & skip > 0) scale <- scale[-(1:skip)]
    if (length(scale) != nsample) stop("Incorrect scale supplied")
    if (any(scale <= 0)) stop("scale values must be positive")    
  }  

  if (skip <= 0) skip <- 0
  if (by > nsample)  stop(paste("by must not be higher than ", nsample, sep=""))
  if (missing(nwrite)) nwrite <- nsample

  if (missing(probs)){
    probs <- 0
    nprobs <- 0
  }else{
    if (any(probs < 0)) stop("probs must be >= 0")
    if (any(probs > 1)) stop("probs must be <= 1")    
    nprobs <- length(probs)
  }

  nsummary <- ngrid*(1+nprobs)
  
  if (nprobs | values) nvalues <- ngrid*nsample
  else                 nvalues <- ngrid

  filesindir <- dir(dir)
  if (!length(filesindir)) stop(paste("Directory ", dir, " is empty", sep=""))    
  wpath <- paste(dir, "/", wfile, sep="")
  temp <- scan(file=wpath, skip=ifelse(header, 1, 0), nlines=1, quiet=TRUE)
  if (length(temp) != nknot) stop(paste("Number of weights in the 1st row of ", wpath, " is not equal to ", nknot, sep=""))

  if (header) skip <- skip + 1
  
  fit <- .C("summary_Gspline1",
            value=double(nvalues),      summary=double(nsummary),
            prob=as.double(probs),      nprob=as.integer(nprobs),   retValue=as.integer(values),
            grid=as.double(x),          ngrid=as.integer(ngrid),    standard=as.integer(standard),
            knots=as.double(mu),        sigma=as.double(sigma),     nknots=as.integer(nknot),
            intcpt=as.double(intcpt),   tau=as.double(scale),       niter=as.integer(nsample),
            wpath=as.character(wpath),  skip=as.integer(skip),      by=as.integer(by),
            logw=as.integer(logw),      nwrite=as.integer(nwrite),  err=integer(1),
            PACKAGE = thispackage)

  if (fit$err) stop("Something went wrong")
  
  RET <- list(summary=data.frame(x=x, Mean=fit$summary[1:ngrid]))
  if (nprobs){
    Quantile <- matrix(fit$summary[-(1:ngrid)], ncol=nprobs)
    colnames(Quantile) <- paste(probs*100, "%", sep="")
    RET$summary <- cbind(RET$summary, Quantile)
  }
  if (values){
    RET$values <- matrix(fit$value, ncol=ngrid, byrow=TRUE)
    colnames(RET$values) <- paste(x)
    rownames(RET$values) <- 1:nsample
  }  

  return(RET)
}
