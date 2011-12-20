#*** summaryGspline2.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              arnost.komarek[AT]mff.cuni.cz
##
##           CREATED:  04/04/2007
##
## PURPOSE:  Summary for sampled BIVARIATE G-spline
##          
##
## FUNCTIONS: summaryGspline2
##
## ******************************************************************************************

summaryGspline2 <-
  function(x1, x2, mu1, mu2, sigma1, sigma2,
           standard=TRUE, intcpt, scale,
           probs, values=FALSE,
           dir=getwd(), wfile="weight.sim", indfile="knotInd.sim", header=TRUE,
           logw=FALSE, is.indfile=TRUE, 
           skip=0, nwrite)
{
  by <- 1
  
  thispackage <- "glmmAK"
  #thispackage <- NULL
  
  ngrid1 <- length(x1)
  if (!ngrid1) stop("Incorrect x1 supplied")

  ngrid2 <- length(x2)
  if (!ngrid2) stop("Incorrect x2 supplied")

  ngrid <- c(ngrid1, ngrid2)
  total.ngrid <- ngrid1 * ngrid2
  
  nknot1 <- length(mu1)
  if (!nknot1) stop("Incorrect mu1 supplied")

  nknot2 <- length(mu2)
  if (!nknot2) stop("Incorrect mu2 supplied")

  nknot <- c(nknot1, nknot2)
  total.length <- nknot1 * nknot2
  
  if (length(sigma1) == 1) sigma1 <- rep(sigma1, nknot1)
  if (length(sigma1) != nknot1) stop("sigma1 and mu1 must be of the same length")

  if (length(sigma2) == 1) sigma2 <- rep(sigma2, nknot2)
  if (length(sigma2) != nknot2) stop("sigma2 and mu2 must be of the same length")

  if (!("iteration.sim" %in% dir(dir))) stop("File iteration.sim not found")
  iters <- scanFH(paste(dir, "/iteration.sim", sep=""))[,1]
  nsample <- length(iters) - skip
  if (nsample <= 0) stop(paste("Too high value of skip (nsample=", length(iters), ", skip=", skip, ")", sep=""))
      
  if (missing(intcpt)){
    intcpt <- matrix(rep(0, 2*nsample), ncol=2)
  }else{
    if (is.null(dim(intcpt))) stop("intcpt must be either matrix or data.frame")
    if (is.data.frame(intcpt)) intcpt <- as.matrix(intcpt)
    if (ncol(intcpt) != 2) stop("intcpt must have 2 columns")
    if (nrow(intcpt) == length(iters) & skip > 0) intcpt <- intcpt[-(1:skip),]
    if (nrow(intcpt) != nsample) stop("Incorrect intcpt supplied")    
  }  
  
  if (missing(scale)){
    scale <- matrix(rep(1, 2*nsample), ncol=2)
  }else{
    if (is.null(dim(scale))) stop("scale must be either matrix or data.frame")
    if (is.data.frame(scale)) scale <- as.matrix(scale)
    if (ncol(scale) != 2) stop("scale must have 2 columns")
    if (nrow(scale) == length(iters) & skip > 0) scale <- scale[-(1:skip),]
    if (nrow(scale) != nsample) stop("Incorrect scale supplied")
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

  nsummary <- total.ngrid*(1+nprobs)
  nsummary1 <- ngrid1*(1+nprobs)
  nsummary2 <- ngrid2*(1+nprobs)  
  
  if (nprobs | values){
    nvalues <- total.ngrid*nsample
    nvalues1 <- ngrid1*nsample
    nvalues2 <- ngrid2*nsample        
  }else{  
    nvalues <- total.ngrid
    nvalues1 <- ngrid1
    nvalues2 <- ngrid2    
  }
    
  filesindir <- dir(dir)
  if (!length(filesindir)) stop(paste("Directory ", dir, " is empty", sep=""))    
  wpath <- paste(dir, "/", wfile, sep="")
  indpath <- paste(dir, "/", indfile, sep="")

  if (!is.indfile){
    temp <- scan(file=wpath, skip=ifelse(header, 1, 0), nlines=1, quiet=TRUE)
    if (length(temp) != total.length) stop(paste("Number of weights in the 1st row of ", wpath, " is not equal to ", total.length, sep=""))
  }  

  if (header) skip <- skip + 1
  
  fit <- .C("summary_BiGspline",
            value=double(nvalues),      value1=double(nvalues1),            value2=double(nvalues2),
            summary=double(nsummary),   summary1=double(nsummary1),         summary2=double(nsummary2),
            prob=as.double(probs),      nprob=as.integer(nprobs),           retValue=as.integer(values),
            grid1=as.double(x1),        grid2=as.double(x2),                ngrid=as.integer(ngrid),      standard=as.integer(standard),
            knots1=as.double(mu1),      knots2=as.double(mu2),
            sigma1=as.double(sigma1),   sigma2=as.double(sigma2),           nknots=as.integer(nknot),
            intcpt=as.double(intcpt),   tau=as.double(scale),               niter=as.integer(nsample),
            wpath=as.character(wpath),  indpath=as.character(indpath),      skip=as.integer(skip),        by=as.integer(by),
            logw=as.integer(logw),      is.indfile=as.integer(is.indfile),  nwrite=as.integer(nwrite),    err=integer(1),
            PACKAGE = thispackage)
  
  if (fit$err) stop("Something went wrong")
  
  RET <- list(summary=list(x1=x1, x2=x2, Mean=matrix(fit$summary[1:total.ngrid], nrow=ngrid1, ncol=ngrid2)),
              summary1=data.frame(x=x1, Mean=fit$summary1[1:ngrid1]),
              summary2=data.frame(x=x2, Mean=fit$summary2[1:ngrid2]))              
  if (nprobs){
    for (ip in 1:nprobs){
      Quantile <- matrix(fit$summary[(ip*total.ngrid+1):((ip+1)*total.ngrid)], nrow=ngrid1, ncol=ngrid2)
      RET$summary[[paste(probs[ip]*100, "%", sep="")]] <- Quantile
      rm(list="Quantile")
    }  
    
    Quantile1 <- matrix(fit$summary1[-(1:ngrid1)], ncol=nprobs)
    colnames(Quantile1) <- paste(probs*100, "%", sep="")
    RET$summary1 <- cbind(RET$summary1, Quantile1)
    rm(list="Quantile1")    

    Quantile2 <- matrix(fit$summary2[-(1:ngrid2)], ncol=nprobs)
    colnames(Quantile2) <- paste(probs*100, "%", sep="")
    RET$summary2 <- cbind(RET$summary2, Quantile2)
    rm(list="Quantile2")    
  }
  if (values){
    RET$values <- matrix(fit$value, ncol=total.ngrid, byrow=TRUE)
    colnames(RET$values) <- paste(rep(x1, ngrid2), "-", rep(x2, each=ngrid1), sep="")
    rownames(RET$values) <- 1:nsample

    RET$values1 <- matrix(fit$value1, ncol=ngrid1, byrow=TRUE)
    colnames(RET$values1) <- paste(x1)
    rownames(RET$values1) <- 1:nsample

    RET$values2 <- matrix(fit$value2, ncol=ngrid2, byrow=TRUE)
    colnames(RET$values2) <- paste(x2)
    rownames(RET$values2) <- 1:nsample    
  }  

  return(RET)
}
