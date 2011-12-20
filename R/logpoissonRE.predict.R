#*** logpoissonRE.predict.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              arnost.komarek[AT]mff.cuni.cz
##
##                   CREATED:  16/02/2007
##
## PURPOSE: Prediction for the Pooisson log-linear model with random effects
##          fitted using Bayesian specification and MCMC
##
## FUNCTIONS: logpoissonRE.predict
## 
#* ********************************************************************************* */

logpoissonRE.predict <- function(nobs, x, xb, offset, cluster,
      intcpt.random=FALSE, hierar.center=FALSE,     
      drandom=c("normal", "gspline"),
      betaF, betaR, varR, is.varR=TRUE,
      prior.gspline,
      probs, values=FALSE,
      dir=getwd(), wfile, indfile, header=TRUE, logw, is.indfile,
      skip=0, nwrite)                               
{
  by <- 1
  
  thispackage <- "glmmAK"
  thispackage <- NULL

  ## Design matrices and check of some input arguments
  y <- rep(0, nobs)  
  DES <- design.logpoisson(y=y, x=x, xb=xb, offset=offset, cluster=cluster, intcpt.random=intcpt.random, hierar.center=hierar.center, predict=TRUE, nobs=nobs)
  if (hierar.center){
    nFixed <- DES$nx
    nRandom <- DES$nxb
  }
  else{
    nFixed <- DES$nx              ## nxb is included in nx, nvb is included in nv
    nRandom <- DES$nxb
  }  

  ## nsample from iteration.sim, nwrite
  if (skip <= 0) skip <- 0  
  if (!("iteration.sim" %in% dir(dir))) stop("File iteration.sim not found")
  iters <- scanFH(paste(dir, "/iteration.sim", sep=""))[,1]
  nsample <- length(iters) - skip
  if (nsample <= 0) stop(paste("Too high value of skip (nsample=", length(iters), ", skip=", skip, ")", sep=""))
  if (missing(nwrite)) nwrite <- nsample

  ## probs and sizes of a needed space
  if (missing(probs)){
    probs <- 0
    nprobs <- 0
  }else{
    if (any(probs < 0)) stop("probs must be >= 0")
    if (any(probs > 1)) stop("probs must be <= 1")    
    nprobs <- length(probs)
  }

  nsummary <- nobs*(1+nprobs)
  
  if (nprobs | values) nvalues <- nobs*nsample
  else                 nvalues <- nobs
  
  ## Check sampled values of betaF, betaR and transpose them
  if (nFixed){
    if (missing(betaF)) stop(paste("nFixed=", nFixed, ", betaF must be given", sep=""))
    if (nFixed == 1){
      if (is.null(dim(betaF))){
        if (length(betaF) == length(iters) & skip > 0) betaF <- betaF[-(1:skip)]
        if (length(betaF) != nsample) stop(paste("nsample=", nsample, ",  length(betaF)=", length(betaF), " is inconsistent", sep=""))        
      }else{
        if (nrow(betaF) == length(iters) & skip > 0) betaF <- betaF[-(1:skip),]
        if (nrow(betaF) != nsample) stop(paste("nsample=", nsample, ",  length(betaF)=", nrow(betaF), " is inconsistent", sep=""))
        betaF <- as.numeric(betaF[,1])
      }  
    }else{
      if (is.null(dim(betaF))) stop(paste("nFixed=", nFixed, ", betaF must be either matrix or data.frame", sep=""))
      if (nrow(betaF) == length(iters) & skip > 0) betaF <- betaF[-(1:skip),]      
      if (nrow(betaF) != nsample) stop(paste("nsample=", nsample, ",  nrow(betaF)=", nrow(betaF), " is inconsistent", sep=""))
      betaF <- as.numeric(t(betaF))
    }  
  }else{
    betaF <- 0
  }  

  if (nRandom & hierar.center){
    if (missing(betaR)) stop(paste("nRandom=", nRandom, " and hierar.center=TRUE, betaR must be given", sep=""))
    if (nRandom == 1){
      if (is.null(dim(betaR))){
        if (length(betaR) == length(iters) & skip > 0) betaR <- betaR[-(1:skip)]        
        if (length(betaR) != nsample) stop(paste("nsample=", nsample, ",  length(betaR)=", length(betaR), " is inconsistent", sep=""))        
      } else{
        if (nrow(betaR) == length(iters) & skip > 0) betaR <- betaR[-(1:skip),]        
        if (nrow(betaR) != nsample) stop(paste("nsample=", nsample, ",  length(betaR)=", nrow(betaR), " is inconsistent", sep=""))
        betaR <- as.numeric(betaR[,1])        
      }  
    }else{
      if (is.null(dim(betaR))) stop(paste("nRandom=", nRandom, ", betaR must be either matrix or data.frame", sep=""))
      if (nrow(betaR) == length(iters) & skip > 0) betaR <- betaR[-(1:skip),]      
      if (nrow(betaR) != nsample) stop(paste("nsample=", nsample, ",  nrow(betaR)=", nrow(betaR), " is inconsistent", sep=""))
      betaR <- as.numeric(t(betaR))
    }  
  }else{
    betaR <- 0
  }  
      
  ## Distribution of random effects and check sampled values of varR (which must be transposed)
  if (nRandom){
    if (missing(varR)) stop(paste("nRandom=", nRandom, ", varR must be given", sep=""))
    
    drandom <- match.arg(drandom)
    if (drandom == "normal"){
      drandomI <- 1

      if (nRandom == 1){
        if (is.null(dim(varR))){
          if (length(varR) == length(iters) & skip > 0) varR <- varR[-(1:skip)]        
          if (length(varR) != nsample) stop(paste("nsample=", nsample, ",  length(varR)=", length(varR), " is inconsistent", sep=""))
        }else{
          if (nrow(varR) == length(iters) & skip > 0) varR <- varR[-(1:skip),]         
          if (nrow(varR) != nsample) stop(paste("nsample=", nsample, ",  length(varR)=", nrow(varR), " is inconsistent", sep=""))
          varR <- as.numeric(varR[,1])          
        }  
      }else{
        lvarR <- 0.5*nRandom*(1+nRandom)
        if (is.null(dim(varR))) stop(paste("nRandom=", nRandom, ", varR must be either matrix or data.frame", sep=""))
        if (nrow(varR) == length(iters) & skip > 0) varR <- varR[-(1:skip),]        
        if (nrow(varR) != nsample) stop(paste("nsample=", nsample, ",  nrow(varR)=", nrow(varR), " is inconsistent", sep=""))
        if (ncol(varR) != lvarR) stop(paste("nrandom=", nRandom, ",  ncol(varR)=", ncol(varR), " is inconsistent", sep=""))
        varR <- as.numeric(t(varR))
      }        
    }else{
      if (drandom == "gspline"){
        drandomI <- 2

        if (nRandom == 1){
          if (is.null(dim(varR))){
            if (length(varR) == length(iters) & skip > 0) varR <- varR[-(1:skip)]        
            if (length(varR) != nsample) stop(paste("nsample=", nsample, ",  length(varR)=", length(varR), " is inconsistent", sep=""))
          }else{
            if (nrow(varR) == length(iters) & skip > 0) varR <- varR[-(1:skip),]         
            if (nrow(varR) != nsample) stop(paste("nsample=", nsample, ",  length(varR)=", nrow(varR), " is inconsistent", sep=""))
            varR <- as.numeric(varR[,1])          
          }  
        }else{
          if (nRandom == 2){
            lvarR <- nRandom
            if (is.null(dim(varR))) stop(paste("nRandom=", nRandom, ", varR must be either matrix or data.frame", sep=""))
            if (nrow(varR) == length(iters) & skip > 0) varR <- varR[-(1:skip),]        
            if (nrow(varR) != nsample) stop(paste("nsample=", nsample, ",  nrow(varR)=", nrow(varR), " is inconsistent", sep=""))
            if (ncol(varR) != lvarR) stop(paste("nrandom=", nRandom, ",  ncol(varR)=", ncol(varR), " is inconsistent", sep=""))
            varR <- as.numeric(t(varR))
          }else{            
            stop("Not (yet) implemented for multivariate G-splines")
          }   
        }
      }  
    }
  }else{
    drandom <- "none"
    drandomI <- 0
    varR <- 0
  }  

  ## Prior for G-spline and check sampled weights
  prior.gspline <- prior.gspline.glmmAK(prior.gspline=prior.gspline, nRandom=nRandom, drandom=drandom, simplified=TRUE)
  Gdim <- attr(prior.gspline, "Dim")               ### nRandom, K[1],...,K[nRandom]
  GdPar <- attr(prior.gspline, "dPar")             ### sigma[1],...,sigma[nRandom], knots

  if (drandom == "gspline"){
    if (nRandom == 1){
      if (missing(logw)) logw <- TRUE
      if (missing(wfile)){
        if (logw) wfile <- "logweight.sim"
        else      wfile <- "weight.sim"
      }  
      indpath <- paste(dir, "/", "knotInd.sim", sep="")
      is.indfile <- FALSE
      
      nknot <- 2*Gdim[2] + 1
      
      filesindir <- dir(dir)
      if (!length(filesindir)) stop(paste("Directory ", dir, " is empty", sep=""))    
      wpath <- paste(dir, "/", wfile, sep="")
      temp <- scan(file=wpath, skip=ifelse(header, 1, 0), nlines=1, quiet=TRUE)
      if (length(temp) != nknot) stop(paste("Number of weights in the 1st row of ", wpath, " is not equal to ", nknot, sep=""))
    }else{
      if (nRandom == 2){
        if (missing(logw)) logw <- FALSE
        if (missing(is.indfile)){
          if (logw) is.indfile <- FALSE
          else      is.indfile <- TRUE
        }          
        if (missing(wfile)){
          if (is.indfile & !logw) wfile <- "weight.sim"
          else if (!is.indfile & logw) wfile <- "logweight.sim"
               else stop("wfile must be given")
        }  
        if (missing(indfile)) indfile <- "knotInd.sim"

        nknot <- 2*Gdim[2:3] + 1
        total.length <- nknot[1]*nknot[2]
        
        filesindir <- dir(dir)
        if (!length(filesindir)) stop(paste("Directory ", dir, " is empty", sep=""))    
        wpath <- paste(dir, "/", wfile, sep="")
        indpath <- paste(dir, "/", indfile, sep="")        

        if (is.indfile){
          temp <- scan(file=indpath, skip=ifelse(header, 1, 0), nlines=1, quiet=TRUE)
          ktemp <- temp[1]
          if (ktemp <= 0) stop(paste("First data value in ", indpath, " is not positive", sep=""))
          if (ktemp > total.length) stop(paste("First data value in ", indpath, " is higher than ", total.length, sep=""))

          temp <- scan(file=wpath, skip=ifelse(header, 1, 0), nlines=1, quiet=TRUE)
          if (length(temp) != ktemp) stop(paste("Number of weights in the 1st row of ", wpath, " is not equal to ", ktemp, sep=""))          
        }else{
          temp <- scan(file=wpath, skip=ifelse(header, 1, 0), nlines=1, quiet=TRUE)
          if (length(temp) != total.length) stop(paste("Number of weights in the 1st row of ", wpath, " is not equal to ", total.length, sep=""))
        }          
      }else{  
        stop("Not (yet) implemented for G-splines of dimension > 2")
      }
    }  
  }else{
    logw <- TRUE
    wpath <- paste(dir, "/", "logweight.sim", sep="")
    indpath <- paste(dir, "/", "knotInd.sim", sep="")
    is.indfile <- FALSE
  }  

  if (header) skip <- skip + 1
  
  fit <- .C("predict_poisson",
            value=double(nvalues),         summary=double(nsummary),
            qprob=as.double(probs),        nqprob=as.integer(nprobs),               retValue=as.integer(values),
            n=as.integer(nobs),            N=as.integer(DES$N),                     ni=as.integer(DES$ni),
            offset=as.double(DES$offset),  X=as.double(t(DES$x)),                   p=as.integer(DES$nx),
            XRE=as.double(t(DES$xb)),      pRE=as.integer(DES$nxb),
            REdist=as.integer(drandomI),   hierarCenter=as.integer(hierar.center),
            Gdim=as.integer(Gdim),         GdPar=as.double(GdPar),
            betaF=as.double(betaF),        betaR=as.double(betaR),                  varR=as.double(varR),         is.varR=as.integer(is.varR),
            niter=as.integer(nsample),
            wpath=as.character(wpath),     indpath=as.character(indpath),           skip=as.integer(skip),        by=as.integer(by),
            logw=as.integer(logw),         is.indfile=as.integer(is.indfile),       nwrite=as.integer(nwrite),    err=integer(1),
            PACKAGE = thispackage)

  if (fit$err) stop("Something went wrong")

  colnaam <- "ecount"
  rownaam <- 1:nobs
  lsumm   <- nobs
  
  RET <- list(Mean=matrix(fit$summary[1:lsumm], ncol=1, nrow=nobs, byrow=TRUE))
  colnames(RET$Mean) <- colnaam
  rownames(RET$Mean) <- rownaam
  if (nprobs){
    for (i in 1:nprobs){
      RET[[i+1]] <- matrix(fit$summary[(i*lsumm+1):((i+1)*lsumm)], ncol=1, nrow=nobs, byrow=TRUE)
      colnames(RET[[i+1]]) <- colnaam
      rownames(RET[[i+1]]) <- rownaam    
    }
    names(RET)[2:(nprobs+1)] <- paste(probs*100, "%", sep="")
  }
  if (values){
    RET$values <- matrix(fit$value, ncol=lsumm, byrow=TRUE)
    colnames(RET$values) <- paste(colnaam, ":", rep(1:nobs, each=1), sep="")
    rownames(RET$values) <- 1:nsample
  }   
  
  return(RET)    
}
