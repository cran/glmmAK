#*** glmmAK.files2coda.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              akom@email.cz
##
##         CREATED:  02/11/2006 as cumlogitRE.files2coda
##                   15/02/2007 renamed to glmmAK.files2coda
##
## PURPOSE: Summaries based on MCMC
##
#* ********************************************************************************* */

glmmAK.files2coda <- function(dir, drandom=c("none", "normal", "gspline"), quiet=FALSE, skip=0,
                               params=list(prob=FALSE, ecount=FALSE, b=FALSE, alloc=FALSE))
{
  require(coda)
  
  filesindir <- dir(dir)
  if (length(filesindir) == 0) stop(paste("No files available in the directory ", dir, sep=""))

  if (skip < 0) stop("skip must be positive")
  
  if (!length(params)) inparams <- "arnost"
  else                 inparams <- names(params)
  
  RET <- list()
  dopar <- list()
  dopar$loglik <- ifelse(sum(match(filesindir, "loglik.sim", nomatch=0)) == 0, FALSE, TRUE)      
  dopar$betaF <- ifelse(sum(match(filesindir, "betaF.sim", nomatch=0)) == 0, FALSE, TRUE)
  dopar$betaR <- ifelse(sum(match(filesindir, "betaR.sim", nomatch=0)) == 0, FALSE, TRUE)
  dopar$varR <- ifelse(sum(match(filesindir, "varR.sim", nomatch=0)) == 0, FALSE, TRUE)
  dopar$gmoment <- ifelse(sum(match(filesindir, "gmoment.sim", nomatch=0)) == 0, FALSE, TRUE)
  
  dopar$betaRadj <- ifelse(sum(match(filesindir, "betaRadj.sim", nomatch=0)) == 0, FALSE, TRUE)
  dopar$varRadj <- ifelse(sum(match(filesindir, "varRadj.sim", nomatch=0)) == 0, FALSE, TRUE)
  
  dopar$b <- ifelse(sum(match(filesindir, "b.sim", nomatch=0)) == 0, FALSE, TRUE)
  tmp <- match("b", inparams, nomatch=0)
  if (tmp & !params$b) dopar$b <- FALSE

  
  ### Read iteration indeces
  cat("* Reading iteration.sim\n")
  if (!sum(match(filesindir, "iteration.sim", nomatch=0))) stop(paste("iteration.sim not available in the directory ", dir, sep=""))
  help <- read.table(file = paste(dir, "/", "iteration.sim", sep = ""), header = TRUE)
  if (skip) RET$iters <- as.numeric(help[-(1:skip),1])
  else      RET$iters <- as.numeric(help[,1])
  niters <- length(RET$iters)
  start <- RET$iters[1]
  end <- RET$iters[niters]
  thin <- (end - start)/(niters-1)
  
  drandom <- match.arg(drandom)

  ### CHAINS
  ### Log-likelihood
  if (dopar$loglik){
    cat("* Reading loglik.sim\n")    
    loglik <- matrix(scan(file = paste(dir, "/", "loglik.sim", sep = ""), skip=1+skip, quiet=quiet), byrow=TRUE, nrow=niters)
    nloglik <- scan(file = paste(dir, "/", "loglik.sim", sep = ""), nlines=1, what=character(), quiet=quiet)
    colnames(loglik) <- nloglik
    rownames(loglik) <- RET$iters
    RET$loglik <- mcmc(loglik, start=start, end=end, thin=thin)
  }  
    
  ### Fixed effects
  if (dopar$betaF){
    cat("* Reading betaF.sim\n")    
    betaF <- matrix(scan(file = paste(dir, "/", "betaF.sim", sep = ""), skip=1+skip, quiet=quiet), byrow=TRUE, nrow=niters)
    nbetaF <- scan(file = paste(dir, "/", "betaF.sim", sep = ""), nlines=1, what=character(), quiet=quiet)
    colnames(betaF) <- nbetaF
    rownames(betaF) <- RET$iters
    RET$betaF <- mcmc(betaF, start=start, end=end, thin=thin)
  }  

  ### Overall means of random effects
  if (dopar$betaR){
    cat("* Reading betaR.sim\n")        
    betaR <- matrix(scan(file = paste(dir, "/", "betaR.sim", sep = ""), skip=1+skip, quiet=quiet), byrow=TRUE, nrow=niters)
    nbetaR <- scan(file = paste(dir, "/", "betaR.sim", sep = ""), nlines=1, what=character(), quiet=quiet)
    colnames(betaR) <- nbetaR
    rownames(betaR) <- RET$iters
    RET$betaR <- mcmc(betaR, start=start, end=end, thin=thin)
  }  

  ### Adjusted betaR (effect of covariates involved in G-spline random effects)
  if (dopar$betaRadj){
    cat("* Reading betaRadj.sim\n")        
    betaRadj <- matrix(scan(file = paste(dir, "/", "betaRadj.sim", sep = ""), skip=1+skip, quiet=quiet), byrow=TRUE, nrow=niters)
    nbetaRadj <- scan(file = paste(dir, "/", "betaRadj.sim", sep = ""), nlines=1, what=character(), quiet=quiet)
    colnames(betaRadj) <- nbetaRadj
    rownames(betaRadj) <- RET$iters
    RET$betaRadj <- mcmc(betaRadj, start=start, end=end, thin=thin)
  }    
  
  ### Overall variances of random effects
  if (dopar$varR){
    cat("Reading varR.sim\n")        
    varR <- matrix(scan(file = paste(dir, "/", "varR.sim", sep = ""), skip=1+skip, quiet=quiet), byrow=TRUE, nrow=niters)
    nvarR <- scan(file = paste(dir, "/", "varR.sim", sep = ""), nlines=1, what=character(), quiet=quiet)
    colnames(varR) <- nvarR
    rownames(varR) <- RET$iters
    RET$varR <- mcmc(varR, start=start, end=end, thin=thin)
  }  

  ### Variance components of G-spline random effects
  if (dopar$varRadj){
    cat("Reading varRadj.sim\n")        
    varRadj <- matrix(scan(file = paste(dir, "/", "varRadj.sim", sep = ""), skip=1+skip, quiet=quiet), byrow=TRUE, nrow=niters)
    nvarRadj <- scan(file = paste(dir, "/", "varRadj.sim", sep = ""), nlines=1, what=character(), quiet=quiet)
    colnames(varRadj) <- nvarRadj
    rownames(varRadj) <- RET$iters
    RET$varRadj <- mcmc(varRadj, start=start, end=end, thin=thin)
  }  
  
  ### Random effects
  if (dopar$b){
    cat("* Reading b.sim\n")        
    bb <- matrix(scan(file = paste(dir, "/", "b.sim", sep = ""), skip=1+skip, quiet=quiet), byrow=TRUE, nrow=niters)
    nbb <- scan(file = paste(dir, "/", "b.sim", sep = ""), nlines=1, what=character(), quiet=quiet)
    colnames(bb) <- nbb
    rownames(bb) <- RET$iters
    RET$b <- mcmc(bb, start=start, end=end, thin=thin)
  }  
  
  ### CHAINS SPECIFIC FOR MODEL WITH NORMAL RANDOM EFFECTS
  if (drandom == "normal"){

  }  

  ### CHAINS SPECIFIC FOR MODEL WITH G-SPLINE RANDOM EFFECTS
  if (drandom == "gspline"){    
    if (!sum(match(filesindir, "gspline.sim", nomatch=0))) stop(paste("gspline.sim not available in the directory ", dir, sep=""))
    gspline <- read.table(file = paste(dir, "/", "gspline.sim", sep = ""), header = TRUE)    
    
    ### Moments of the G-spline
    if (dopar$gmoment){
      cat("* Reading gmoment.sim\n")                
      gmoment <- matrix(scan(file = paste(dir, "/", "gmoment.sim", sep = ""), skip=1+skip, quiet=quiet), byrow=TRUE, nrow=niters)
      ngmoment <- scan(file = paste(dir, "/", "gmoment.sim", sep = ""), nlines=1, what=character(), quiet=quiet)
      colnames(gmoment) <- ngmoment
      rownames(gmoment) <- RET$iters
      RET$gmoment <- mcmc(gmoment, start=start, end=end, thin=thin)

      cat("* Reading lambda.sim\n")                
      lambda <- matrix(scan(file = paste(dir, "/", "lambda.sim", sep = ""), skip=1+skip, quiet=quiet), byrow=TRUE, nrow=niters)
      nlambda <- scan(file = paste(dir, "/", "lambda.sim", sep = ""), nlines=1, what=character(), quiet=quiet)
      colnames(lambda) <- nlambda
      rownames(lambda) <- RET$iters
      RET$lambda <- mcmc(lambda, start=start, end=end, thin=thin)      
    }            
  }

  return(RET)
}



