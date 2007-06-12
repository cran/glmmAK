#*** logpoisson.help.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              akom@email.cz
##
##         CREATED:  08/02/2007
##
## PURPOSE: Poisson regression
##          * low level functions
##          * not to be called by the users
##
## FUNCTIONS: design.logpoisson
##            fit.logpoisson
## 
#* ********************************************************************************* */

### To be used inside 'logpoisson'
### -> no input checks
### ------------------------------------
###
### Adjust input design matrices and response vector
### and performs some checks
###
### \item{predict}{if TRUE then arguments \code{y} is ignored. This is intended to be set to TRUE for predict functions.
###          In that case, argument \code{nobs} must indicate number of rows in matrices \code{x}, \code{xb}.
###     }
### \item{nobs}{number of rows in matrices \code{x}, \code{xb}. Only needed if \code{predict=TRUE}, ignored otherwise.}
###
design.logpoisson <- function(y, x, xb, offset, cluster, intcpt.random=FALSE, hierar.center=FALSE, predict=FALSE, nobs)
{
  intercept <- !intcpt.random
  
  if (predict){
    ny <- nobs    
    y <- rep(0, ny)
  } else{
    ny <- length(y)
  }  
  if (!ny) stop("No observations available")

  
  ## Covariates for random effects 
  if (missing(xb) & !intcpt.random){
    xb <- numeric(0)
    nxb <- 0
  }
  else{
    if (intcpt.random){
      if (missing(xb)){
        nxb <- 1
        xb <- matrix(rep(1, ny), ncol=1)
        colnames(xb) <- "(Intercept)"
      }
      else{
        if(is.null(dim(xb))){
          if (length(xb) != ny) stop("length(xb) not consistent with length(y)")
          temp <- deparse(substitute(xb))
          nxb <- 2
          xb <- as.matrix(cbind(rep(1, ny), xb))
          colnames(xb) <- c("(Intercept)", temp)
        }  
        else{
          if (nrow(xb) != ny) stop("nrow(xb) not consistent with length(y)")          
          nxb <- 1 + ncol(xb)
          temp <- colnames(xb)
          xb <- as.matrix(cbind(rep(1, ny), xb))
          colnames(xb) <- c("(Intercept)", temp)
        }  
      }  
    }
    else{    ## else (intcpt.random)
      if(is.null(dim(xb))){
        if (length(xb) != ny) stop("length(xb) not consistent with length(y)")
        temp <- deparse(substitute(xb))        
        nxb <- 1
        xb <- matrix(xb, ncol=1)
        colnames(xb) <- temp
      }
      else{
        if (nrow(xb) != ny) stop("nrow(xb) not consistent with length(y)")        
        nxb <- ncol(xb)
        xb <- as.matrix(xb)
      }  
    }  
  }  ## end of else (missing(xb) & !intcpt.random)


  ## Fixed effects covariates not proportional w.r.t. odds
  if (missing(x) & !intercept){
    x <- numeric(0)
    nx <- 0
  }
  else{
    if (intercept){
      if (missing(x)){
        nx <- 1
        x <- matrix(rep(1, ny), ncol=1)
        colnames(x) <- "(Intercept)"
      }
      else{
        if(is.null(dim(x))){
          if (length(x) != ny) stop("length(x) not consistent with length(y)")
          temp <- deparse(substitute(x))
          nx <- 2
          x <- as.matrix(cbind(rep(1, ny), x))
          colnames(x) <- c("(Intercept)", temp)
        }  
        else{
          if (nrow(x) != ny) stop("nrow(x) not consistent with length(y)")          
          nx <- 1 + ncol(x)
          temp <- colnames(x)
          x <- as.matrix(cbind(rep(1, ny), x))
          colnames(x) <- c("(Intercept)", temp)
        }  
      }  
    }
    else{    ## else (intercept)
      if(is.null(dim(x))){
        if (length(x) != ny) stop("length(x) not consistent with length(y)")
        temp <- deparse(substitute(x))
        nx <- 1
        x <- matrix(x, ncol=1)
        colnames(x) <- temp
      }
      else{
        if (nrow(x) != ny) stop("nrow(x) not consistent with length(y)")        
        nx <- ncol(x)
        x <- as.matrix(x)        
      }  
    }  
  }  ## end of else (missing(x) & !intercept)

  ## Cluster indicator
  if (nxb){
    if (missing(cluster)) stop("cluster indicator must be given")
    if (length(cluster) != ny) stop("length(cluster) not consistent with length(y)")
  }
  else{
    cluster <- 1:ny
  }  

  ## Offset term
  if (missing(offset)){
    offset <- rep(0, ny)
  }
  else{
    if (length(offset) == 1) offset <- rep(offset, ny)
    if (length(offset) != ny) stop("length(offset) not consistent with length(y)")
  }  
  
  
  ## Put everything together
  RESP <- data.frame(y=y, offset=offset, cluster=cluster)

  if (nx){
    RESP  <- cbind(RESP, x)
    bindx <- 1:nx             ## indeces of X in the regression coefficient vector
    indx  <- 3 + (1:nx)       ## indeces of X in RESP
  }
  else{
    bindx <- numeric(0)    
    indx  <- numeric(0)    
  }  
  if (nxb){
    RESP   <- cbind(RESP, xb)
    bindxb <- nx + (1:nxb)           ## indeces of X(b) in the regression coefficient vector    
    indxb  <- 3 + nx + (1:nxb)        ## indeces of X(b) in RESP
  }
  else{
    bindxb <- numeric(0)    
    indxb  <- numeric(0)
  }  

  ## Remove missing observations
  not.NA <- !apply(is.na(RESP), 1, sum)
  RESP <- RESP[not.NA,]

  ## Sort it with respect to the cluster
  RESP <- RESP[order(RESP$cluster),]

  ## Extract various components
  y       <- as.numeric(RESP[, 1])
  offset  <- as.numeric(RESP[, 2])
  cluster <- as.numeric(RESP[, 3])
  ni <- as.numeric(table(cluster))
  N <- length(ni)

  if (nx){
    x <- as.matrix(RESP[, indx])
    colnames(x) <- colnames(RESP)[indx]
  }
  if (nxb){
    xb <- as.matrix(RESP[, indxb])
    colnames(xb) <- colnames(RESP)[indxb]
  }    
  ny <- length(y)

  ## Number of observations and other checks
  if (!(nx+nxb)) stop("Neither covariates nor intercept included in the model")

  if (any(y < 0)) stop("Response may not be negative")

  xv.qr <- qr(RESP[,-(1:3)])
  if (xv.qr$rank < nx+nxb & !predict){
    warning(paste("Rank of the matrix (X,X(b)) is ", xv.qr$rank, " < ", nx+nxb, sep=""))
  }  

  if (!hierar.center){
    ## Design matrices for fixed effects -> put together X, X(b)
    if (nx & nxb){
      x <- cbind(x, xb)
      nx <- nx + nxb
    }  
    else{
      if (nxb){        
        x <- xb
        nx <- nxb
      }
    }
  }  
  
  RET <- list(y=y, offset=offset, cluster=cluster, x=x, xb=xb,
              ny=ny, N=N, ni=ni, nx=nx, nxb=nxb,
              bindx=bindx, bindxb=bindxb)
  return(RET)  
}  


### ******************************************************************

fit.logpoisson <- function(y, x, offset, epsilon=1e-08, maxit=25, trace=FALSE)
{
  thispackage <- "glmmAK"
  #thispackage <- NULL

  ny <- length(y)
  if (length(x)){
    nx <- ifelse(is.null(dim(x)), 1, ncol(x))
  }
  else{
    nx <- 0
  }  
  
  ## Initial values
  intcpt.init <- mean(log(y+0.1) - offset)
  if (nx == 1)  beta.init <- intcpt.init
  else{
    if ("(Intercept)" %in% colnames(x)){      
      beta.init <- rep(0, nx)
      names(beta.init) <- colnames(x)
      beta.init["(Intercept)"] <- intcpt.init      
    }
    else{
      beta.init <- c(intcpt.init, rep(0, nx-1))
    }
  }  
  names(beta.init) <- colnames(x)

  ## Fit the model
  LTnx <- nx*(nx+1)/2
  fit <- .C("fit_poisson",
            theta=as.double(beta.init),    eta=double(ny),              mu=double(ny),
            ll=double(1),                  U=double(nx),                I=double(LTnx),
            Y=as.integer(y),               offset=as.double(offset),    X=as.double(t(x)),
            p=as.integer(nx),              n=as.integer(ny),            niter=as.integer(maxit),
            toler=as.double(epsilon),      trace=as.integer(trace),     err=integer(1),
            PACKAGE=thispackage)

  Imat <- matrix(NA, nrow=nx, ncol=nx)
  Imat[lower.tri(Imat, diag=TRUE)] <- fit$I
  tImat <- t(Imat)
  Imat[upper.tri(Imat, diag=FALSE)] <- tImat[upper.tri(tImat, diag=FALSE)]
  fit$I <- Imat
  if (fit$err)
    if (fit$err == 1) warning("Convergence not reached. Maximal number of step-halfing steps performed.")
    else
      if (fit$err == 2) warning("Convergence not reached. Maximal number of iterations performed.")
      else
        if (fit$err == 3) warning("Convergence not reached. Minus Hessian not positive definite.")
        else
          if (fit$err == 100) warning("Convergence criterion satisfied but minus Hessian is not positive definite.")
          else
            if (fit$err == 101) warning("Convergence not reached. Maximal number of step-halfing steps performed. Additionally: minus Hessian is not positive definite.")
            else            
              if (fit$err == 102) warning("Convergence not reached. Maximal number of iterations performed. Additionally: minus Hessian is not positive definite.")
                else{
                  cat("Error ", fit$err, "\n", sep="")
                  stop("Numerical or other problems. No results produced.")
                }  

  names(fit$theta) <- names(fit$U) <- colnames(fit$I) <- rownames(fit$I) <- names(beta.init)
  

  ### NULL MODEL (intercept only)
  x.NULL <- matrix(rep(1, ny), ncol=1)
  colnames(x.NULL) <- "(Intercept)"    
  fit.NULL <- .C("fit_poisson",
                 theta=as.double(intcpt.init),  eta=double(ny),                mu=double(ny),
                 ll=double(1),                  U=double(1),                   I=double(1),
                 Y=as.integer(y),               offset=as.double(offset),      X=as.double(t(x.NULL)),
                 p=as.integer(1),               n=as.integer(ny),              niter=as.integer(maxit),
                 toler=as.double(epsilon),      trace=as.integer(trace),       err=integer(1),
                 PACKAGE=thispackage)

  LogLik <- c(fit$ll, fit.NULL$ll)
  names(LogLik) <- c("model", "null")

  RET <- list(coefficients=fit$theta,
              loglik=LogLik,
              score=fit$U,
              vcov=fit$I,
              linear.predictors=fit$eta,
              fitted.values=fit$mu,
              converged=(fit$err==0),
              iter=fit$niter,
              y=y,
              x=x)

  class(RET) <- "logpoisson"
  return(RET)    
}

