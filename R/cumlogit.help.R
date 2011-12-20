#*** cumlogit.help.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              arnost.komarek[AT]mff.cuni.cz
##
##         CREATED:  25/09/2006
##
## PURPOSE: Cumulative logit model
##          * low level functions
##          * not to be called by the users
##
## FUNCTIONS: design.cumlogit
##            fit.cumlogit
##            init00.cumlogit
##            prob00.cumlogit
##            linear.predictors00.cumlogit
## 
#* ********************************************************************************* */


### To be used inside 'cumlogit'
### -> no input checks
### ------------------------------------
###
### Adjust input design matrices and response vector
### and performs some checks
###
### \item{predict}{if TRUE then arguments \code{y} is ignored. This is intended to be set to TRUE for predict functions.
###          In that case, argument \code{nobs} must indicate number of rows in matrices \code{v}, \code{x}, \code{vb}, \code{xb}.
###     }
### \item{nobs}{number of rows in matrices \code{v}, \code{x}, \code{vb}, \code{xb}. Only needed if \code{predict=TRUE}, ignored otherwise.}
###
design.cumlogit <- function(y, v, x, vb, xb, cluster, intcpt.random=FALSE, hierar.center=FALSE, C, predict=FALSE, nobs)
{
  intercept <- !intcpt.random
  
  if (predict){
    ny <- nobs    
    y <- rep(0, ny)
  } else{
    ny <- length(y)
  }  
  if (!ny) stop("No observations available")

  if (C <= 0) stop("Incorrect 'C' argument supplied")
  
  ## Covariates for random effects, not proportional w.r.t. odds
  if (missing(vb) & !intcpt.random){
    vb <- numeric(0)
    nvb <- 0
  }
  else{
    if (intcpt.random){
      if (missing(vb)){
        nvb <- 1
        vb <- matrix(rep(1, ny), ncol=1)
        colnames(vb) <- "(Intercept)"
      }
      else{
        if(is.null(dim(vb))){
          if (length(vb) != ny) stop("length(vb) not consistent with length(y)")
          temp <- deparse(substitute(vb))
          nvb <- 2
          vb <- as.matrix(cbind(rep(1, ny), vb))
          colnames(vb) <- c("(Intercept)", temp)
        }  
        else{
          if (nrow(vb) != ny) stop("nrow(vb) not consistent with length(y)")          
          nvb <- 1 + ncol(vb)
          temp <- colnames(vb)
          vb <- as.matrix(cbind(rep(1, ny), vb))
          colnames(vb) <- c("(Intercept)", temp)
        }  
      }  
    }
    else{    ## else (intcpt.random)
      if(is.null(dim(vb))){
        if (length(vb) != ny) stop("length(vb) not consistent with length(y)")
        temp <- deparse(substitute(vb))        
        nvb <- 1
        vb <- matrix(vb, ncol=1)
        colnames(vb) <- temp
      }
      else{
        if (nrow(vb) != ny) stop("nrow(vb) not consistent with length(y)")        
        nvb <- ncol(vb)
        vb <- as.matrix(vb)
      }  
    }  
  }  ## end of else (missing(vb) & !intcpt.random)

  ## Covariates for random effects, proportional w.r.t. odds
  if (missing(xb)){
    xb <- numeric(0)
    nxb <- 0
  }
  else{
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

  ## Fixed effects covariates proportional w.r.t. odds
  if (missing(x)){
    x <- numeric(0)
    nx <- 0
  }
  else{
    if (is.null(dim(x))){
      if (length(x) != ny) stop("length(x) not consistent with length(y)")
      temp <- deparse(substitute(x))
      nx <- 1
      x <- matrix(x, ncol=1)
      colnames(x) <- temp      
    }
    else{
      nx <- ncol(x)
      x <- as.matrix(x)
    }  
  }

  ## Fixed effects covariates not proportional w.r.t. odds
  if (missing(v) & !intercept){
    v <- numeric(0)
    nv <- 0
  }
  else{
    if (intercept){
      if (missing(v)){
        nv <- 1
        v <- matrix(rep(1, ny), ncol=1)
        colnames(v) <- "(Intercept)"
      }
      else{
        if(is.null(dim(v))){
          if (length(v) != ny) stop("length(v) not consistent with length(y)")
          temp <- deparse(substitute(v))
          nv <- 2
          v <- as.matrix(cbind(rep(1, ny), v))
          colnames(v) <- c("(Intercept)", temp)
        }  
        else{
          if (nrow(v) != ny) stop("nrow(v) not consistent with length(y)")          
          nv <- 1 + ncol(v)
          temp <- colnames(v)
          v <- as.matrix(cbind(rep(1, ny), v))
          colnames(v) <- c("(Intercept)", temp)
        }  
      }  
    }
    else{    ## else (intercept)
      if(is.null(dim(v))){
        if (length(v) != ny) stop("length(v) not consistent with length(y)")
        temp <- deparse(substitute(v))
        nv <- 1
        v <- matrix(v, ncol=1)
        colnames(v) <- temp
      }
      else{
        if (nrow(v) != ny) stop("nrow(v) not consistent with length(y)")        
        nv <- ncol(v)
        v <- as.matrix(v)        
      }  
    }  
  }  ## end of else (missing(v) & !intercept)

  ## Cluster indicator
  if (nxb + nvb){
    if (missing(cluster)) stop("cluster indicator must be given")
    if (length(cluster) != ny) stop("length(cluster) not consistent with length(y)")
  }
  else{
    cluster <- 1:ny
  }  

  ## Put everything together
  RESP <- data.frame(y=y, cluster=cluster)
  
  if (nv){
    RESP <- cbind(RESP, v)
    bindv <- rep(1:nv, C) + rep(0:(C-1), each=nv)*rep(nv+nvb, nv*C)  ## indeces of V in the regression coefficient vector    
    indv <- 2 + (1:nv)                                               ## indeces of V in RESP
  }
  else{
    bindv <- numeric(0)    
    indv <- numeric(0)    
  }  
  if (nvb){
    RESP <- cbind(RESP, vb)
    bindvb <- nv + (rep(1:nvb, C) + rep(0:(C-1), each=nvb)*rep(nv+nvb, nvb*C))  ## indeces of V(b) in the reg. coef. vector    
    indvb <- 2 + ((nv+1):(nv+nvb))                                              ## indeces of V(b) in RESP

  }
  else{
    bindvb <- numeric(0)    
    indvb <- numeric(0)
  }    
  
  if (nx){
    RESP <- cbind(RESP, x)
    bindx <- C*(nv + nvb) + (1:nx)               ## indeces of X in the regression coefficient vector    
    indx <- 2 + nv + nvb + (1:nx)                ## indeces of X in RESP

  }
  else{
    bindx <- numeric(0)    
    indx <- numeric(0)
  }  
  if (nxb){
    RESP <- cbind(RESP, xb)
    bindxb <- C*(nv + nvb) + nx + (1:nxb)       ## indeces of X(b) in the regression coefficient vector    
    indxb <- 2 + nv + nvb + nx + (1:nxb)        ## indeces of X(b) in RESP

  }
  else{
    bindxb <- numeric(0)    
    indxb <- numeric(0)
  }  

  ## Remove missing observations
  not.NA <- !apply(is.na(RESP), 1, sum)
  RESP <- RESP[not.NA,]

  ## Sort it with respect to the cluster
  RESP <- RESP[order(RESP$cluster),]

  ## Extract various components
  y <- as.numeric(RESP[, 1])
  cluster <- as.numeric(RESP[, 2])
  ni <- as.numeric(table(cluster))
  N <- length(ni)

  if (nv){
    v <- as.matrix(RESP[,indv])
    colnames(v) <- colnames(RESP)[indv]
  } 
  if (nvb){
    vb <- as.matrix(RESP[,indvb])
    colnames(vb) <- colnames(RESP)[indvb]
  }  
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
  if (!(nx+nv+nxb+nvb)) stop("Neither covariates nor intercept included in the model")

  yTAB <- table(y)
  yCAT <- as.numeric(names(yTAB))
  if (any(!(yCAT %in% 0:C))) stop("Response not from {0,...,C} found")

  xv.qr <- qr(RESP[,-(1:2)])
  if (xv.qr$rank < nv+nvb+nx+nxb & !predict){
    warning(paste("Rank of the matrix (V,V(b),X,X(b)) is ", xv.qr$rank, " < ", nv+nvb+nx+nxb, sep=""))
  }  

  ind.logit1 <- 0
  nvb4C <- nvb  
  nxb4C <- nxb
  bindvb4C <- 0   ## to be passed to C++ 
  bindxb4C <- 0   ## it may not be numeric(0) like bindvb or bindxb!
  if (!hierar.center){
    ## Design matrices for fixed effects -> put together V, V(b) and X, X(b)
    if (nv & nvb){
      v <- cbind(v, vb)
      nv <- nv + nvb
    }  
    else{
      if (nvb){
        v <- vb
        nv <- nvb
      }  
    }    
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
    
    ## Design matrices for random effects -> put together V(b), X(b)
    ## -> there are no "random effects" covariates not proportional w.r.t. odds
    if (nxb & nvb){
      xb <- cbind(vb, xb)
      vb <- numeric(0)
      ind.logit1 <- c(bindvb[1:nvb], bindxb)
      bindvb4C <- bindvb
      bindxb4C <- bindxb
      nxb <- nxb + nvb
      nvb <- 0
    }
    else{
      if (nvb){
        xb <- vb
        vb <- numeric(0)
        ind.logit1 <- bindvb[1:nvb]
        bindvb4C <- bindvb        
        nxb <- nvb
        nvb <- 0
      }
      else{
        ind.logit1 <- bindxb
        bindxb4C <- bindxb        
      }  
    }
  }  

  #cat("nvb4C=", nvb4C, "   bindvb = ", bindvb4C, "\n", sep="  ")
  #cat("nxb4C=", nxb4C, "   bindxb = ", bindxb4C, "\n", sep="  ")  
  
  RET <- list(y=y, cluster=cluster, C=C, x=x, v=v, xb=xb, vb=vb,
              ny=ny, N=N, ni=ni, nx=nx, nv=nv, nxb=nxb, nvb=nvb,
              bindx=bindx, bindv=bindv, bindxb=bindxb, bindvb=bindvb,
              nxb4C=nxb4C, nvb4C=nvb4C,
              bindxb4C=bindxb4C, bindvb4C=bindvb4C,
              ind.logit1=ind.logit1)
  return(RET)  
}  

### To be used inside 'cumlogit'
### -> no input checks
### ------------------------------------
###
### y, x, v should already be adjusted and must be matrices
###  (missing values removed etc.)
###
fit.cumlogit <- function(y, v, x, C=1, logit.order=c("decreasing", "increasing"), epsilon=1e-08, maxit=25, trace=FALSE)
{
  thispackage <- "glmmAK"
  #thispackage <- NULL

  n <- length(y)
  if (length(v)){
    nv <- ifelse(is.null(dim(v)), 1, ncol(v))
  }
  else{
    nv <- 0
  }    
  if (length(x)){
    nx <- ifelse(is.null(dim(x)), 1, ncol(x))
  }
  else{
    nx <- 0
  }  
  
  ## Order of logits
  logit.order <- match.arg(logit.order)  
  if (logit.order == "increasing"){
    yFIT <- C - y
  }
  else{
    yFIT <- y
  }  
  
  ## Initial values
  betaGamma <- init00.cumlogit(y=yFIT, v=v, x=x, C=C, only.intcpt=FALSE, logit.order=logit.order)

  ## Check whether some predicted category probabilities are not negative
  ## (might happen with betaGamma)
  probs <- prob00.cumlogit(coef=betaGamma, v=v, x=x, C=C, logit.order=logit.order)
  if (any(probs < 0)){
    betaGamma <- init00.cumlogit(y=yFIT, v=v, x=x, C=C, only.intcpt=TRUE, logit.order=logit.order)
    probs <- prob00.cumlogit(coef=betaGamma, v=v, x=x, C=C, logit.order=logit.order)
    if (any(probs < 0)) stop("Cannot determine initial values of regression coefficients leading to non-negative probabilities")
  }  

  ### Fit the full model
  nregr <- C*nv + nx
  FIT <- .C("fit_cumlogit",
               betaGamma=as.double(betaGamma),   ll=double(1),             U=double(nregr),
               I.obs=double(nregr^2),            I.exp=double(nregr^2),
               Y=as.integer(yFIT),               X=as.double(t(x)),        V=as.double(t(v)),
               C=as.integer(C),                  p=as.integer(nx),         q=as.integer(nv),        n=as.integer(n),
               niter=as.integer(maxit),          toler=as.double(epsilon), trace=as.integer(trace),
               err=integer(1), PACKAGE=thispackage
            )
  
  FIT$I.obs <- matrix(FIT$I.obs, nrow=nregr, ncol=nregr)
  FIT$I.exp <- matrix(FIT$I.exp, nrow=nregr, ncol=nregr)
  if (FIT$err){
    FIT$I.exp <- matrix(NA, nrow=nregr, ncol=nregr)
    if (FIT$err == 1) warning("Convergence not reached. Maximal number of step-halfing steps performed.")
    else
      if (FIT$err == 2) warning("Convergence not reached. Maximal number of iterations performed.")
      else
        if (FIT$err == 3) warning("Convergence not reached. Minus Hessian not positive definite.")
        else{
          FIT$I.obs <- matrix(NA, nrow=nregr, ncol=nregr)
          if (FIT$err == 100) warning("Convergence criterion satisfied but minus Hessian is not positive definite.")
          else
            if (FIT$err == 101) warning("Convergence not reached. Maximal number of step-halfing steps performed. Additionally: minus Hessian is not positive definite.")
            else            
              if (FIT$err == 102) warning("Convergence not reached. Maximal number of iterations performed. Additionally: minus Hessian is not positive definite.")
              else
                if (FIT$err == 106) warning("Convergence reached but the expected information matrix contains +/-Inf.")
                else
                  if (FIT$err == 107) warning("Convergence reached but the expected information matrix is not positive definite.")        
                  else{
                    cat("Error ", FIT$err, "\n", sep="")
                    stop("Numerical or other problems. No results produced.")
                  }  
        }           
  }   
  names(FIT$betaGamma) <- names(FIT$U) <- names(betaGamma)
  colnames(FIT$I.obs) <- rownames(FIT$I.obs) <- colnames(FIT$I.exp) <- rownames(FIT$I.exp)<- names(betaGamma)

  
  ### NULL MODEL (intercepts only)
  v.NULL <- matrix(rep(1, n), ncol=1)
  colnames(v.NULL) <- "(Intercept)"  
  x.NULL <- NULL
  betaGamma.NULL <- init00.cumlogit(y=y, v=v, x=x, C=C, only.intcpt=TRUE, logit.order=logit.order)
  betaGamma.NULL <- betaGamma.NULL[seq(1, (C-1)*nv+1, by=nv)]
  nregr.NULL <- C
  FIT.NULL <- .C("fit_cumlogit",
                 betaGamma=as.double(betaGamma.NULL),   ll=double(1),              U=double(nregr.NULL),
                 I.obs=double(nregr.NULL^2),            I.exp=double(nregr.NULL^2),                 
                 Y=as.integer(yFIT),                    X=double(0),               V=as.double(t(v.NULL)),
                 C=as.integer(C),                       p=as.integer(0),           q=as.integer(1),            n=as.integer(n),
                 niter=as.integer(maxit),               toler=as.double(epsilon),  trace=as.integer(trace),
                 err=integer(1), PACKAGE=thispackage
            )

  LogLik <- c(FIT$ll, FIT.NULL$ll)
  names(LogLik) <- c("model", "null")
  
  RET <- list(coefficients=FIT$betaGamma,
              loglik=LogLik,
              score=FIT$U,
              vcov=FIT$I.obs,
              expect.vcov=FIT$I.exp,
              logit.order=logit.order,
              linear.predictors=linear.predictors00.cumlogit(coef=FIT$betaGamma, v=v, x=x, C=C, logit.order=logit.order),
              fitted.values=prob00.cumlogit(FIT$betaGamma, v=v, x=x, C=C, logit.order=logit.order),
              converged=(FIT$err==0),
              iter=FIT$niter,
              C=C,
              y=y,
              x=x,
              v=v)

  class(RET) <- "cumlogit"
  return(RET)  
}  


### To be used inside 'cumlogit'
### -> no input checks
### ------------------------------------
##
## \item{only.intcpt}{logical, if TRUE only intercepts are estimated and all other regression coefficients are set to zero}
## \item{logit.order}{here only used to give correct labels to the coefficients. All initials are computed for
##     decreasing logits}
##
init00.cumlogit <- function(y, v, x, C, only.intcpt=FALSE, logit.order=c("decreasing", "increasing"))
{
  logit.order <- match.arg(logit.order)

  nv <- if (!length(v)) 0 else ncol(v)  
  nx <- if (!length(x)) 0 else ncol(x)
  
  if (nx & nv) XX <- cbind(v, x)
  if (nx & !nv) XX <- x
  if (!nx & nv) XX <- v
  
  COEF <- numeric()
  for (i in 1:C){
    YY <- (y >= i)
    if (only.intcpt) COEF <- rbind(COEF, glm(YY~1, family=binomial(link=logit))$coefficients)
    else             COEF <- rbind(COEF, glm.fit(x=XX, y=YY, family=binomial(link=logit), intercept=FALSE)$coefficients)
  }

  if (nv){
    if (only.intcpt){
      GAMMA <- cbind(COEF, matrix(0, nrow=C, ncol=nv-1))
      GAMMA <- as.numeric(t(GAMMA))
    }
    else{
      GAMMA <- COEF[, 1:nv]
      GAMMA <- as.numeric(t(GAMMA))
    }
    if (logit.order == "decreasing"){
      if (C == 1) names(GAMMA) <- colnames(v)[1:nv]
      else        names(GAMMA) <- paste(colnames(v)[1:nv], ":", rep(1:C, each=nv), sep="")
    }  
    else{
      if (C == 1) names(GAMMA) <- colnames(v)[1:nv]
      else        names(GAMMA) <- paste(colnames(v)[1:nv], ":", rep(C:1, each=nv), sep="")
    }  
  }
  else{
    GAMMA <- numeric(0)
  }  
  
  if (nx){
    if (only.intcpt){
      BETA <- rep(0, nx)
    }
    else{
      BETA <- matrix(COEF[, nv+(1:nx)], ncol=nx)
      BETA <- apply(BETA, 2, mean)
    }
    names(BETA) <- colnames(x)    
  }  
  else{
    BETA <- numeric(0)
  }
  
  COEF <- c(GAMMA, BETA)
  
  return(COEF)
}  


### To be used inside 'cumlogit'
### -> no input checks
### ------------------------------------
linear.predictors00.cumlogit <- function(coef, v, x, C, logit.order=c("decreasing", "increasing"))
{
  logit.order <- match.arg(logit.order)

  nv <- if (!length(v)) 0 else ncol(v)
  nnv <- if (!length(v)) 0 else nrow(v)
  
  nx <- if (!length(x)) 0 else ncol(x)
  nnx <- if (!length(x)) 0 else nrow(x)
  
  nn <- max(nnx, nnv)
  if (!nn) stop("No covariates supplied")

  etaV <- matrix(0, nrow=nn, ncol=C)
  if (nv){
    gamma <- coef[1:(C*nv)]
    for (i in 1:C){
      etaV[,i] <- v %*% gamma[((i-1)*nv+1):(i*nv)]
    }  
  }  
  
  if (nx){
    beta <- coef[C*nv + (1:nx)]
    etaX <- x %*% beta    
  }
  else{
    etaX <- rep(0, nn)
  }
  
  eta <- matrix(rep(etaX, C), ncol=C) + etaV
  if (logit.order == "decreasing"){
    if (C == 1) colnames(eta) <- "eta"
    else        colnames(eta) <- paste("eta", 1:C, sep="")
  }  
  else{
    if (C == 1) colnames(eta) <- "eta"
    else        colnames(eta) <- paste("eta", C:1, sep="")
  }  

  return(eta)  
}  

### To be used inside 'cumlogit'
### -> no input checks
### ------------------------------------
prob00.cumlogit <- function(coef, v, x, C, logit.order=c("decreasing", "increasing"))
{
  logit.order <- match.arg(logit.order)  
  eta <- linear.predictors00.cumlogit(coef=coef, v=v, x=x, C=C, logit.order=logit.order)
  
  hh <- exp(eta)/(1 + exp(eta))
  pi <- 1 - hh[,1]
  if (C > 1){
    for (i in 1:(C-1)){
      pi <- cbind(pi, hh[,i] - hh[,i+1])
    }
  }
  pi <- cbind(pi, hh[,C])
  if (logit.order == "decreasing"){
    colnames(pi) <- paste("prob", 0:C, sep="")
  }
  else{
    colnames(pi) <- paste("prob", C:0, sep="")
  }  
  
  return(pi)  
}  
