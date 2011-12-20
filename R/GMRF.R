#*** GMRF.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              arnost.komarek[AT]mff.cuni.cz
##
##           CREATED:  13/11/2006
##
## PURPOSE:  Gaussian Markov random fields
##           - random numbers, densities etc.
##
## FUNCTIONS:   momentsGMRF
##              rGMRF
##              dGMRF
##              dGMRF2
##
#* ********************************************************************


##########################################################################################
### momentsGMRF:
###
##########################################################################################
momentsGMRF <- function(mean=0, Q=1, Sigma, A, b=0)
{
  ## Q from Sigma
  if (!missing(Sigma)){
    if (is.null(dim(Sigma))){
      if (length(Sigma) > 1) stop("Sigma must be either a number of a square matrix")
      nx <- 1
    }  
    else{
      if (ncol(Sigma) != nrow(Sigma)) stop("Sigma must be a square matrix")
      nx <- ncol(Sigma)
    }  
    LSigma <- chol(Sigma)
    Q <- chol2inv(LSigma)
  }

  ## Decomposition of Q and Sigma
  if (is.null(dim(Q))){
    if (length(Q) > 1) stop("Q must be either a number of a square matrix")
    nx <- 1
    Q <- matrix(Q, nrow=1, ncol=1)
  }  
  else{
    if (ncol(Q) != nrow(Q)) stop("Sigma must be a square matrix")
    nx <- ncol(Q)
  }    
  Li <- chol(Q)
  Sigma <- chol2inv(Li)
  Li <- t(Li)                      ## I want to have lower triangle

  ## Unconstrained mean
  if (length(mean) == 1){
    mean <- rep(mean, nx)
    names(mean) <- paste("x", 1:nx, sep="")
  }
  else{
    if (length(mean) != nx) stop(paste("mean must be of length ", nx, sep=""))
    if (is.null(names(mean))) names(mean) <- paste("x", 1:nx, sep="")
  }
  rownames(Sigma) <- colnames(Sigma) <- rownames(Li) <- colnames(Li) <- names(mean)

  ## Possible constraints
  if (missing(A)){
    nc <- 0
  }  
  else{
    if (is.null(dim(A))){
      if (length(A) != nx) stop(paste("A must be a vector of length ", nx, " or a matrix with ", nx, " columns", sep=""))
      nc <- 1
      A <- matrix(A, nrow=1, ncol=nx)
    }
    else{
      if (ncol(A) != nx) stop(paste("A must be a vector of length ", nx, " or a matrix with ", nx, " columns", sep=""))
      nc <- nrow(A)      
    }
    colnames(A) <- names(mean)    

    if (length(b) == 1){
      b.nonZERO <- (b != 0)
      b <- rep(b, nc)    
    }
    else{
      if (length(b) != nc) stop(paste("b must be of length ", nc, sep=""))
    }
    rownames(A) <- names(b) <- paste("constraint", 1:nc, sep="")    
  }

  ## Compute moments of the constrained GMRF
  if (nc){
    Sigma.tA <- Sigma %*% t(A)
    AQiA <- A %*% Sigma.tA
    iAQiA <- chol2inv(chol(AQiA))
    Hat <- Sigma.tA %*% iAQiA

    mean.star <-  mean - Hat %*% (A %*% mean - b)
    Sigma.star <- Sigma - Hat %*% t(Sigma.tA)

    mean.star <- as.numeric(mean.star)
    names(mean.star) <- rownames(Sigma.star) <- colnames(Sigma.star) <- names(mean)
  }
  else{
    mean.star <-  mean
    Sigma.star <- Sigma    
  }  

  RET <- list(mean=mean.star,  Sigma=Sigma.star)
  attr(RET, "mean.unconstr") <- mean
  attr(RET, "Sigma.unconstr") <- Sigma
  attr(RET, "Q.unconstr.cholesky") <- Li
  attr(RET, "nconstraint") <- nc
  if (nc){
    attr(RET, "A") <- A
    attr(RET, "b") <- b
  }  

  return(RET)
}


##########################################################################################
### rGMRF: Random number generation from GMRF
###
##########################################################################################
rGMRF <- function(n, mean=0, Q=1, Sigma, A, b=0)
{
  thispackage <- "glmmAK"
  #thispackage <- NULL

  ## moments and related stuff
  MOM <- momentsGMRF(mean=mean, Q=Q, Sigma=Sigma, A=A, b=b)  
  nx <- length(MOM$mean)
  
  ## number of sampled points
  if (n <= 0) stop("n must be positive")

  ## parameters of GMRF
  mean <- attr(MOM, "mean.unconstr")
  Li <- attr(MOM, "Q.unconstr.cholesky")
  nc <- attr(MOM, "nconstraint")

  Li.LT <- Li[lower.tri(Li, diag=TRUE)]  
  mu.nonZERO <- any(mean != 0)  
  if (nc){
    A <- attr(MOM, "A")  
    b <- attr(MOM, "b")    
    b.nonZERO <- any(b != 0)
  }
  else{
    A <- 0
    b <- 0
    b.nonZERO <- 0
  }  
  
  ## Sample
  SAMPLE <- .C("rGMRFR", x=double(n*nx),
                         log.dens=double(n),
                         mu=as.double(mean),
                         Li=as.double(Li.LT),
                         A=as.double(A),
                         e=as.double(b),
                         nx=as.integer(nx),
                         nc=as.integer(nc),
                         nrandom=as.integer(n),
                         mu.nonZERO=as.integer(mu.nonZERO),
                         e.nonZERO=as.integer(b.nonZERO),
		PACKAGE=thispackage)

  if (n == 1){
    names(SAMPLE$x) <-  names(mean)
    names(SAMPLE$log.dens) <- 1:n
  }
  else{
    SAMPLE$x <- matrix(SAMPLE$x, byrow=TRUE, ncol=nx, nrow=n)
    colnames(SAMPLE$x) <- names(mean)
    names(SAMPLE$log.dens) <- rownames(SAMPLE$x) <- 1:n    
  }

  return(list(x=SAMPLE$x, log.dens=SAMPLE$log.dens))
}


##########################################################################################
### dGMRF: (Log-)density of GMRF, computed using the factorization of pi(x|Ax)
###
##########################################################################################
dGMRF <- function(x, mean=0, Q=1, Sigma, A, b=0, log=FALSE)
{
  thispackage <- "glmmAK"
  #thispackage <- NULL

  ## moments and related stuff
  MOM <- momentsGMRF(mean=mean, Q=Q, Sigma=Sigma, A=A, b=b)
  nx <- length(MOM$mean)
  
  ## number of points where to evaluate density
  if (is.null(dim(x))){
    if (nx == 1){
      npoints <- length(x)
    }
    else{
      if (length(x) != nx) stop(paste("Dimension of GMRF is ", nx, " and x has length ", length(x), sep=""))
      npoints <- 1
    }  
  }
  else{
    if (ncol(x) != nx) stop(paste("Dimension of GMRF is ", nx, " and x has ", ncol(x), " columns", sep=""))
    npoints <- nrow(x)
  }  

  ## parameters of GMRF  
  mean <- attr(MOM, "mean.unconstr")
  Li <- attr(MOM, "Q.unconstr.cholesky")
  nc <- attr(MOM, "nconstraint")

  Li.LT <- Li[lower.tri(Li, diag=TRUE)]  
  mu.nonZERO <- any(mean != 0)  
  if (nc){
    A <- attr(MOM, "A")  
    b <- attr(MOM, "b")    
    b.nonZERO <- any(b != 0)
  }
  else{
    A <- 0
    b <- 0   
    b.nonZERO <- 0
  }  
  
  ## Evaluate (log-)density
  DENSITY <- .C("dGMRFR", value=double(npoints),
                          x=as.double(t(x)),
                          unlog=as.integer(!log),
                          mu=as.double(mean),
                          Li=as.double(Li.LT),
                          A=as.double(A),
                          e=as.double(b),
                          nx=as.integer(nx),
                          nc=as.integer(nc),
                          npoints=as.integer(npoints),
                          mu.nonZERO=as.integer(mu.nonZERO),
                          e.nonZERO=as.integer(b.nonZERO), 
		PACKAGE=thispackage)

  return(as.numeric(DENSITY$value))
}


##########################################################################################
### dGMRF2: (Log-)density of GMRF, computed using the eigen-value decomposition of
###         Sigma(star)
###
##########################################################################################
dGMRF2 <- function(x, mean=0, Q=1, Sigma, A, b=0, log=FALSE)
{
  thispackage <- "glmmAK"
  #thispackage <- NULL

  ## moments and related stuff
  MOM <- momentsGMRF(mean=mean, Q=Q, Sigma=Sigma, A=A, b=b)
  nx <- length(MOM$mean)
  
  ## number of points where to evaluate density
  if (is.null(dim(x))){
    if (nx == 1){
      npoints <- length(x)
    }
    else{
      if (length(x) != nx) stop(paste("Dimension of GMRF is ", nx, " and x has length ", length(x), sep=""))
      npoints <- 1
    }  
  }
  else{
    if (ncol(x) != nx) stop(paste("Dimension of GMRF is ", nx, " and x has ", ncol(x), " columns", sep=""))
    npoints <- nrow(x)
  }  

  ## parameters of GMRF  
  mean <- attr(MOM, "mean.unconstr")
  Li <- attr(MOM, "Q.unconstr.cholesky")
  nc <- attr(MOM, "nconstraint")

  Li.LT <- Li[lower.tri(Li, diag=TRUE)]  
  mu.nonZERO <- any(mean != 0)  
  if (nc){
    A <- attr(MOM, "A")  
    b <- attr(MOM, "b")    
    b.nonZERO <- any(b != 0)
  }
  else{
    A <- 0
    b <- 0   
    b.nonZERO <- 0
  }  
  
  ## Evaluate (log-)density
  DENSITY <- .C("dGMRF2R", value=double(npoints),
                           mustar=double(nx),
                           LiSigmastar=double(nx*nx),
                           x=as.double(t(x)),
                           unlog=as.integer(!log),
                           mu=as.double(mean),
                           Li=as.double(Li.LT),
                           A=as.double(A),
                           e=as.double(b),
                           nx=as.integer(nx),
                           nc=as.integer(nc),
                           npoints=as.integer(npoints),
                           mu.nonZERO=as.integer(mu.nonZERO),
                           e.nonZERO=as.integer(b.nonZERO),
		PACKAGE=thispackage)

  return(as.numeric(DENSITY$value))
}


