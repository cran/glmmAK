#*** copula.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              akom@email.cz
##
##         CREATED:  26/10/2006
##
## PURPOSE: Several copulas
##
## FUNCTIONS: 
## 
#* ********************************************************************************* */

### =========================================
### Copula's joint distribution functions
### =========================================
Cplackett <- function(u, v, theta=1)
{
  if (any(u < 0) | any(u > 1)) stop("All components of u must lie between 0 and 1")
  if (any(v < 0) | any(v > 1)) stop("All components of v must lie between 0 and 1")  
  if (theta <= 0) stop("theta must be positive")
  if (length(u) != length(v)) stop("u and v must have the same length")
  
  if (theta == 1){
    RES <- u*v
  }
  else{
    th1 <- theta - 1
    sumx <- u + v
    P <- 1 + th1*sumx
    sD <- sqrt(P^2 - 4*u*v*theta*th1)
    RES <- (P - sD)/(2*th1)
  }  

  return(RES)
}  

Cgauss <- function(u, v, theta=0)
{
  require(mvtnorm)

  if (any(u < 0) | any(u > 1)) stop("All components of u must lie between 0 and 1")
  if (any(v < 0) | any(v > 1)) stop("All components of v must lie between 0 and 1")
  if (theta < -1 | theta > 1) stop("theta must be >= -1 and <= 1")  
  if (length(u) != length(v)) stop("u and v must have the same length")

  iPhiu <- qnorm(u)
  iPhiv <- qnorm(v)
  Sigma <- matrix(c(1, theta, theta, 1), nrow=2)
  RES <- numeric(length(u))
  for (i in 1:length(u)) RES[i] <- pmvnorm(lower=rep(-Inf, 2), upper=c(iPhiu[i], iPhiv[i]), mean=rep(0, 2), sigma=Sigma)
  if (is.matrix(u)){
    RES <- matrix(RES, ncol=ncol(u), nrow=nrow(u))
  }  
  return(RES)
}  

Cclayton <- function(u, v, theta=0)
{
  if (any(u < 0) | any(u > 1)) stop("All components of u must lie between 0 and 1")
  if (any(v < 0) | any(v > 1)) stop("All components of v must lie between 0 and 1")
  if (theta <= -1) stop("theta must be > -1")  
  if (length(u) != length(v)) stop("u and v must have the same length")
  
  if (theta == 0){
    RES <- u*v
  }
  else{
    RES <- (u^(-theta) + v^(-theta) - 1)^(-1/theta)
    RES[RES <= 0] <- 0
  }   

  return(RES)
}


### =========================================
### Copula's joint densities
### =========================================
cplackett <- function(u, v, theta=1)
{
  if (any(u < 0) | any(u > 1)) stop("All components of u must lie between 0 and 1")
  if (any(v < 0) | any(v > 1)) stop("All components of v must lie between 0 and 1")  
  if (theta <= 0) stop("theta must be positive")
  if (length(u) != length(v)) stop("u and v must have the same length")
  
  if (theta == 1){
    RES <- rep(1, length(u))
    if (is.matrix(u)){
      RES <- matrix(RES, ncol=ncol(u), nrow=nrow(u))
    }      
  }
  else{
    th1 <- theta - 1
    sumx <- u + v
    P <- 1 + th1*sumx
    sD <- sqrt(P^2 - 4*u*v*theta*th1)
    RES <- (th1*(P - 2*u*theta)*(P - 2*v*theta))/(2*sD^3) + (theta+1)/(2*sD)
  }  

  return(RES)
}  

cgauss <- function(u, v, theta=0)
{
  require(mvtnorm)

  if (any(u < 0) | any(u > 1)) stop("All components of u must lie between 0 and 1")
  if (any(v < 0) | any(v > 1)) stop("All components of v must lie between 0 and 1")
  if (theta < -1 | theta > 1) stop("theta must be >= -1 and <= 1")  
  if (length(u) != length(v)) stop("u and v must have the same length")

  iPhiu <- qnorm(u)
  iPhiv <- qnorm(v)
  DENOM1 <- dnorm(iPhiu)
  DENOM2 <- dnorm(iPhiv)
  DENOM <- DENOM1 * DENOM2
  Sigma <- matrix(c(1, theta, theta, 1), nrow=2)
  RES <- numeric(length(u))
  for (i in 1:length(u)){
    NUMER <- dmvnorm(x=c(iPhiu[i], iPhiv[i]), mean=rep(0, 2), sigma=Sigma)
    RES[i] <- NUMER/DENOM[i]
  }  
  if (is.matrix(u)){
    RES <- matrix(RES, ncol=ncol(u), nrow=nrow(u))
  }  
  return(RES)
}  

cclayton <- function(u, v, theta=0)
{
  if (any(u < 0) | any(u > 1)) stop("All components of u must lie between 0 and 1")
  if (any(v < 0) | any(v > 1)) stop("All components of v must lie between 0 and 1")
  if (theta <= -1) stop("theta must be > -1")  
  if (length(u) != length(v)) stop("u and v must have the same length")
  
  if (theta == 0){
    RES <- rep(1, length(u))
    if (is.matrix(u)){
      RES <- matrix(RES, ncol=ncol(u), nrow=nrow(u))
    }      
  }
  else{
    itheta <- 1/theta
    RES <- (itheta+1)*theta*u^(-(theta+1))*v^(-(theta+1))*(1/v^theta + 1/u^theta)^(-(itheta+2))
    RES[RES <= 0] <- 0
  }   

  return(RES)
}




