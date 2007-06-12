#*** logpoisson.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              akom@email.cz
##
##           CREATED:  08/02/2007
##
## PURPOSE: Poisson regression (log-linear model)
##
## FUNCTIONS: logpoisson
##            print.logpoisson
##            summary.logpoisson
##
#* ********************************************************************************* */

logpoisson <- function(y, x, offset=0, epsilon=1e-08, maxit=25, trace=FALSE)
{  
  ## Design matrices and check of some input arguments
  DES <- design.logpoisson(y=y, x=x, offset=offset, intcpt.random=FALSE, hierar.center=FALSE)
  
  ## xnames
  if (missing(x)) colnames(DES$x) <- "(Intercept)"
  else{
    if (is.null(dim(x))) xname <- deparse(substitute(x))
    else                 xname <- colnames(x)
    colnames(DES$x) <- c("(Intercept)", xname)
  }  

  RET <- fit.logpoisson(y=DES$y, x=DES$x, offset=DES$offset, epsilon=epsilon, maxit=maxit, trace=trace)
  return(RET)
}


print.logpoisson <- function(x, ...)
{
  VCOV <- x$vcov

  ncoef <- length(x$coefficients)
  if (ncoef == 1) sdCoef <- sqrt(VCOV)
  else{
    sdCoef <- sqrt(VCOV[cbind(1:ncoef, 1:ncoef)])
  }  

  Zvalue <- x$coefficients/sdCoef
  pvalue <- 2 * pnorm(-abs(Zvalue))
  coef.table <- cbind(x$coefficients, sdCoef, Zvalue, pvalue)
  dimnames(coef.table) <- list(names(x$coefficients), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  
  LR.test <- 2*(x$loglik["model"] - x$loglik["null"])
  df.LR.test <- length(x$coefficients) - 1
  if (df.LR.test) p.LR.test <- pchisq(LR.test, df.LR.test, lower.tail=FALSE)
  else            p.LR.test <- NA

  cat("\n      Poisson log-linear model")
  cat("\n\n")
  print(coef.table)
  cat("\n")
  cat("\n             Log-likelihood: ", x$loglik["model"], sep="")
  cat("\nLog-likelihood (null model): ", x$loglik["null"], sep="")
  cat("\n\n    -2 Log-Likelihood ratio: ", 2*(x$loglik["model"] - x$loglik["null"]), sep="")
  cat("\n    on ", df.LR.test, " degrees of freedom,  P-value = ", p.LR.test, "\n", sep="")
  
##  cat("\nScore:\n")
##  print(x$score)

##  cat("\nvcov:\n")
##  print(x$vcov)
  
  if (x$converged) cat("\nConverged after ", x$iter, " iterations.\n", sep="")
  else             cat("\nNot converged (", x$iter, " iterations performed).\n", sep="")
  
  return(invisible(x))      
}


summary.logpoisson <- function(object, ...)
{
  print(object, ...)
  return(invisible(object))
}  
