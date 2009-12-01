#*** cumlogit.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              akom@email.cz
##
##           CREATED:   07/08/2006
##   WORKING VERSION:   16/08/2006
## MAJOR MODIFICATION:  11/10/2006
##                      (internally, regression coefficients are stored in the order 'v' covariates first, 'x' covariates then)
##                      01/12/2009  bug in the print method fixed
##
## PURPOSE: Cumulative logit model
##
## FUNCTIONS: cumlogit
##            print.cumlogit
##            summary.cumlogit
## 
#* ********************************************************************************* */

cumlogit <- function(y, v, x, C=1, logit.order=c("decreasing", "increasing"), epsilon=1e-08, maxit=25, trace=FALSE)
{
  ## Design matrices and check of some input arguments
  DES <- design.cumlogit(y=y, v=v, x=x, intcpt.random=FALSE, hierar.center=FALSE, C=C)
  
  ## Fit the model
  RET <- fit.cumlogit(y=DES$y, v=DES$v, x=DES$x, C=C, logit.order=logit.order, epsilon=epsilon, maxit=maxit, trace=trace)
  return(RET)  
}  


print.cumlogit <- function(x, vcov=c("observed", "expected"), ...)
{
  type.vcov <- match.arg(vcov)  
  if (type.vcov == "observed") VCOV <- x$vcov
  else
    if (type.vcov == "expected") VCOV <- x$expect.vcov
  
  ncoef <- length(x$coefficients)
  if (ncoef == 1) sdCoef <- sqrt(VCOV)
  else{
    sdCoef <- sqrt(VCOV[cbind(1:ncoef, 1:ncoef)])
  }  
  
  Zvalue <- x$coefficients/sdCoef
  pvalue <- 2 * pnorm(-abs(Zvalue))
  coef.table <- cbind(x$coefficients, sdCoef, Zvalue, pvalue)
  dimnames(coef.table) <- list(names(x$coefficients), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

  nv <- ifelse(length(x$v), ncol(x$v), 0)
  nx <- ifelse(length(x$x), ncol(x$x), 0)
  
  LR.test <- 2*(x$loglik["model"] - x$loglik["null"])
  df.LR.test <- length(x$coefficients) - x$C
  if (df.LR.test) p.LR.test <- pchisq(LR.test, df.LR.test, lower.tail=FALSE)
  else            p.LR.test <- NA

  if (x$C > 1){
    cat("\n      Cumulative logit model")
    cat("\n      Response categories: ")
    cat(0:x$C, sep=", ")
    cat("\n      Logit order: ")
    cat(x$logit.order)
    cat("\n\n")
  }
  else{
    cat("\n      Logit model for binary response\n\n")
  }

  if (nv){
    cat("\nCOVARIATES WHOSE EFFECT IS NOT PROPORTIONAL w.r.t. ODDS\n")
#    cat("-------------------------------------------------------\n")
    cat("\n")
    for (s in 1:x$C){
      if (x$logit.order == "decreasing"){      
        cat("     Logit: log[P(Y >= ", s, ")/P(Y <= ", s-1, ")]\n", sep="")
        indNow <- if (nv > 1) ((s-1)*nv+1):(s*nv) else s
      }  
      else{
        cat("     Logit: log[P(Y <= ", s-1, ")/P(Y >= ", s, ")]\n", sep="")
        indNow <- if (nv > 1) ((x$C-s)*nv+1):((x$C-s+1)*nv) else (x$C-s+1)
      }  
      if (nv > 1){
        PRINT <- coef.table[indNow,]
      }  
      else{
        PRINT <- matrix(coef.table[indNow,], nrow=1)
        colnames(PRINT) <- colnames(coef.table)
        rownames(PRINT) <- rownames(coef.table)[indNow]
      }  
      print(PRINT)
      cat("\n")
    }    
  }

  if (nx){
    cat("\nCOVARIATES WHOSE EFFECT IS PROPORTIONAL w.r.t. ODDS\n")
#    cat("---------------------------------------------------\n")
    cat("\n")
    if (x$logit.order == "decreasing")
      cat("     Logit: log[P(Y >= ", "s", ")/P(Y <= ", "s-1", ")]\n", sep="")
    else
      cat("     Logit: log[P(Y <= ", "s-1", ")/P(Y >= ", "s", ")]\n", sep="")      
    if (nx > 1){
      PRINT <- coef.table[x$C*nv + (1:nx),]      
    }
    else{
      PRINT <- matrix(coef.table[x$C*nv + 1,], nrow=1)
      colnames(PRINT) <- colnames(coef.table)
      rownames(PRINT) <- rownames(coef.table)[1]      
    }  
    print(PRINT)
  }    

  cat("\n")
  cat("\n             Log-likelihood: ", x$loglik["model"], sep="")
  cat("\nLog-likelihood (null model): ", x$loglik["null"], sep="")
  cat("\n\n    -2 Log-Likelihood ratio: ", 2*(x$loglik["model"] - x$loglik["null"]), sep="")
  cat("\n    on ", df.LR.test, " degrees of freedom,  P-value = ", p.LR.test, "\n", sep="")

##  cat("\nScore:\n")
##  print(x$score)

##  cat("\nvcov (observed):\n")
##  print(x$vcov)

##  cat("\nvcov (expected):\n")
##  print(x$expect.vcov)
  
  if (x$converged) cat("\nConverged after ", x$iter, " iterations.\n", sep="")
  else             cat("\nNot converged (", x$iter, " iterations performed).\n", sep="")
  
  return(invisible(x))  
}  

summary.cumlogit <- function(object, vcov=c("observed", "expected"), ...)
{
  print(object, vcov=vcov, ...)
  return(invisible(object))
}  
