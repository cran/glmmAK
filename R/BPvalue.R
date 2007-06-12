#*** BPvalue.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              akom@email.cz
##
##           CREATED:  31/05/2007
##
## PURPOSE: Univariate Bayesian P-value
##          (pseudo-contour probability)       
##
## FUNCTIONS: BPvalue
## 
#* ********************************************************************************* */

BPvalue <- function(sample)
{
  if (is.null(dim(sample))){
    negative <- sum(sample < 0)/length(sample)
    positive <- sum(sample > 0)/length(sample)
    p.value <- 2*min(negative, positive)
  }
  else{
    negative <- apply(sample < 0, 2, sum)/dim(sample)[1]
    positive <- apply(sample > 0, 2, sum)/dim(sample)[1]
    p.value <- 2*apply(rbind(negative, positive), 2, min)
  }

  return(p.value)
}  
