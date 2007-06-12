/*** Mvtdist2.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  13/09/2006
//
//                 dmvtnorm2:  14/09/2006
//    dmvtnorm2a, dmvtnorm2b:  02/10/2006
//
// PURPOSE: Utilities for several multivariate distributions
//          
//     * partially taken from mvtdist.cpp and mvtdist.h of the bayesSurv package
//
//
/* ********************************************************************************* */

#include "Mvtdist2.h"

namespace Mvtdist2 {

/* ********************************************************************************* */
/* rmvtnorm2: Sample from N(m, V^{-1})                                               */
/*               where m = V^{-1}*b                                                  */
/*                                                                                   */
/* INPUT: V ..... inverse variance of the normal distribution                        */
/*        b ..... vector b (as above)                                                */
/*                  its length should be equal to V->nrow()                          */
/*                                                                                   */
/* OUTPUT: z .... vector sampled from N(m, V^{-1})                                   */
/*         V .... Cholesky decomposition of V(IN), that is L where V(IN) = L*t(L)    */
/*         b ....                                                                    */
/*                                                                                   */
/*           * in the case that rank(V) < nrow(V), nothing is sampled                */
/*             and neither V nor b are modified                                      */
/*                                                                                   */
/* RETURN: rank of V                                                                 */
/*                                                                                   */
/* ALGORITHM: V = L*t(L)                                                             */
/*            V^{-1} = t(L^{-1})*L^{-1}                                              */
/*            U ~ N(0, eye)                                                          */
/*                                                                                   */
/*            Z = m + t(L^{-1})*U                 ~ N(m, V^{-1})                     */
/*              = t(L^{-1}) * (L^{-1}*b + U)                                         */
/*                                                                                   */
/*    1) solve L*x = b                                                               */
/*               x = L^{-1}*b                                                        */
/*                                                                                   */
/*    2) solve t(L)*Z = x + U                                                        */
/*                  Z = t(L^{-1}) * (L^{-1}*b + U)                                   */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
/** NOT IMPLEMNETED, but the algorithm is used on several places                    **/
/** ------------------------------------------------------------------------------- **/


/* ********************************************************************************* */
/* dmvtnorm2: Normal density or its logarithm                                        */
/*                                                                                   */
/* RETURN: computed value                                                            */
/* INPUT:                                                                            */
/*        x_m:  vector holding (x - mean), can be both row or column                 */
/*         Li:  Cholesky decomposition of the inverse variance matrix                */
/*              that is Sigma^{-1} = Li*t(Li)                                        */
/*               Li.nrow() determines the dimension of the normal distribution       */
/*       dlog:  1 when the logarithm of the density is to be computed                */
/*              0 otherwise                                                          */
/*              DEFAULT is 0                                                         */
/*   fullDens:  1 when the density is computed including (2*pi)^{-q/2} factor        */
/*              0 when this factor is ignored                                        */
/*              DEFAULT is 1                                                         */
/*                                                                                   */
/* OUTPUT:                                                                           */
/*        x_m:  x_m(OUT) = t(x_m(IN))*L                                              */
/*              (always either row or column vector in the same way as in INPUT)     */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
double
dmvtnorm2(MatrixRect<double> &x_m,  const MatrixLT<double> &Li,  const int &dlog,  const int &fullDens)
{
  static double VALUE, logDETER;
  static int i, j;
  double *resultP, *xP;
  const double *LiP;
  
  /*** (x - mean)'*Li             ***/
  /*** and log(det(Sigma)^{-1/2}) ***/
  resultP = x_m.a();
  LiP = Li.aconst();
  logDETER = 0;
  for (i = 0; i < Li.nrow(); i++){
    *resultP *= *LiP;

    if (*LiP < _AK_ZERO){
      VALUE = (dlog ? R_NegInf : 0);
      return VALUE;
    }
    logDETER += log(*LiP);

    LiP++;
    xP = resultP + 1;
    for (j = i+1; j < Li.nrow(); j++){
      *resultP += (*xP)*(*LiP);
      xP++;
      LiP++;
    }
    resultP++;
  }

  /*** -(1/2)*(x - mean)'Sigma^{-1}(x - mean) ***/
  VALUE = x_m.xTx();
  VALUE *= (-0.5);

  /*** log(density) ***/
  VALUE += logDETER;

  /*** final adjustements ***/
  if (fullDens) VALUE -= Li.nrow()*_AK_LOG_SQRT_2PI;
  if (!dlog) VALUE = exp_AK(VALUE);

  return(VALUE);
}


/* ********************************************************************************* */
/* dmvtnorm2a: Normal density or its logarithm.                                      */
/*             Neither (2*pi)^(-q/2), nor |Sigma|^(-1/2) factor are computed.        */
/*             That is -0.5*(x-mu)'Sigma^(-1)*(x-mu) or its exp() is computed        */
/*                                                                                   */
/* RETURN: computed value                                                            */
/* INPUT:                                                                            */
/*        x_m:  working vector of length Li.nrow()                                   */
/*          x:  vector holding x                                                     */
/*         mu:  vector holding mean                                                  */
/*         Li:  Cholesky decomposition of the inverse variance matrix                */
/*              that is Sigma^{-1} = Li*t(Li)                                        */
/*              Li.nrow() determines the dimension of the normal distribution        */
/*       dlog:  1 when the logarithm of the density is to be computed                */
/*              0 otherwise                                                          */
/*              DEFAULT is 0                                                         */
/*                                                                                   */
/* OUTPUT:                                                                           */
/*        x_m:  x_m(OUT) = t(x-mu)*L                                                 */
/*              (always either row or column vector in the same way as in INPUT)     */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
double
dmvtnorm2a(MatrixRect<double> &x_m,  const double *x,  const double *mu,  const MatrixLT<double> *Li,  const int &dlog)
{
  static double VALUE;
  static int i, j;
  double *resultP, *x_mP;
  const double *xP, *muP, *LiP;

  /*** x - mean ***/
  x_mP = x_m.a();
  xP = x;
  muP = mu;
  for (i = 0; i < x_m.length(); i++){
    *x_mP = (*xP) - (*muP);
    x_mP++;
    xP++;
    muP++;
  }

  /*** (x - mean)'*Li             ***/
  resultP = x_m.a();
  LiP = Li->aconst();
  for (i = 0; i < Li->nrow(); i++){
    *resultP *= *LiP;

    if (*LiP < _AK_ZERO){
      VALUE = (dlog ? R_NegInf : 0);
      return VALUE;
    }

    LiP++;
    x_mP = resultP + 1;
    for (j = i+1; j < Li->nrow(); j++){
      *resultP += (*x_mP)*(*LiP);
      x_mP++;
      LiP++;
    }
    resultP++;
  }

  /*** -(1/2)*(x - mean)'Sigma^{-1}(x - mean) ***/
  VALUE = x_m.xTx();
  VALUE *= (-0.5);

  /*** final adjustements ***/
  if (!dlog) VALUE = exp_AK(VALUE);

  return(VALUE);
}


/* ********************************************************************************* */
/* dmvtnorm2b: Normal density or its logarithm.                                      */
/*             Neither (2*pi)^(-q/2), nor |Sigma|^(-1/2) factor are computed.        */
/*             That is -0.5*(x-mu)'Sigma^(-1)*(x-mu) or its exp() is computed        */
/*                                                                                   */
/*    It is assumed that Sigma^(-1) is diagonal                                      */
/*                                                                                   */
/* RETURN: computed value                                                            */
/* INPUT:                                                                            */
/*          x:  vector holding x                                                     */
/*         mu:  vector holding mean                                                  */
/*     invVar:  diagonal elements of the inverse of the covariance matrix            */
/*        dim:  dimension of the normal distribution                                 */      
/*       dlog:  1 when the logarithm of the density is to be computed                */
/*              0 otherwise                                                          */
/*              DEFAULT is 0                                                         */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
double
dmvtnorm2b(const double *x,  const double *mu,  const double *invVar,  const int &dim,  const int &dlog)
{
  static double VALUE, x_mu;
  static int i;
  const double *xP, *muP, *iSigmaP;  
  
  xP = x;
  muP = mu;
  iSigmaP = invVar;
  x_mu = (*xP) - (*muP);
  VALUE = x_mu * (*iSigmaP) * x_mu;
  for (i = 1; i < dim; i++){
    xP++;
    muP++;
    iSigmaP++;
    x_mu = (*xP) - (*muP);   
    VALUE += x_mu * (*iSigmaP) * x_mu;
  }
  VALUE *= (-0.5);

  if (!dlog) VALUE = exp_AK(VALUE);

  return(VALUE);
}

}    /*** end of the namespace Mvtdist2   ***/

