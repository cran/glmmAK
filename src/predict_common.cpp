/*** predict_common.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:          16/02/2007
//
//
//    PURPOSE:  Some common utilities for prediction in GALMMs
//
//           30/01/2007:   predict_Nrandom_effect_Zero
//           30/01/2007:   predict_Nrandom_effect
//           31/01/2007:   predict_G1random_effect_Zero
//           31/01/2007:   predict_G1random_effect
//           12/04/2007:   predict_GBirandom_effect_Zero
//           12/04/2007:   predict_GBirandom_effect
//           30/01/2007:   print_iter_info2
//
/* ********************************************************************************* */

#include "predict_common.h"

namespace predict_common{

/*** ============================================================================================= ***/
/*** predict_Nrandom_effects_Zero:  Predict the values of zero-mean normal random effects          ***/
/*** predict_Nrandom_effects:       Predict the values of general-mean normal random effects       ***/
/***                                                                                               ***/
/***                                                                                               ***/
/*** ============================================================================================= ***/
//
// bb[nbb, N]:              Sampled values of the random effects
// ivar_varbb[LT(nbb)]:     (Inverse) covariance matrix of random effects
//                          It is decomposed on exit!
// meanbb[nbb]:             Mean of the random effects
// nbb[1]:                  Dimension of random effects
// N[1]:                    Number of sampled values
// is_varbb[1]:             Indicates whether ivar_varbb contains covariance matrix of random effects 
//
void 
predict_Nrandom_effect_Zero(double *bb,  double *ivar_varbb,  const int *nbb,  const int *N,  const int *is_varbb)
{
  static int i;
  static int err[1];
  static double *bP;

  /*** Decompose the (inverse) covariance matrix ***/
  AK_BLAS_LAPACK::chol_dpptrf(ivar_varbb, nbb, err);
  if (*err){
    throw returnR("Error in predict_common.cpp: predict_Nrandom_effect_Zero. Supplied covariance/precision matrix is not positive definite", 1);
  }

  bP = bb;
  if (*is_varbb){
    for (i = 0; i < *nbb; i++){
      Mvtdist3::rmvnormZero2006(bP, ivar_varbb, nbb);
      bP += *nbb;
    }
  }
  else{
    for (i = 0; i < *nbb; i++){
      Mvtdist3::rmvnormQZero2006(bP, ivar_varbb, nbb);
      bP += *nbb;
    }
  }
  return;
}


void 
predict_Nrandom_effect(double *bb,  double *ivar_varbb,  const double *meanbb, const int *nbb,  const int *N,  const int *is_varbb)
{
  static int i;
  static int err[1];
  static double *bP;

  /*** Decompose the (inverse) covariance matrix ***/
  AK_BLAS_LAPACK::chol_dpptrf(ivar_varbb, nbb, err);
  if (*err){
    throw returnR("Error in predict_common.cpp: predict_Nrandom_effect. Supplied covariance/precision matrix is not positive definite", 1);
  }

  bP = bb;
  if (*is_varbb){
    for (i = 0; i < *nbb; i++){
      Mvtdist3::rmvnorm2006(bP, meanbb, ivar_varbb, nbb);
      bP += *nbb;
    }
  }
  else{
    for (i = 0; i < *nbb; i++){
      Mvtdist3::rmvnormQ2006(bP, meanbb, ivar_varbb, nbb);
      bP += *nbb;
    }
  }
  return;
}


/*** ================================================================================================================ ***/
/*** predict_G1random_effects_Zero:  Predict the values of zero-intercept UNIVARIATE G-spline random effects          ***/
/*** predict_G1random_effects:       Predict the values of general-intercept UNIVARIATE G-spline random effects       ***/
/***                                                                                                                  ***/
/***                                                                                                                  ***/
/*** ================================================================================================================ ***/
//
// bb[N]:                   Sampled values of the random effects
// weight[nknots]:          INPUT:   (Log-)weights (up to a proportionality constant)
//                          OUTPUT:  Cumulative weights (up to a proportionality constant)
// ivar_varbb[1]:           (Inverse) overall variance of random effects (tau^{-2} or tau^{2})
//                          Square root (and inverse) is computed on exit!
// intcptbb[1]:             G-spline intercept
// N[1]:                    Number of sampled values
// is_varbb[1]:             Indicates whether ivar_varbb contains overall variance of G-spline random effects 
//
// knots[nknots]:           Knots (mixture means)
// sigma[nknots]:           Standard deviations of mixture components
// nknots:                  Number of knots
// logw:                    If <> 0 then it is assumed that 'weight' contains log-weights
//
void 
predict_G1random_effect_Zero(double *bb,              double *weight,       double *ivar_varbb,   
                             const int *N,            const int *is_varbb,
                             const double *knots,     const double *sigma,  const int *nknots,    const int *logw)
{
  const double _ZERO_ = 0.0;

  /*** Decompose the (inverse) variance ***/
  if (*ivar_varbb <= 0){
    throw returnR("Error in predict_common.cpp: predict_G1random_effect. Supplied variance/precision not positive", 1);
  }
  if (*is_varbb){
    *ivar_varbb = sqrt(*ivar_varbb);            /*** => tau  ***/
  }
  else{
    *ivar_varbb = 1/sqrt(*ivar_varbb);          /*** => tau  ***/
  }

  Random::rGspline1(bb, weight, N, knots, sigma, nknots, &_ZERO_, ivar_varbb, logw);
 
  return;
}


void 
predict_G1random_effect(double *bb,              double *weight,       double *ivar_varbb,   
                        const double *intcptbb,  const int *N,         const int *is_varbb,
                        const double *knots,     const double *sigma,  const int *nknots,    const int *logw)
{
  /*** Decompose the (inverse) variance ***/
  if (*ivar_varbb <= 0){
    throw returnR("Error in predict_common.cpp: predict_G1random_effect. Supplied variance/precision not positive", 1);
  }
  if (*is_varbb){
    *ivar_varbb = sqrt(*ivar_varbb);            /*** => tau  ***/
  }
  else{
    *ivar_varbb = 1/sqrt(*ivar_varbb);          /*** => tau  ***/
  }

  Random::rGspline1(bb, weight, N, knots, sigma, nknots, intcptbb, ivar_varbb, logw);
 
  return;
}


/*** ================================================================================================================ ***/
/*** predict_GBirandom_effects_Zero:  Predict the values of zero-intercept BIVARIATE G-spline random effects          ***/
/*** predict_GBirandom_effects:       Predict the values of general-intercept BIVARIATE G-spline random effects       ***/
/***                                                                                                                  ***/
/***                                                                                                                  ***/
/*** ================================================================================================================ ***/
//
// bb[2*N]:                 Sampled values of the random effects
// weight[total_length]:   INPUT:   (Log-)weights (up to a proportionality constant)
//                                   Places that correspond to zero weights are filled by zeros and ignored
//                         OUTPUT:  Cumulative weights (up to a proportionality constant)
//                                  for non-zero components on first k_effect places
//
// ind_w_effect[total_length]:   Used only if (is_indfile <> 0)
//                               in that case INPUT:   0/1 indicators of non-zero components
//                                            OUTPUT:  indeces of non-zero components
//                                                     (on scale 0,...,nknots[0]*nknots[1]-1)
//                               
// ivar_varbb[2]:           (Inverse) overall variances of random effects (tau^{-2} or tau^{2})
//                           Square root (and inverse) is computed on exit!
// intcptbb[2]:             G-spline intercepts
//
// N[1]:                    Number of sampled values
// is_varbb[1]:             Indicates whether ivar_varbb contains overall variance of G-spline random effects 
//
// knots0[nknots[0]]:    Knots (mixture means) for margin 0
// knots1[nknots[1]]:    Knots (mixture means) for margin 1
// sigma0[nknots[0]]:    Standard deviations of mixture components for margin 0
// sigma1[nknots[1]]:    Standard deviations of mixture components for margin 1
//
// nknots[2]:            Numbers of knots in each margin
// total_length[1]:      = nknots[0] * nknots[1]
//
// logw:                 If <> 0 then it is assumed that 'weight' contains log-weights
// is_indfile:           If = 0 then it is assumed that k_effect = total_length and 
//                       ind_w_effect is ignored
void 
predict_GBirandom_effect_Zero(double *bb,            double *weight,        int *ind_w_effect,  double *ivar_varbb,
                              const int *N,          const int *is_varbb,
                              const double *knots0,  const double *knots1,  
                              const double *sigma0,  const double *sigma1,  const int *nknots,  const int *total_length,
                              const int *logw,       const int *is_indfile)
{
  const double _ZERO_[2] = {0.0, 0.0};

  /*** Decompose the (inverse) variance ***/
  if (ivar_varbb[0] <= 0 || ivar_varbb[1] <= 0){
    throw returnR("Error in predict_common.cpp: predict_GBirandom_effect. Supplied variance/precision not positive", 1);
  }
  if (*is_varbb){
    ivar_varbb[0] = sqrt(ivar_varbb[0]);            /*** => tau  ***/
    ivar_varbb[1] = sqrt(ivar_varbb[1]);
  }
  else{
    ivar_varbb[0] = 1/sqrt(ivar_varbb[0]);          /*** => tau  ***/
    ivar_varbb[1] = 1/sqrt(ivar_varbb[1]);
  }

  Random::rBiGspline(bb, weight, ind_w_effect, N, knots0, knots1, sigma0, sigma1, nknots, total_length, _ZERO_, ivar_varbb, logw, is_indfile);
 
  return;
}


void 
predict_GBirandom_effect(double *bb,              double *weight,        int *ind_w_effect,  double *ivar_varbb,
                         const double *intcptbb,  const int *N,          const int *is_varbb,
                         const double *knots0,    const double *knots1,  
                         const double *sigma0,    const double *sigma1,  const int *nknots,  const int *total_length,
                         const int *logw,         const int *is_indfile)
{
  /*** Decompose the (inverse) variance ***/
  if (ivar_varbb[0] <= 0 || ivar_varbb[1] <= 0){
    throw returnR("Error in predict_common.cpp: predict_GBirandom_effect. Supplied variance/precision not positive", 1);
  }
  if (*is_varbb){
    ivar_varbb[0] = sqrt(ivar_varbb[0]);            /*** => tau  ***/
    ivar_varbb[1] = sqrt(ivar_varbb[1]);
  }
  else{
    ivar_varbb[0] = 1/sqrt(ivar_varbb[0]);          /*** => tau  ***/
    ivar_varbb[1] = 1/sqrt(ivar_varbb[1]);
  }

  Random::rBiGspline(bb, weight, ind_w_effect, N, knots0, knots1, sigma0, sigma1, nknots, total_length, intcptbb, ivar_varbb, logw, is_indfile);
 
  return;
}


/*** ====================================================================== ***/
/*** print_iter_info2:  Print iteration information                         ***/
/***                                                                        ***/
/*** ====================================================================== ***/
void
print_iter_info2(int &backs,  const int &iter,  const int &nwrite,  const int &lastIter)
{
  static int i;

  if (!((iter+1) % nwrite) || (iter+1) == lastIter){
    for (i = 0; i < backs; i++) Rprintf("\b");
    Rprintf("%d", iter+1);
    backs = int(log10(double(iter+1))) + 1;
  }
}

}    /*** end of the namespace predict_common ***/
