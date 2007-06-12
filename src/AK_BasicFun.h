/*** AK_BasicFun.h ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  <16/03/2005 (as constants.h)
// REVISION 1:  04/08/2006
//
// PURPOSE: Some useful constants and basic functions
//
/* ********************************************************************************* */

#ifndef _AK_BASIC_FUN_H_
#define _AK_BASIC_FUN_H_


const int _AK_INT_MAX = 999999;
const int _AK_ONE_INT = 1;
const int _AK_ZERO_INT = 0;

const double _AK_LOG_SQRT_2PI = 0.918938533204672741780329736406;	     // log(sqrt(2*pi)) 
const double _AK_NORM_ZERO = 1e-16;                                          // qnorm(1 - 1e-16) is still non infty in R
const double _AK_ZERO = 1e-50;
const double _AK_invFLT_MAX = 1e-50;
const double _AK_LOG_2 = 0.6931472;
const double _AK_SCALE_ZERO = 1e-20;
const double _AK_invSQRT_TWO_PI = 0.39894228040143270286;
const double _AK_invTWO_PI = 0.1591549430918953;

const double _AK_EMIN = -115;                                          // exp(-115) will be 0 (just > 1e-50)
const double _AK_EMAX = 115;                                           // exp(115) will be Infty (just < 1e50)
const double _AK_expEMIN = exp(_AK_EMIN);

const double _AK_TOL_CHOL = 1e-10;           // tolerance for the Cholesky decomposition

/* exp(x) */
inline double
exp_AK(const double &x){
  return (x < _AK_EMIN ? 0.0 : (x > _AK_EMAX ? R_PosInf : exp(x)));
};

/* exp(x)/(1 + exp(x)) */
inline double
invlogit_AK(const double &x){
  const double exp_x=exp(x);
  return (x < _AK_EMIN ? 0.0 : (x > _AK_EMAX ? 1.0 : exp_x/(1 + exp_x)));
};

/* log(x) */
inline double
log_AK(const double &x){
  return(x < _AK_ZERO ? R_NegInf : log(x));
}

/* max(vector) */
inline int
maxArray(const int *x, const int &nx){
  int VALUE = *x;
  const int *xP = x;
  for (int i = 1; i < nx; i++){
    xP++;
    if (*xP > VALUE) VALUE = *xP;
  }
  return(VALUE);
}

inline double
maxArray(const double *x, const int &nx){
  double VALUE = *x;
  const double *xP = x;
  for (int i = 1; i < nx; i++){
    xP++;
    if (*xP > VALUE) VALUE = *xP;
  }
  return(VALUE);
}

#endif
