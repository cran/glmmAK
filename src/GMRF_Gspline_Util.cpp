/*** GMRF_Gspline_Util.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  22/11/2006
//
//                    d2a:  22/11/2006
//                   d2a2:  22/11/2006
//                    a2d:  22/11/2006
//                   a2d2:  22/11/2006
//     copy_within_update:  17/11/2006
//
//                 NRstep:  22/11/2006
//            NRstephalve:  22/11/2006
//
//            update4_ll0:  08/11/2006
//           update4_ll12:  08/11/2006
//
//           rltruncGamma:  29/11/2006  (as (almost) copy of rltruncGamma from Random.cpp, resp. Random2.cpp)
//                          
//
// PURPOSE: Functions that implement utilities for (intrinsic) Gaussian random fields
//          in the G-spline context
//
/* ********************************************************* */

#include "GMRF_Gspline_Util.h"

namespace GMRF_Gspline_Util {

/*** ============================================================================= ***/
/*** Create a(d) from d                                                            ***/
/*** *** 2 versions ***                                                            ***/
/*** ============================================================================= ***/
//
// coef[na]:  INPUT:   d coefficients on first (na-1) places
//            OUTPUT:  a(d) coefficients
// a[na]:     a(d) coefficients
// d[na-1]:   d coefficients
//
extern "C"{
  void
  d2a(double *coef,  const int *constraint,  const int *iref,  const int *na)
  {
    static int j;
    static double *coefP, *arefP;
    static double asum;

    j = 0;
    coefP = coef;

    switch (*constraint){
    case _Mean_:
      asum = 0;
      while (j < *iref){
        asum += *coefP;
        coefP++;
        j++;
      }
     
      arefP = coefP;

      j = *na - 1;
      coefP += j - (*iref);
      while (j > *iref){
        *coefP = *(coefP - 1);
        asum += *coefP;
        coefP--;
        j--;
      }

      *arefP = -asum;
      return;

    case _Reference_:
      while (j < *iref){
        coefP++;
        j++;
      }
     
      arefP = coefP;

      j = *na - 1;
      coefP += j - (*iref);
      while (j > *iref){
        *coefP = *(coefP - 1);
        coefP--;
        j--;
      }

      *arefP = 0;
      return;

    default:
      REprintf("constraint = %d\n", *constraint);
      throw returnR("Error in GMRF_Gspline_Util.cpp: d2a(). Unknown value of the identifiability constraint.", 1);
    }
  }


  void
  d2a2(double *a,  const double *d,  const int *constraint,  const int *iref,  const int *na)
  {
    static int j;
    static double *aP, *arefP;
    static double asum;
    static const double *dP;

    j = 0;
    aP = a;
    dP = d;

    switch (*constraint){
    case _Mean_:
      asum = 0;
      while (j < *iref){
        *aP = *dP;
        asum += *aP;
        aP++;
        dP++;
        j++;
      }
     
      arefP = aP;
      aP++;
      j++;

      while (j < *na){
        *aP = *dP;
        asum += *aP;
        aP++;
        dP++;
        j++;
      }

      *arefP = -asum;
      return;

    case _Reference_:
      while (j < *iref){
        *aP = *dP;
        aP++;
        dP++;
        j++;
      }
     
      *aP = 0;
      aP++;
      j++;

      while (j < *na){
        *aP = *dP;
        aP++;
        dP++;
        j++;
      }
      return;

    default:
      REprintf("constraint = %d\n", *constraint);
      throw returnR("Error in GMRF_Gspline_Util.cpp: d2a2(). Unknown value of the identifiability constraint.", 1);
    }
  }
}


/*** ============================================================================= ***/
/*** Create d from a (omit reference a from the vector)                            ***/
/*** ============================================================================= ***/
//
// coef[na]:  INPUT:   a(d) coefficients
//            OUTPUT:  d coefficients on first (na-1) places
//
extern "C"{
  void
  a2d(double *coef,  const int *iref,  const int *na)
  {
    static int j;
    static double *coefP;

    j = 0;
    coefP = coef;
    while (j < *iref){
      coefP++;
      j++;
    }
    j++;
    while (j < *na){
      *coefP = *(coefP+1);
      coefP++;
      j++;
    }

    return;
  }

  void
  a2d2(double *d,  const double * a,  const int *iref,  const int *na)
  {
    static int j;
    static double *dP;
    static const double *aP;

    j = 0;
    aP = a;
    dP = d;
    while (j < *iref){
      *dP = *aP;
      dP++;
      aP++;
      j++;
    }
    aP++;
    j++;
    while (j < *na){
      *dP = *aP;
      dP++;
      aP++;
      j++;
    }

    return;
  }
}


/*** ============================================================================= ***/
/*** Sub-part of update_GMRF_Gspline                                               ***/
/***   copy some arrays                                                            ***/
/*** ============================================================================= ***/
//
// aTo[na]:
// aFrom[na]:
// na_1:       na - 1
// aatype:     Type of transformation aFrom -> aTo (see GMRF_Gspline.h for union of possible values)
//             0 = a(d) -> a(d)
//             1 = a(d) -> d
//             2 = d -> a(d)
//             3 = d -> d
//
void
copy_within_update
  (double *expaTo,          double *sumexpaTo,          double *wTo,                     double *minwTo,  
   double *DaTo,            double *QaTo,               double *min_half_aQaTo,          double *aTo,
   const double *expaFrom,  const double *sumexpaFrom,  const double *wFrom,             const double *minwFrom,  
   const double *DaFrom,    const double *QaFrom,       const double *min_half_aQaFrom,  const double *aFrom,
   const int *constraint,   const int *iref,            const int *na,                   const int *na_1,
   const int &aatype)
{
  AK_BLAS_LAPACK::copyArray(expaTo, expaFrom, *na);
  *sumexpaTo = *sumexpaFrom;
  AK_BLAS_LAPACK::copyArray(wTo, wFrom, *na);
  *minwTo = *minwFrom;
  AK_BLAS_LAPACK::copyArray(DaTo, DaFrom, *na);
  AK_BLAS_LAPACK::copyArray(QaTo, QaFrom, *na);
  *min_half_aQaTo = *min_half_aQaFrom;

  switch (aatype){
  case _a2a_:
    AK_BLAS_LAPACK::copyArray(aTo, aFrom, *na);
    return;

  case _a2d_:
    a2d2(aTo, aFrom, iref, na);
    return;

  case _d2a_:
    d2a2(aTo, aFrom, constraint, iref, na);
    return;

  case _d2d_:
    AK_BLAS_LAPACK::copyArray(aTo, aFrom, *na_1);
    return;

  default:
    REprintf("aatype = %d\n", aatype);
    throw returnR("Error in GMRF_Gspline_Util.cpp: copy_within_update(). Unknown value of aatype argument.", 1);
  }
}


/*** ============================================================================= ***/
/*** Sub-part of ML_est_GMRF_Gspline                                               ***/
/***   Update a coefficient by adding Newton-Raphson step                          ***/
/*** ============================================================================= ***/
//
//  a[na]
//  NR_step[na-1]
//
void
NRstep(double *a,  const double *NR_step,  const int *constraint,  const int *iref,  const int *na)
{
  static int j;
  static double asum;
  static double *aP, *arefP;
  static const double *NR_stepP;

  aP       = a;
  NR_stepP = NR_step;
  j = 0;

  switch (*constraint){
  case _Mean_:
    asum = 0;
    while (j < *iref){
      *aP += *NR_stepP;
      asum += *aP;
      aP++;
      NR_stepP++;
      j++;
    }

    arefP = aP;
    aP++;
    j++;

    while (j < *na){
      *aP += *NR_stepP;
      asum += *aP;
      aP++;
      NR_stepP++;
      j++;
    }

    *arefP = -asum;
    return;

  case _Reference_:
    while (j < *iref){
      *aP += *NR_stepP;
      aP++;
      NR_stepP++;
      j++;
    }

    *aP = 0;
    aP++;
    j++;

    while (j < *na){
      *aP += *NR_stepP;
      aP++;
      NR_stepP++;
      j++;
    }
    return;

  default:
    REprintf("constraint = %d\n", *constraint);
    throw returnR("Error in GMRF_Gspline_Util.cpp: NRstep(). Unknown value of the identifiability constraint.", 1);
  }
}


/*** ============================================================================= ***/
/*** Sub-part of ML_est_GMRF_Gspline                                               ***/
/***   Update a coefficient by subtracting 0.5*Newton-Raphson step                 ***/
/***   Multiply Newton-Raphson step by 0.5                                         ***/
/*** ============================================================================= ***/
//
//  a[na]
//  NR_step[na-1]
//
void
NRstephalve(double *a,  double *NR_step,  const int *constraint,  const int *iref,  const int *na)
{
  static int j;
  static double asum;
  static double *aP, *arefP, *NR_stepP;

  aP       = a;
  NR_stepP = NR_step;
  j = 0;

  switch (*constraint){
  case _Mean_:
    asum = 0;
    while (j < *iref){
      *NR_stepP *= 0.5;
      *aP -= *NR_stepP;
      asum += *aP;
      aP++;
      NR_stepP++;
      j++;
    }

    arefP = aP;
    aP++;
    j++;

    while (j < *na){
      *NR_stepP *= 0.5;
      *aP -= *NR_stepP;
      asum += *aP;
      aP++;
      NR_stepP++;
      j++;
    }

    *arefP = -asum;
    return;

  case _Reference_:
    while (j < *iref){
      *NR_stepP *= 0.5;
      *aP -= *NR_stepP;
      aP++;
      NR_stepP++;
      j++;
    }

    *aP = 0;
    aP++;
    j++;

    while (j < *na){
      *NR_stepP *= 0.5;
      *aP -= *NR_stepP;
      aP++;
      NR_stepP++;
      j++;
    }
    return;

  default:
    REprintf("constraint = %d\n", *constraint);
    throw returnR("Error in GMRF_Gspline_Util.cpp: NRstep(). Unknown value of the identifiability constraint.", 1);
  }
}


/*** ============================================================================= ***/
/*** Functions that update some quantities needed to compute log-likelihood        ***/
/***  or its derivatives                                                           ***/
/***                                                                               ***/
/***   FUNCTIONS: update4_ll0, update4_ll12                                        ***/
/***                                                                               ***/
/*** ============================================================================= ***/
//
// expa[na]:           array with exp(a[j])
// sumexpa[1]:         sum(exp(a[j]))
// Da[na]:             difference of order 'order' in the first (na - order) places
//                     for working purposes, array of length 'na' is needed
// min_half_aQa[1]:  
// a[na]:
// order[1]:           order of the differences  
// diffOper[order+1]:  vector defining the difference operator, e.g.,
//                        order = 0:   1
//                        order = 1:   -1, 1
//                        order = 2:   1, -2, 1
//                        order = 3:   -1, 3, -3, 1
//                         etc.
// na[1]:              length of a
//
void
update4_ll0(double *expa,     double *sumexpa,   double *Da,           double *min_half_aQa,
            const double *a,  const int *order,  const int *na)
{
  static int j;
  static double *expaP, *DaP;
  static const double *aP;

  /*** expa, sumexpa, copy a to Da ***/
  expaP    = expa;
  DaP      = Da;
  aP       = a;
  *sumexpa = 0.0;
  for (j = 0; j < *na; j++){
    *expaP = exp_AK(*aP);
    *sumexpa += *expaP;
    *DaP = *aP;

    aP++;
    expaP++;
    DaP++;
  }

  /*** Da ***/
  GMRF::diff(Da, order, na);

  /*** min_half_aQa ***/
  *min_half_aQa = 0.0;
  DaP = Da;
  for (j = 0; j < *na - (*order); j++){
    *min_half_aQa += (*DaP) * (*DaP);
    DaP++;
  }
  *min_half_aQa *= -0.5;

  return;
}


void
update4_ll12(double *expa,     double *sumexpa,   double *Da,           double *min_half_aQa,
             double *Qa,       double *w,         double *minw,
             const double *a,  const int *order,  const int *diffOper,  const int *na)
{
  static int j;
  static double *wP, *expaP;

  /*** expa, sumexpa, Da, min_half_aQa ***/
  update4_ll0(expa, sumexpa, Da, min_half_aQa, a, order, na);

  /*** w, minw ***/
  *minw = 1.0;
  wP    = w;
  expaP = expa;
  for (j = 0; j < *na; j++){
    *wP = *expaP/(*sumexpa);
    if (*wP < GMRF_Gspline_Util::_null_weight) *wP = _null_weight;
    if (*wP < *minw) *minw = *wP;
    wP++;
    expaP++;
  }

  /*** Qa = t(D)*D*a ***/
  GMRF::tdiff(Qa, Da, diffOper, order, na);

  return;
}


/*** ============================================================================================= ***/
/*** rltruncGamma:  Sample from left-truncated gamma distribution                                  ***/
/*                                                                                                   */
/* PARAMETERS:                                                                                       */
/*   shape .......... shape of the gamma distribution                                                */
/*   scale .......... scale (= 1/rate) of the gamma distribution                                     */
/*   minx ........... a value > 0 indicating where the gamma distribution is to be truncated         */
/*                                                                                                   */
/*** ============================================================================================= ***/
void
rltruncGamma(double *x,  const double *shape,  const double *scale,   const double *minx)
{
  double u;  

  /* Rprintf("minx = %e, shape = %e, scale = %e\n", *minx, *shape, *scale);  */
  double Flower = pgamma(*minx, *shape, *scale, 1, 0);
  if (Flower >= 1 - _AK_NORM_ZERO){   // truncation time irrealistic large => sampled values = minx
    *x = *minx;
  }
  else{
    if (Flower <= _AK_NORM_ZERO){     // truncation time = 0, sample from gamma distribution
      *x = rgamma(*shape, *scale);
    }
    else{
      u = runif(0, 1) * (1 - Flower) + Flower;
      *x = qgamma(u, *shape, *scale, 1, 0);
    }
  }

  return;
}


}    /*** end of namespace GMRF_Gspline_Util ***/
