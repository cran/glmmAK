/*** GMRF_Gspline.h ***/

#ifndef _GAUSSIAN_MARKOV_RANDOM_FIELD_AND_G_SPLINE_UTILITIES_H_
#define _GAUSSIAN_MARKOV_RANDOM_FIELD_AND_G_SPLINE_UTILITIES_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_BLAS_LAPACK.h"
#include "GMRF.h"

namespace GMRF_Gspline_Util {

enum aIdentifiabilityConstraint {_Mean_, _Reference_, _General_};
enum priorForLambdaTypes {_Fixed_, _Gamma_, _SDUnif_};                         // types of priors for lambda
enum aaType {_a2a_, _a2d_, _d2a_, _d2d_};

const double _null_weight = 1e-6;             /** all weights being smaller than this will be considered as being equal to _null_weight **/
                                              /**  changes are performed by 'update4_ll12' function                                     **/

extern "C"{
  void
  d2a(double *coef,  const int *constraint,  const int *iref,  const int *na);

  void
  d2a2(double *a,  const double *d,  const int *constraint,  const int *iref,  const int *na);
}

extern "C"{
  void
  a2d(double *coef,  const int *iref,  const int *na);

  void
  a2d2(double *d,  const double * a,  const int *iref,  const int *na);
}

void
copy_within_update
  (double *expaTo,          double *sumexpaTo,          double *wTo,                     double *minwTo,  
   double *DaTo,            double *QaTo,               double *min_half_aQaTo,          double *aTo,
   const double *expaFrom,  const double *sumexpaFrom,  const double *wFrom,             const double *minwFrom,  
   const double *DaFrom,    const double *QaFrom,       const double *min_half_aQaFrom,  const double *aFrom,
   const int *constraint,   const int *iref,            const int *na,                   const int *na_1,
   const int &aatype);


void
NRstep(double *a,  const double *NR_step,  const int *constraint,  const int *iref,  const int *na);

void
NRstephalve(double *a,  double *NR_step,  const int *constraint,  const int *iref,  const int *na);

void
update4_ll0(double *expa,     double *sumexpa,   double *Da,    double *min_half_aQa,
            const double *a,  const int *order,  const int *na);

void
update4_ll12(double *expa,     double *sumexpa,   double *Da,           double *min_half_aQa,
             double *Qa,       double *w,         double *minw,
             const double *a,  const int *order,  const int *diffOper,  const int *na);


void
rltruncGamma(double *x,  const double *shape,  const double *scale,   const double *minx);

}    /*** end of namespace GMRF_Gspline_Util ***/

#endif

