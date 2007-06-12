/*** GMRF_Gspline.h ***/

#ifndef _GAUSSIAN_MARKOV_RANDOM_FIELD_AND_G_SPLINE_H_
#define _GAUSSIAN_MARKOV_RANDOM_FIELD_AND_G_SPLINE_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_BLAS_LAPACK.h"
#include "GMRF.h"
#include "GMRF_Gspline_Util.h"

namespace GMRF_Gspline {

const int    _maxiter      = 1;          /** Maximal number of Newton-Raphson iterations when constructing the normal approximation   **/
const int    _max_nAttempt = 10;         /** Maximal number of attempts to make the Hessian PD by adding epsw to diagonal             **/
const double _toler        = 1e-3;       /** Tolerance to detect convergence of the Newton-Raphson                                    **/
const int    _max_stephalf = 10;         /** Maximal number of stephalving steps within Newton-Raphson                                **/

const double _null_mass = 1e-6;          /** Constant added to the diagonal of the Hessian when trying to make it PD                  **/


void
update(int *accept,               double *a,                    double *lambda,  
       double *expa,              double *sumexpa,              double *w,                  double *minw,
       double *Da,                double *Qa,                   double *min_half_aQa,   
       double *workML,            double *worka,                double *workGMRF,
       const int *allocN,         const int *prior_for_lambda,  const double *par_lambda,   const double* par_rscale,
       const double *Q,           const int *order,             const int *diffOper,        const double *epsw,
       const int *constraint,     const int *iref,              const int *na,              const int *nobs,
       const int *lambda_a_block);

void
ML_est(double *ll,             double *dll,           double *ddll,
       double *a,              double *workD,         int *niter,          int *err,                 
       const int *allocN,      const double *lambda,  const double *Q,     const int *order,    const int *diffOper,
       const int *constraint,  const int *iref,
       const int *na,          const int *nobs,       const int *maxiter,  const double *epsw);

void
ll0(double *ll,             
    const double *a,        const int *allocN,           const double *lambda,
    const double *sumexpa,  const double *min_half_aQa,  
    const int *na,          const int *nobs);

void
ll1(double *ll,             double *dll,
    const double *a,        const int *allocN,           const double *lambda,
    const double *sumexpa,  const double *min_half_aQa,  const double *Qa,      const double *w,
    const int *constraint,  const int *iref,        
    const int *na,          const int *nobs);

void
ll2(double *ll,             double *dll,                 double *ddll,          double *workll2,
    const double *a,        const int *allocN,           const double *lambda,
    const double *sumexpa,  const double *min_half_aQa,  const double *Qa,      const double *w,
    const double *Q,        const int *order,
    const int *constraint,  const int *iref,
    const int *na,          const int *nobs);

extern "C"{
  void
  mcmc_GMRF_Gspline(int *acceptSample,   double *aSample,              double *wSample,   double *lambdaSample,  
                    int *iter,
                    const int *allocN,   const int *prior_for_lambda,  const double *par_lambda,  const double *F,  
                    const int *order,    const int *constraint,        const int *iref,
                    const int *na,       const int *nobs,              const int *lambda_a_block,
                    const int *nsimul);
}

}  /*** end of namespace GMRF_Gspline ***/

#endif
