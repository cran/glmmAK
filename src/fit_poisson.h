/*** fit_poisson.h ***/

#ifndef _FIT_POISSON_H_
#define _FIT_POISSON_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_Optim.h"

#include "MatrixLT.h"
#include "AK_BLAS_LAPACK.h"

#include "ll_poisson.h"

namespace Fit_poisson{

extern "C"{

  const int _AK_MAX_STEPHALF = 10;         /*** Maximal number of step-halfing steps ***/
   
  const int _AK_MAX_PD_ATTEMPTS = 10;      /*** Maximal number of attempts to make the minus Hessian positive definite ***/
  const double _AK_EPS_PD_ATTEMPT = 0.01;  /*** Epsilon for 'chol_solvePD' function                                    ***/

  void
  fit_poisson(double *theta,       double *eta,           double *mu,
              double *ll,          double *U,             double *I,
              const int *Y,        const double *offset,  const double *X,
              const int *p,        const int *n,
              int* niter,          const double *toler,   const int *trace,
              int *err);
}

}  /* end of the namespace Fit_poisson */

#endif
