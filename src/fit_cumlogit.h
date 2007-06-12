/*** fit_cumlogit.h ***/

#ifndef _FIT_CUMLOGIT_H_
#define _FIT_CUMLOGIT_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_Optim.h"

#include "MatrixLT.h"
#include "MatrixRect.h"

#include "ll_cumlogit.h"
#include "util_cumlogit.h"

namespace Fit_cumlogit{

extern "C"{

  const int _AK_MAX_STEPHALF = 10;         /*** Maximal number of step-halfing steps ***/
   
  const int _AK_MAX_PD_ATTEMPTS = 10;      /*** Maximal number of attempts to make the minus Hessian positive definite ***/
  const double _AK_EPS_PD_ATTEMPT = 0.01;  /*** Epsilon for 'chol_solvePD' function                                    ***/

  /* Types of the optimization algorithm, 0 = NR = Newton-Raphson, 1 = FS = Fisher scoring */
  enum optimAlgorithm {NR, FS};  

  void
  fit_cumlogit(double *theta,
               double *ll,          double *U,           double *I_obs,      double *I_exp,
               const int *Y,        const double *X,     const double *V,  
               const int *C,        const int *p,        const int *q,       const int *n,
               int* niter,          const double *toler, const int *trace,
               int *err);

}  /* end of extern "C" */

}  /* end of the namespace Fit_cumlogit */

#endif
