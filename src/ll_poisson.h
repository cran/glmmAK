/*** ll_poisson.h ***/

#ifndef _LL_POISSON_H_
#define _LL_POISSON_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

#include "MatrixLT.h"
#include "AK_BLAS_LAPACK.h"

namespace Ll_poisson{

void
ll_poisson(double *ll,
           double *U,
           double *I,
           double *eta,
           double *mu,
           const double *offset,
           const double *theta,
           const int *y,
           const double *log_y_factor,
           const double *x,
           const MatrixLT<double> *xx,
           const int &nObs,
           const int &p,
           const int &order);

}  /** end of the namespace Ll_poisson **/

#endif
