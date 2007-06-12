/*** Mvtdist2.h ***/

#ifndef _MULTI_VARIATE_DISTRIBUTIONS_H_
#define _MULTI_VARIATE_DISTRIBUTIONS_H_

#include <cmath>

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

#include "MatrixRect.h"
#include "MatrixLT.h"

namespace Mvtdist2{

double
dmvtnorm2(MatrixRect<double> &x_m,  const MatrixLT<double> &Li,  const int &dlog=0,  const int &fullDens=1);

double
dmvtnorm2a(MatrixRect<double> &x_m,  const double *x,  const double *mu,  const MatrixLT<double> *Li,  const int &dlog=0);

double
dmvtnorm2b(const double *x,  const double *mu,  const double *invVar,  const int &dim,  const int &dlog=0);

}

#endif
