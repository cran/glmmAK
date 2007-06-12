/*** Quantile2.h ***/

#ifndef _QUANTILE_TWEE_H_
#define _QUANTILE_TWEE_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"

namespace Quantile{

extern "C"{

  void
  Quantile(double *quantile,
           const double *sample,  const int *ngrid,  const int *nsample,
           const double *prob,    const int *nprob);

}

}  /*** end of the namespace Quantile ***/

#endif
