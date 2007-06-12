/*** util_BiGspline2.h ***/

#ifndef _UTIL_BI_G_SPLINE_TWEE_H_
#define _UTIL_BI_G_SPLINE_TWEE_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

#include "MatrixRect.h"
#include "BiGspline2.h"
#include "Random2.h"

namespace util_BiGspline2{

  void
  array2allocBi(MatrixRect<int> *allocM,  const int *allocA,  const int *K,  const int &nObs);

  void
  alloc2arrayBi(int *allocA,  const MatrixRect<int> *allocM,  const int *K,  const int &nObs);

  void
  alloc2mixtureNBi(int *mixtureN,  const int *allocS,  const int &total_length,  const int &nObs);

  void
  updateAllocBi(int *allocS,        int *mixtureN,       BiGspline2 *Gspl,
                const double *Obs,  const double *EObs,  const int &nObs);

}

#endif
