/*** util_Gspline2.h ***/

#ifndef _UTIL_G_SPLINE_TWEE_H_
#define _UTIL_G_SPLINE_TWEE_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

#include "MatrixRect.h"
#include "Gspline2.h"
#include "Random2.h"

namespace util_Gspline2{

extern "C"{

void
max_poster_prob(int *Alloc,               double *BasisStdDev,   double *Knots,      const int *K,
                const double *Data,       const int *Dim,        const int *nData,
                const double *Intcpt,     const double *StdDev);

}    /* end of extern "C" */

void
array2alloc(MatrixRect<int> *allocM,  const int *allocA,  const int *K,  const int &dim,  const int &nObs);

void
alloc2array(int *allocA,  const MatrixRect<int> *allocM,  const int *K,  const int &dim,  const int &nObs);

void
alloc2mixtureN(MatrixRect<int> *mixtureN,  const int *alloc,  const int *K,  const int &dim,  const int &nObs);

void
updateAlloc1(int *Alloc,         MatrixRect<int> *MixtureN,  Gspline2 *Gspl,
             const double *Obs,  const double *EObs,         const int &nObs);

}

#endif
