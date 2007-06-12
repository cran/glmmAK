/*** predict_poisson.h ***/

#ifndef _PREDICT_POISSON_H_
#define _PREDICT_POISSON_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_BLAS_LAPACK.h"

#include "util_poisson.h"
#include "predict_common.h"
#include "mcmc_Random.h"
#include "In_Output2.h"
#include "Quantile2.h"
#include "summary_Gspline1.h"
#include "summary_BiGspline.h"

namespace predict_poissonA{

extern "C"{

void
predict_poisson(double *ecount_value,             double *ecount_summary,                
                const double *qprob,              const int *nqprob,              const int *retValue,
                const int *n,                     const int *N,                   const int *ni,
                const double *offset,             const double *X,                const int *p,
                const double *XRE,                const int *pRE,
                const int *REdist,                const int *hierarCenter,
                const int *G_dim,                 const double *G_dPar,
                const double *betaF,              const double *betaR,            double *ivarR_varR,             const int *is_varR,
                const int *niter,                 
                char **wpath,                     char **indpath,                 const int *skip,                const int *by,
                const int *logw,                  const int *is_indfile,          const int *nwrite,              int *err);

}    /*** end of extern "C" ***/

void
print_iter_infoP2(int &backs,  const int &iter,  const int &nwrite,  const int &lastIter);

}    /*** end of namespace predict_cumlogitA ***/

#endif
