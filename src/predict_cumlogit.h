/*** predict_cumlogit.h ***/

#ifndef _PREDICT_CUMLOGIT_H_
#define _PREDICT_CUMLOGIT_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_BLAS_LAPACK.h"

#include "util_cumlogit.h"
#include "predict_common.h"
#include "mcmc_Random.h"
#include "In_Output2.h"
#include "Quantile2.h"
#include "summary_Gspline1.h"
#include "summary_BiGspline.h"

namespace predict_cumlogitA{

extern "C"{

void
predict_cumlogit(double *prob_value,               double *prob_summary,                
                 const double *qprob,              const int *nqprob,              const int *retValue,
                 const int *C,                     const int *n,                   const int *N,                   const int *ni,
                 const double *X,                  const double *V,                const int *p,                   const int *q,
                 const double *XRE,                const double *VRE,              const int *pRE,                 const int *qRE,
                 const int *REdist,                const int *hierarCenter,
                 const int *G_dim,                 const double *G_dPar,
                 const double *betaF,              const double *betaR,            double *ivarR_varR,             const int *is_varR,
                 const int *niter,                 
                 char **wpath,                     char **indpath,                 const int *skip,                const int *by,
                 const int *logw,                  const int *is_indfile,          const int *nwrite,              int *err);

}    /*** end of extern "C" ***/

}    /*** end of namespace predict_cumlogitA ***/

#endif
