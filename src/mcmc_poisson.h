/*** mcmc_poisson.h ***/

#ifndef _MCMC_POISSON_H_
#define _MCMC_POISSON_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "In_Output2.h"
#include "TemplateFun2.h"

#include "AK_BLAS_LAPACK.h"

#include "MatrixLT.h"
#include "MatrixRect.h"
#include "PredictorPoiss.h"
#include "RandomPoiss.h"
#include "Gspline2.h"
#include "BiGspline2.h"
#include "mcmc_Random.h"
#include "util_Gspline2.h"
#include "util_BiGspline2.h"

namespace mcmc_poissonA{

extern "C"{

void
mcmc_poisson(char **dirP,
             const int *Y,                 const double *offset,           const int *n,                   const int *N,              const int *ni,
             const double *X,              const int *p,
             double *Beta,                 const double *priorBeta,
             const double *XRE,            const int *pRE,
             double *RE,                   const int *REdist,              const int *hierarCenter,
             double *ParRE,                const int *priorREMean_InvVar,  const double *priorParREMean,   const double *priorParREInvVar,
             double *G_lambda_a,           const int *G_dim,               const int *G_ipar,              const double *G_dpar,      
             const double *G_lambdaPrior,  const int *G_apar,              const double *G_aContrast,       
             int *allocRE,
             int *iter,                    const int *nsimul,              const int *store,               const int *mainSimul,
             const int *precision,         int *err);

void
print_iter_infoP(int &writeAll,  int &backs,  const int &iter,  const int &nwrite,  const int &lastIter);

}    /* end of extern "C" */

}  /*** end of the namespace mcmc_poissonA ***/

#endif

