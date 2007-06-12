/*** predict_common.h ***/

#ifndef _PREDICT_COMMON_H_
#define _PREDICT_COMMON_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_BLAS_LAPACK.h"

#include "RandomCL.h"
#include "Mvtdist3.h"

namespace predict_common{

void 
predict_Nrandom_effect_Zero(double *bb,  double *ivar_varbb,  const int *nbb,  const int *N,  const int *is_varbb);

void 
predict_Nrandom_effect(double *bb,  double *ivar_varbb,  const double *meanbb, const int *nbb,  const int *N,  const int *is_varbb);

void 
predict_G1random_effect_Zero(double *bb,              double *weight,       double *ivar_varbb,   
                             const int *N,            const int *is_varbb,
                             const double *knots,     const double *sigma,  const int *nknots,    const int *logw);

void 
predict_G1random_effect(double *bb,              double *weight,       double *ivar_varbb,   
                        const double *intcptbb,  const int *N,         const int *is_varbb,
                        const double *knots,     const double *sigma,  const int *nknots,    const int *logw);

void 
predict_GBirandom_effect_Zero(double *bb,            double *weight,        int *ind_w_effect,  double *ivar_varbb,
                              const int *N,          const int *is_varbb,
                              const double *knots0,  const double *knots1,  
                              const double *sigma0,  const double *sigma1,  const int *nknots,  const int *total_length,
                              const int *logw,       const int *is_indfile);

void 
predict_GBirandom_effect(double *bb,              double *weight,        int *ind_w_effect,  double *ivar_varbb,
                         const double *intcptbb,  const int *N,          const int *is_varbb,
                         const double *knots0,    const double *knots1,  
                         const double *sigma0,    const double *sigma1,  const int *nknots,  const int *total_length,
                         const int *logw,         const int *is_indfile);

void
print_iter_info2(int &backs,  const int &iter,  const int &nwrite,  const int &lastIter);

}    /*** end of the namespace predict_common ***/

#endif
