/*** Random2.h ***/

#ifndef _RANDOM_TWEE_H_
#define _RANDOM_TWEE_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

namespace Random {

void
rltruncGamma(double* x,  
             const double* shape,  const double* rate,   const double* minx,
             const int* n);

void
discreteSampler2(int* sampledj,  double* propA,  const int* kP,  const int* nP,  const int* cumul);

int 
findUniformIndex(const double u,  const int startInd,  const int endInd,  const int k);

int 
findIndex(const double u,  const int startInd,  const int endInd,  const double* ValuesA);

void
rGspline1(double *x,             double *weight, 
          const int *nx,
          const double *knots,   const double *sigma,  const int *nknots,
          const double *intcpt,  const double *tau,    const int *logw);

extern "C"{
void
rGspline1R(double *x,             double *weight, 
           const int *nx,
           const double *knots,   const double *sigma,  const int *nknots,
           const double *intcpt,  const double *tau,    const int *logw);
}

void
rBiGspline(double *x,              double *weight,             int *ind_w_effect,
           const int *nx,
           const double *knots0,   const double *knots1,
           const double *sigma0,   const double *sigma1,       const int *nknots,   const int *total_length,
           const double *intcpt,   const double *tau,      
           const int *logw,        const int *is_indfile);

extern "C"{
void
rBiGsplineR(double *x,              double *weight,             int *ind_w_effect,
            const int *nx,
            const double *knots0,   const double *knots1,
            const double *sigma0,   const double *sigma1,       const int *nknots,   const int *total_length,
            const double *intcpt,   const double *tau,      
            const int *logw,        const int *is_indfile);
}

}   /*** end of namespace Random ***/

#endif
