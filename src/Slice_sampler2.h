/*** Slice_sampler2.h ***/

#ifndef _SLICE_SAMPLER_H_
#define _SLICE_SAMPLER_H_

#include <R.h>
#include <Rmath.h>

#include <iostream>
#include <climits>
#include <cmath>

#include "AK_Error.h"

#ifdef __cplusplus
extern "C" {
#endif

namespace Slice_sampler{

void
ss_stepping_out(double* interv,       double* g_interv,
                const double* x0,     const double* z,     const double* w,  const int* m,
                const double* bound,  const int* is_bound, 
                void (*eval)(const double*, double*, const double*, const int*),
                const double* deval,  const int* ieval);

void
ss_doubling(double* interv,       double* g_interv,
            const double* x0,     const double* z,      const double* w,     const int* p,
            const double* bound,  const int* is_bound,  const int* unimodal,
            void (*eval)(const double*, double*, const double*, const int*),
            const double* deval,  const int* ieval);

void
ss_accept_doubling(int* accept,           const double* x0,        const double* x1,     
                   const double* z,       const double* w,         const double* interv,
                   void (*eval)(const double*, double*, const double*, const int*),
                   const double* deval,   const int* ieval);

void
ss_shrinkage_sample(double* x1,           double* interv,       double* g_interv,     const double* x0,     
                    const double* z,      const double* w,      const int* doubling,  const int* unimodal,
                    void (*eval)(const double*, double*, const double*, const int*),
                    const double* deval,  const int* ieval);

void
ss_bisection_overrelax(double* x1,           double* interv,      const double* x0,     const double* z,      
                       const double* w,      const int* a,        const int* doubling,
                       void (*eval)(const double*, double*, const double*, const int*),
                       const double* deval,  const int* ieval);

void
ss_exact_sample(double* x1,           double* interv,       double* g_interv,     const double* x0,
                const double* z,
                void (*eval)(const double*, double*, const double*, const int*),
                const double* deval,  const int* ieval);

void
ss_exact_overrelax(double* x1,           double* interv,  const double* x0, const double* z,
                   void (*eval)(const double*, double*, const double*, const int*),
                   const double* deval,  const int* ieval);

}  /*** end of the namespace Slice_sampler ***/

#ifdef __cplusplus
}
#endif

#endif
