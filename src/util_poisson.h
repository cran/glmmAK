/*** util_poisson.h ***/

#ifndef _UTIL_POISSON_H_
#define _UTIL_POISSON_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

namespace util_poisson{

void
ecount_etas_poisson(double *ecount,   double *eta,   const double *offset,  const double *theta,  const double *x,
		    const int &nObs,  const int &p);

void
ecount_etas_etaREs_poisson(double *ecount,        double *etaF,          double *etaRE,          double *eta, 
                           const double *offset,  const double *theta,   const double *thetaRE,  
                           const double *x,       const double *xRE,     const int &nCluster,    const int *ni,
                           const int &p,          const int &pRE);

}   /*** end of the namespace util_poisson ***/

#endif
