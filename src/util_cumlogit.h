/*** util_cumlogit.h ***/

#ifndef _UTIL_CUMLOGIT_H_
#define _UTIL_CUMLOGIT_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

#include "MatrixLT.h"
#include "MatrixRect.h"

namespace util_cumlogit{

void
prob_cumlogit(MatrixRect<double> *prob,
              const MatrixRect<double> *offset,
              const MatrixRect<double> *etaX,
              const MatrixRect<double> *etaV);

void
prob01_cumlogit(int *anyZero,
                double *prob,
                double *etaX,
                double *etaV,
                double *h,
                double *dh,
                const double *offset,
                const double *theta,
                const double *x,
                const double *v,
                const int &p,
                const int &q,
                const int &C);


/*****************************************************************************************/


void
prob_cumlogit2(int *anyZero,
               double *prob,
               const double *offset,
               const double *etaX,
               const double *etaV,
               const int &nObs,
               const int &C);

void
etas_cumlogit2(double *etaX,
               double *etaV,
               double *eta,
               const double *theta,
               const double *x,
               const double *v,
               const int &nObs,
               const int &C,
               const int &p,
               const int &q);

void
prob_etas_cumlogit2(int *anyZero,
                    double *prob,
                    double *etaX,
                    double *etaV,
                    double *eta,
                    const double *offset,
                    const double *theta,
                    const double *x,
                    const double *v,
                    const int &nObs,
                    const int &C,
                    const int &p,
                    const int &q);

void
prob01_cumlogit2(double *h,
                 double *dh,
                 const double *prob,
                 const int &C);

void
prob_etas_etaREs_cumlogit2(int *anyZero,
                           double *prob,
                           double *etaX,
                           double *etaV,
                           double *etaXRE,
                           double *etaVRE,
                           double *eta,
                           const double *offset,
                           const double *theta,
                           const double *thetaRE,
                           const double *x,
                           const double *v,
                           const double *xRE,
                           const double *vRE,
                           const int &nCluster,
                           const int *ni,
                           const int &C,
                           const int &p,
                           const int &q,
                           const int &pRE,
                           const int &qRE);

}  /*** end of the namespace util_cumlogit ***/

#endif
