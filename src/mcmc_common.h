/*** mcmc_common.h ***/

#ifndef _MCMC_COMMON_H_
#define _MCMC_COMMON_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

#include "AK_BLAS_LAPACK.h"

#include "MatrixLT.h"
#include "MatrixRect.h"
#include "MatrixUtil.h"
#include "Mvtdist2.h"
#include "Mvtdist3.h"
#include "Random2.h"
#include "util_cumlogit.h"
#include "ll_poisson.h"

namespace mcmc_common{

const int _AK_MAX_PD_ATTEMPTS_2 = 10;      /*** Maximal number of attempts to make the minus Hessian positive definite ***/
const double _AK_EPS_PD_ATTEMPT_2 = 0.01;  /*** Epsilon for 'chol_solvePD' function                                    ***/

void
update_reg_gamermanCL(int *accept,                       
                      double *ll,
                      MatrixRect<double> *U,
                      MatrixLT<double> *I,
                      double *etaX,                      double *PropetaX,
                      double *etaV,                      double *PropetaV,
                      double *eta,                       double *Propeta,
                      double *prob,                      double *Propprob,
                      double *work,
                      const double *offset,              
                      double *Theta,                     double *PropTheta,
                      const int *y,
                      const double *x,
                      const double *v,
                      const MatrixLT<double> *xx,
                      const MatrixLT<double> *vv,
                      const MatrixRect<double> *xv,
                      const int &nObs,
                      const int &C,
                      const int &p,
                      const int &q,
                      const double *PriorMean,
                      const double *PriorWMean,
                      const double *PriorInvVar,
                      const MatrixLT<double> *PriorInvVarL,
                      const bool &diagInvVar,
                      MatrixRect<double> *PropMean,
                      void (*ll_FS2)(double*, MatrixRect<double>*, MatrixLT<double>*, 
                                     double*, double*, double*, double*, double*, int*,
                                     const double*, const double*, const int*, const double*, const double*,
                                     const MatrixLT<double>*, const MatrixLT<double>*, const MatrixRect<double>*,
                                     const int&, const int&, const int&, const int&, const int&, const int&),
                      void (*ll_NR2)(double*, MatrixRect<double>*, MatrixLT<double>*, 
                                     double*, double*, double*, double*, double*, int*,
                                     const double*, const double*, const int*, const double*, const double*,
                                     const MatrixLT<double>*, const MatrixLT<double>*, const MatrixRect<double>*,
                                     const int&, const int&, const int&, const int&, const int&, const int&),
                      const char *caller);

void
update_reg_gamermanPoiss(int *accept,                       
                         double *ll,
                         double *U,
                         double *I,
                         double *eta,                       double *Propeta,
                         double *mu,                        double *Propmu,
                         double *work,
                         const double *offset,              
                         double *Theta,                     double *PropTheta,
                         const int *y,
                         const double *log_y_factor,
                         const double *x,
                         const MatrixLT<double> *xx,
                         const int &nObs,
                         const int &p,
                         const double *PriorMean,
                         const double *PriorWMean,
                         const double *PriorInvVar,
                         const double *PriorInvVarL,
                         const bool &diagInvVar,
                         double *PropMean,
                         const char *caller);

void
update_norm_mean1(double *alpha, 
                  MatrixLT<double> *I,     
                  MatrixRect<double> *U,
                  MatrixRect<double> *PropMean,
                  const double *ObsBar,
                  const double *InvVar,
                  const double *PriorWMean,
                  const double *PriorInvVar,
                  const bool &diagPriorVar,
                  const int &nObs,
                  const int &dim,
                  const char *caller);

void
update_norm_invvar_wishart(double *InvVar,
                           MatrixLT<double> *InvS_Full,
                           double *Obs_alpha,
                           double *work_rwishart,
                           const double *Obs,
                           const double *alpha,
                           const double *df_Wishart,
                           const double *InvS_Wishart,
                           const int &nObs,
                           const int &dim,
                           const char *caller);

void
update_norm_invvar_sduniform(double *InvVar,
                             double *rate,
                             const double *Obs,
                             const double *alpha,
                             const double *invS2,
                             const int &nObs,
                             const int &dim,
                             const char *caller);

}

#endif

