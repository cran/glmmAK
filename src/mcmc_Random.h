/*** mcmc_Random.h ***/

#ifndef _MCMC_RANDOM_H_
#define _MCMC_RANDOM_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_Optim.h"

#include "MatrixRect.h"
#include "MatrixLT.h"

#include "Gspline2.h"
#include "BiGspline2.h"

#include "mcmc_common.h"
#include "Slice_sampler2.h"

namespace mcmc_Random {

  enum _Distribution_Random_Effect_ {_None, _Normal, _Gspline};
  enum priorForInvVarTypes {_Fixed, _Wishart, _SDUnif, _GammaIndep};   // types of priors for inverse of the covariance matrix  of random effects
  enum priorForMeanTypes {_Fixed_, _Normal_};                          // types of priors for means of random effects

  const double _zero_invvariance = 1e-5;        // lower limit for inverse variance parameters 
                                                // *  used by updateInvVarRE2 function

  /** constants for Newton-Raphson solver within the slice sampler **/
  const int _maxiter_solver_nr = 10;           // maximum number of NR iterations
  const double _toler_solver_nr = 1e-3;        // tolerance to detect convergence of NR


/***** =============================================== *****/
/*****                                                 *****/
/***** MODELS WITH NORMALLY DISTRIBUTED RANDOM EFFECTS *****/
/*****                                                 *****/
/***** =============================================== *****/
  void
  updateMeanRE1(MatrixRect<double> *REMean,                  MatrixRect<double> *ThetaBar,      
                MatrixRect<double> *PropMean,                MatrixRect<double> *U,                         MatrixLT<double> *I,
                const MatrixRect<double> *Theta,             const MatrixLT<double> *REInvVar,
                const MatrixRect<double> *REMeanPriorWMean,  const MatrixRect<double> *REMeanPriorInvVar,   const int &prior_for_REMean,
                const int &nTheta,                           const int &N);

  void
  updateInvVarRE1(MatrixLT<double> *REInvVar,                 MatrixLT<double> *REInvVarL,
                  MatrixRect<double> *Theta_REMean,           MatrixRect<double> *U,        MatrixLT<double> *I,  MatrixRect<double> *work_rwishart,
                  const MatrixRect<double> *Theta,            const MatrixRect<double> *REMean,
                  const MatrixRect<double> *REInvVarPriorDF,  const MatrixLT<double> *REInvVarPriorInvScale,      const int &prior_for_REInvVar,
                  const int &nTheta,                          const int &N);


/***** =============================================== *****/
/*****                                                 *****/
/***** MODELS WITH G-SPLINE DISTRIBUTED RANDOM EFFECTS *****/
/*****                                                 *****/
/***** =============================================== *****/

  /*** Update of G-spline intercept parameters (alpha) ***/
  /*** ----------------------------------------------- ***/
  void
  updateMeanRE2(MatrixRect<double> *REMean,                  MatrixRect<double> *ThetaBar,
                const MatrixRect<double> *Theta,
                const MatrixRect<double> *REMeanPriorWMean,  const MatrixRect<double> *REMeanPriorInvVar,   const int &prior_for_REMean,
                const Gspline2 *Gspl,                        const MatrixRect<int> *alloc,
                const int &nTheta,                           const int &N);

  void
  updateMeanRE2Bi(MatrixRect<double> *REMean,                  MatrixRect<double> *ThetaBar,
                  MatrixRect<double> *PropMean,                MatrixRect<double> *U,                         MatrixLT<double> *I,
                  const MatrixRect<double> *Theta,
                  const MatrixRect<double> *REMeanPriorWMean,  const MatrixRect<double> *REMeanPriorInvVar,   const int &prior_for_REMean,
                  const BiGspline2 *Gspl,                      const MatrixRect<int> *alloc,
                  const int &nTheta,                           const int &N);



  /*** Update of G-spline scale parameters (tau) ***/
  /*** ----------------------------------------- ***/
  void
  updateInvVarRE2(MatrixLT<double> *REInvVar,                 MatrixLT<double> *REInvVarL,                    MatrixLT<double> *REVar,
                  Gspline2 *Gspl,                             MatrixRect<double> *work_invVarSlice,           
                  const MatrixRect<double> *Theta,            const MatrixRect<double> *REMean,
                  const MatrixRect<double> *REInvVarPriorDF,  const MatrixLT<double> *REInvVarPriorInvScale,  const int &prior_for_REInvVar,
                  const MatrixRect<int> *alloc,
                  const int &nTheta,                          const int &N,                                   const int &overrelax);

  void
  updateInvVarRE2Bi(MatrixLT<double> *REInvVar,                 MatrixLT<double> *REInvVarL,                    MatrixLT<double> *REVar,
                    BiGspline2 *Gspl,                           MatrixRect<double> *work_invVarSlice,           
                    const MatrixRect<double> *Theta,            const MatrixRect<double> *REMean,
                    const MatrixRect<double> *REInvVarPriorDF,  const MatrixLT<double> *REInvVarPriorInvScale,  const int &prior_for_REInvVar,
                    const MatrixRect<int> *alloc,
                    const int &nTheta,                          const int &N,                                   const int &overrelax);

  void
  full_Gspline_InvVar_pars(double *pars, 
                           const int &BiGspline,             const MatrixRect<int> *alloc,  
                           const MatrixRect<double> *knots,  const int *K,                   const double *invsigma2,
                           const MatrixRect<double> *Theta,  const int &nTheta,              const int &N,          
                           const double *REMean,
                           const int &prior_for_REInvVar,    const double *REInvVarPriorDF,  const double *REInvVarPriorInvScale);

  void
  full_Gspline_InvVar_logdens0(const double* x,  double* yu,  const double* pars,  const int* ipars);

  void
  full_Gspline_InvVar_logdens3(const double* x,     double* yu,        double* ypu,     double* yppu,  
                               const double* pars,  const int* ipars,  const int& what);

}

#endif
