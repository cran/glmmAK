/*** RandomPoiss.h ***/

#ifndef _RANDOM_POISSON_H_
#define _RANDOM_POISSON_H_

#include <cmath>

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_Optim.h"

#include "MatrixRect.h"
#include "MatrixLT.h"
#include "Gspline2.h"
#include "BiGspline2.h"
#include "mcmc_Random.h"
#include "mcmc_common.h"
#include "ll_poisson.h"

class RandomPoiss
{
  private:
  int _p;                // number of random effects
  int _nTheta;           // number of random effects (= _p)
  int _N;                // number of clusters  
  MatrixRect<int> _ni;   // number of observations within each clusters                              /* 1 x _N     */
  int _max_ni;           // maximal number of observations within each cluster
  int _n;                // number of observations                                                   /* = sum(_ni) */

  int _REdist;    // distribution of random effects (_Normal = 0, _Gspline = 1)

  /*** Main values ***/
  MatrixRect<double> _Theta;          // current values of random effects                            /* _nTheta x _N */

  /*** Parameters of the distribution of random effects ***/
  MatrixRect<double> _REMean;         // t(a) = vector of means of random effects                               /* 1 x _nTheta       */
  MatrixLT<double> _REInvVar;         // D^{-1} = inverse variance matrix of random effects                     /* LT(_nTheta)       */
  MatrixLT<double> _REInvVarL;        // L = Cholesky decomposition of D^{-1}, i.e., D^{-1} = L * t(L)          /* LT(_nTheta)       */

  /*** Hyperparameters of the distribution of random effects ***/
  int _prior_for_REMean;                 // prior for RE mean (_Fixed_ or _Normal_)                                        /* 1           */
  MatrixRect<double> _REMeanPriorMean;   // vector of prior means for means of random effects                              /* 1 x _nTheta */
  MatrixRect<double> _REMeanPriorInvVar; // vector of prior inverse variances of means of random effects                   /* 1 x _nTheta */
  MatrixRect<double> _REMeanPriorWMean;  // vector of prior inverse variances of means of random effects weighted by means /* 1 x _nTheta */

  int _prior_for_REInvVar;                  // prior for RE inverse variance (_Fixed, _Wishart or _SDUnif)          /* 1                 */
  MatrixRect<double> _REInvVarPriorDF;      // _prior_for_REInvVar = 
                                            //    _Fixed, _SDUnif: arbitrary
                                            //    _Wishart:        prior degrees of freedom of the Wishart prior  [1 x 1]
                                            //    _GammaIndep:     prior shape parameters of the gamma priors     [1 x_dim]
  MatrixLT<double> _REInvVarPriorInvScale;  // _prior_for_REInvVar 
                                            //    _Fixed: arbitary
                                            //    _Wishart:        prior inverse scale matrix of the Wishart prior           [LT(_nTheta)]
                                            //    _SDUnif:         1/S[i]^2 as the first _dim elements of _REInvVarPriorInvScale._a,
                                            //                     where S[i] is the upper limit of the uniform distribution
                                            //    _GammaIndep:     inverse scale parameters (rates) of the gamma priors as the first _dim
                                            //                     elements of _REInvVarPriorInvScale._a

  /*** Working spaces (updated only when needed) ***/
    /* _nTheta x _N */
  MatrixRect<double> _Theta_REMean;   // random effects minus their mean

    /* 1 x _nTheta */
  MatrixRect<double> _ThetaBar;       // average value (over clusters) of random effects
  MatrixRect<double> _REWMean;        // t(D^{-1}*a)
  MatrixRect<double> _REStdDev;       // standard deviations of random effects (square roots of diag(D))
  MatrixRect<double> _PropTheta;      // space to store proposed value of the cluster-specific random effect etc.
  MatrixRect<double> _PropMean;       // space to store (partial) mean of the proposal etc. 
  MatrixRect<double> _U;              // score of the log-likelihood etc.

    /* _nTheta x _nTheta */
  MatrixLT<double> _REVar;            // D = inversion of _REInvVar                              /* _nTheta x _nTheta */
  MatrixLT<double> _I;                // Hessian of the log-likelihood, inverse variance of the proposal and its decomposition etc.

  /*** Linear predictors and related working spaces ***/
  MatrixRect<double> _eta;        // linear predictor                                                    /* 1 x _n        */
  MatrixRect<double> _Propeta;    // proposed linear predictor                                           /* 1 x _max_ni   */

  /*** Design matrices ***/
  MatrixRect<double> _x;          // design matrix _x[_p, _n] of covariates                              /* _p x _n */
  MatrixLT<double> *_xx;          // array of length _n with matrices _x[*,i]) * t(_x[*,i])              /* _n, _p x _p */

  /*** Additional working spaces ***/
  MatrixRect<double> _work;             // working array needed by Fisher scoring          /* 1 x LT(_p) */
  MatrixRect<double> _work_rwishart;    // working array needed by rwishart                /* 1 x 2*_nTheta*_nTheta                                              */
                                                                                           /* originally (till 12/01/2007): 1 x 2*LT(_nTheta) + _nTheta*_nTheta  */
  MatrixRect<double> _work_invVarSlice; // working array needed by update_InvVarRE2        /* 1 x 4*_nTheta */


  public:

  /***** Constructors and destructors *****/
  RandomPoiss();
  RandomPoiss(const int &p,           const int &N,             const int *ni,
              const double *theta,    const double *mean_ivar,
              const int &REdist,      const int *mean_ivarPrior,  const double *meanPrior,  const double *ivarPrior, 
              const double *x);
  RandomPoiss(const RandomPoiss &P);
  RandomPoiss& operator=(const RandomPoiss &P);
  ~RandomPoiss();

  /*****  Access functions to various components *****/
  inline int
  p() const { return _p;}

  inline int
  nTheta() const { return _nTheta;}

  inline int
  max_ni() const { return _max_ni;}

  inline int
  n() const { return _n;}

  inline int
  N() const { return _N;}

  inline int
  N_nTheta() const { return _N * _nTheta;}

  inline int
  lInvVar() const { return _REInvVar.length();}

  inline int
  prior_for_REMean() const { return _prior_for_REMean;}

  inline int
  prior_for_REInvVar() const { return _prior_for_REInvVar;}


  inline const double*
  ThetaAconst() const { return _Theta.aconst();}

  inline MatrixRect<double>*
  Theta() { return &_Theta;}

  inline const MatrixRect<double>*
  Thetaconst() const { return &_Theta;}


  inline const double*
  REMeanAconst() const { return _REMean.aconst();}

  inline MatrixRect<double>*
  REMean() { return &_REMean;}

  inline const MatrixRect<double>*
  REMeanconst() const { return &_REMean;}


  inline const MatrixRect<double>*
  REMeanPriorWMeanconst() const { return &_REMeanPriorWMean;}


  inline const MatrixRect<double>*
  REMeanPriorInvVarconst() const { return &_REMeanPriorInvVar;}


  inline const MatrixRect<double>*
  REInvVarPriorDFconst() const { return &_REInvVarPriorDF;}


  inline const MatrixLT<double>*
  REInvVarPriorInvScaleconst() const { return &_REInvVarPriorInvScale;}


  inline const double*
  REInvVarAconst() const { return _REInvVar.aconst();}

  inline MatrixLT<double>*
  REInvVar() { return &_REInvVar;}

  inline const MatrixLT<double>*
  REInvVarconst() const { return &_REInvVar;}


  inline const double*
  REInvVarLAconst() const { return _REInvVarL.aconst();}

  inline MatrixLT<double>*
  REInvVarL() { return &_REInvVarL;}


  inline const double*
  REVarAconst() const { return _REVar.aconst();}


  inline MatrixLT<double>*
  REVar() { return &_REVar;}


  inline const double*
  REStdDevAconst() const { return _REStdDev.aconst();}


  inline MatrixRect<double>*
  ThetaBar() { return &_ThetaBar;}


  inline MatrixRect<double>*
  Theta_REMean() { return &_Theta_REMean;}


  inline MatrixRect<double>*
  PropMean() { return &_PropMean;}


  inline MatrixRect<double>*
  U() { return &_U;}

  inline MatrixLT<double>*
  I() { return &_I;}

  inline MatrixRect<double>*
  work_rwishart() { return &_work_rwishart;}

  inline MatrixRect<double>*
  work_invVarSlice() { return &_work_invVarSlice;}


  inline const MatrixRect<double>*
  etaconst() const { return &_eta;}

  inline const double*
  etaAconst() const { return _eta.aconst();}

  inline MatrixRect<double>*
  eta() { return &_eta;}

  inline double*
  etaA() { return _eta.a();}


  inline const double*
  xAconst() const { return _x.aconst();}

  inline const MatrixRect<double>*
  xconst() const { return &_x;}


  inline const MatrixLT<double>*
  xxconst() const { return _xx;}



  /***** Model related functions *****/
  void
  ll_poisson(double *ll,            double *mu,
	     const double *offset,  const int *y,  const double *log_y_factor,  const int &order);


  /***** MCMC related functions for cumulative logit model with NORMAL random effects *****/
  void
  updateRE1(int *accept,           double *ll, 
            double *mu,            double *Propmu,
            const double *offset,  const int *y,   const double *log_y_factor);


  /***** MCMC related functions for cumulative logit model with UNIVARIATE G-SPLINE random effects *****/
  void
  updateRE2(int *accept,           double *ll, 
            double *mu,            double *Propmu,
            const Gspline2 *Gspl,
            const double *offset,  const int *y,    const double *log_y_factor,
            const MatrixRect<int> *alloc);


  /***** MCMC related functions for cumulative logit model with BIVARIATE G-SPLINE random effects *****/
  void
  updateRE2Bi(int *accept,             double *ll, 
              double *mu,              double *Propmu,
              const BiGspline2 *Gspl,
              const double *offset,    const int *y,    const double *log_y_factor,
              const MatrixRect<int> *alloc);


  /***** Utilities *****/
  int
  invert_REInvVar();

  void
  invert_REInvVarFromL();

  void
  print(const int &colmax) const;


};    /*** end of class RandomPoiss ***/


namespace RandomPoissA{

/***** Related function *****/
void
full_Gspline_InvVar_logdens0(const double* x,  double* yu,  const double* pars,  const int* ipars);

void
full_Gspline_InvVar_logdens3(const double* x,     double* yu,        double* ypu,     double* yppu,  
                             const double* pars,  const int* ipars,  const int& what);

}    /*** end of the namespace RandomPoissA ***/

#endif


