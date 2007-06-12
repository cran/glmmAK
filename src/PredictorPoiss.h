/*** PredictorPoiss.h ***/

#ifndef _PREDICTOR_POISS_H_
#define _PREDICTOR_POISS_H_

#include <cmath>

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

#include "MatrixRect.h"
#include "MatrixLT.h"
#include "mcmc_common.h"
#include "ll_poisson.h"

class PredictorPoiss
{
  private:
  int _p;         // number of covariates
  int _nTheta;    // number of regression parameters (= _p)
  int _n;         // number of observations
   
  MatrixRect<double> _Theta;          // current values of regression parameters                             /* 1 x _nTheta  */
  MatrixRect<double> _PropTheta;      // proposed value of the regression parameter                          /* 1 x _nTheta  */

  MatrixRect<double> _PriorMean;      // vector of prior means                                               /* 1 x _nTheta  */
  MatrixRect<double> _PriorInvVar;    // vector of prior inverse variances                                   /* 1 x _nTheta  */
  MatrixRect<double> _PriorWMean;     // vector of R^{-1}*a, where R = diag(PriorVar), a = PriorMean         /* 1 x _nTheta  */
  
  MatrixRect<double> _PropMean;   // space to store (partial) mean of the proposal                                          /* 1 x _nTheta */  
  MatrixRect<double> _U;          // score of the log-likelihood, used also as a working space                              /* 1 x _nTheta */
  MatrixLT<double> _I;            // Hessian of the log-likelihood, inverse variance of the proposal and its decomposition  /* _nTheta x _nTheta */  

  MatrixRect<double> _eta;        // linear predictor                                                    /* 1 x _n     */
  MatrixRect<double> _Propeta;    // proposed linear predictor                                           /* 1 x _n     */

  MatrixRect<double> _x;          // design matrix _x[_p, _n] of covariates proportional w.r.t. odds     /* _p x _n */
  MatrixLT<double> *_xx;          // array of length _n with matrices _x[*,i]) * t(_x[*,i])              /* _n, _p x _p */

  MatrixRect<double> _work;       // working array                                                       /* 1 x LT(_p)  */

  public:

  /***** Constructors and destructors *****/
  PredictorPoiss();
  PredictorPoiss(const int &p,  const int &n,  const double *theta,  const double *priorp,  const double *x);
  PredictorPoiss(const PredictorPoiss &P);
  PredictorPoiss& operator=(const PredictorPoiss &P);
  ~PredictorPoiss();

  /*****  Access functions to various components *****/
  inline int
  p() const { return _p;}

  inline int
  nTheta() const { return _nTheta;}

  inline int
  n() const { return _n;}


  inline const double*
  ThetaAconst() const { return _Theta.aconst();}

  inline MatrixRect<double>*
  Theta() { return &_Theta;}

  inline const MatrixRect<double>*
  Thetaconst() const { return &_Theta;}


  inline MatrixRect<double>*
  U() { return &_U;}

  inline MatrixLT<double>*
  I() { return &_I;}


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
  xconst() const { return &_x; }


  inline const MatrixLT<double>*
  xxconst() const { return _xx;}


  inline double*
  workA() { return _work.a();}


  /***** Useful utilities *****/
  void
  change_Theta(const double *newTheta);


  /***** Model related functions *****/
  void
  ll_poisson(double *ll,            double *mu,
	     const double *offset,  const int *y,  const double *log_y_factor,  const int &order);


  /***** MCMC related functions *****/
  void
  update1(int *accept,           double *ll, 
          double *mu,            double *Propmu,
          const double *offset,  const int *y,     const double *log_y_factor);

  
  /***** Utilities *****/
  void
  print(const int &colmax) const;


};  /* end of class PredictorCL */

#endif

