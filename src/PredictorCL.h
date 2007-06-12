/*** PredictorCL.h ***/

#ifndef _PREDICTOR_CL_H_
#define _PREDICTOR_CL_H_

#include <cmath>

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

#include "MatrixRect.h"
#include "MatrixLT.h"
#include "mcmc_common.h"
#include "ll_cumlogit.h"

class PredictorCL
{
  private:
  int _p;         // number of covariates proportional w.r.t. odds
  int _q;         // number of covariates not proportional w.r.t. odds
  int _nTheta;    // number of regression parameters (= _p + _C*_q)
  int _n;         // number of observations
  int _C;         // number of logits (= number of response categories minus 1)
   
  MatrixRect<double> _Theta;          // current values of regression parameters                             /* 1 x _nTheta  */
  MatrixRect<double> _PropTheta;      // proposed value of the regression parameter                          /* 1 x _nTheta  */

  MatrixRect<double> _PriorMean;      // vector of prior means                                               /* 1 x _nTheta  */
  MatrixRect<double> _PriorInvVar;    // vector of prior inverse variances                                   /* 1 x _nTheta  */
  MatrixRect<double> _PriorWMean;     // vector of R^{-1}*a, where R = diag(PriorVar), a = PriorMean         /* 1 x _nTheta  */
  
  MatrixRect<double> _PropMean;   // space to store (partial) mean of the proposal                                          /* 1 x _nTheta */  
  MatrixRect<double> _U;          // score of the log-likelihood, used also as a working space                              /* 1 x _nTheta */
  MatrixLT<double> _I;            // Hessian of the log-likelihood, inverse variance of the proposal and its decomposition  /* _nTheta x _nTheta */  

  MatrixRect<double> _etaX;       // parts of the predictor proportional w.r.t. odds                     /* 1 x _n      */
  MatrixRect<double> _etaV;       // parts of the predictor not proportional w.r.t. odds                 /* _C x _n     */
  MatrixRect<double> _eta;        // _etaX + _etaV                                                       /* _C x _n     */
  MatrixRect<double> _PropetaX;   // proposed parts of the predictor proportional w.r.t. odds            /* 1 x _n      */
  MatrixRect<double> _PropetaV;   // proposed parts of the predictor not proportional w.r.t. odds        /* _C x _n     */
  MatrixRect<double> _Propeta;    // proposed _etaX + _etaV                                              /* _C x _n     */

  MatrixRect<double> _x;          // design matrix _x[_p, _n] of covariates proportional w.r.t. odds     /* _p x _n */
  MatrixRect<double> _v;          // design matrix _v[_q, _n] of covariates proportional w.r.t. odds     /* _q x _n */
  MatrixLT<double> *_xx;          // array of length _n with matrices _x[*,i]) * t(_x[*,i])              /* _n, _p x _p */
  MatrixLT<double> *_vv;          // array of length _n with matrices _v[*,i]) * t(_v[*,i])              /* _n, _q x _q */
  MatrixRect<double> *_xv;        // array of length _n with matrices _x[*,i]) * t(_v[*,i])              /* _n, _p x _q */

  MatrixRect<double> _work;       // working array needed by Fisher scoring                              /* 1 x 5*_C            */

  public:

  /***** Constructors and destructors *****/
  PredictorCL();
  PredictorCL(const int &p,         const int &q,          const int &n,           const int &C,
              const double *theta,  const double *priorp,
              const double *x,      const double *v);
  PredictorCL(const PredictorCL &P);
  PredictorCL& operator=(const PredictorCL &P);
  ~PredictorCL();

  /*****  Access functions to various components *****/
  inline int
  p() const { return _p;}

  inline int
  q() const { return _q;}

  inline int
  nTheta() const { return _nTheta;}

  inline int
  n() const { return _n;}

  inline int
  C() const { return _C;}


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
  etaXconst() const { return &_etaX;}

  inline const double*
  etaXAconst() const { return _etaX.aconst();}

  inline MatrixRect<double>*
  etaX() { return &_etaX;}

  inline double*
  etaXA() { return _etaX.a();}


  inline const MatrixRect<double>*
  etaVconst() const { return &_etaV;}

  inline const double*
  etaVAconst() const { return _etaV.aconst();}

  inline MatrixRect<double>*
  etaV() { return &_etaV;}

  inline double*
  etaVA() { return _etaV.a();}


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


  inline const double*
  vAconst() const { return _v.aconst();}

  inline const MatrixRect<double>*
  vconst() const { return &_v;}


  inline const MatrixLT<double>*
  xxconst() const { return _xx;}

  inline const MatrixLT<double>*
  vvconst() const { return _vv;}

  inline const MatrixRect<double>*
  xvconst() const { return _xv;}


  inline double*
  workA() { return _work.a();}


  /***** Useful utilities *****/
  void
  change_Theta(const double *newTheta);


  /***** Model related functions *****/
  void
  ll_cumlogit2(double *ll,           double *prob,   int *anyZero,  
	       const double *offset,  const int *y,   const int &update_etas,  const int &order,
               void (*ll_2)(double*, MatrixRect<double>*, MatrixLT<double>*, 
                            double*, double*, double*, double*, double*, int*,
                            const double*, const double*, const int*, const double*, const double*,
                            const MatrixLT<double>*, const MatrixLT<double>*, const MatrixRect<double>*,
                            const int&, const int&, const int&, const int&, const int&, const int&)
              );


  /***** MCMC related functions *****/
  void
  update1(int *accept,                       double *ll, MatrixRect<double> *prob,    MatrixRect<double> *Propprob,
          const MatrixRect<double> *offset,  const MatrixRect<int> *y);

  
  /***** Utilities *****/
  void
  print(const int &colmax) const;


};  /* end of class PredictorCL */

#endif

