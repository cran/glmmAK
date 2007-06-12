/*** PredictorCL.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  12/09/2006
//     MAJOR MODIFICATIONS:  13/10/2006  (order of 'p' and 'q' parameters in _Theta changed)
//                            
//                  BASICS:  12/09/2006
//            change_Theta:  29/01/2007
//            ll_cumlogit2:  27/09/2006
//                 update1:  14/09/2006
//                   print:  27/09/2006
//
// PURPOSE: Class to hold parameters related to the (linear) predictor in the cumulative logits GALMM
//
//          * it also holds design matrices and current values of the predictor itself  
//          * it implements also MCMC updates of these parameters
//
/* ********************************************************************************* */

#include "PredictorCL.h"

/* ********************************************************************************* */
/* Constructors and destructors                                                      */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */

/***** Nonparametric constructor *****/
PredictorCL::PredictorCL()
  : _p(0), _q(0), _nTheta(0), _n(0), _C(1), _xx(NULL), _vv(NULL), _xv(NULL)
{
}


/***** Parametric constructor 1 *****/
/*                                                                                              */
/* priorp[2*(C*q+p)]: parameters of the prior of theta                                          */
/*               priorp[0,...,C*q+p-1]: prior means                                             */
/*       priorp[c*q+p,...,2*(C*q+p)-1]: prior inverse variances                                 */
/*                                                                                              */
/*                                                                                              */
PredictorCL::PredictorCL(const int &p,         const int &q,         const int &n,           const int &C,
                         const double *theta,  const double *priorp,
                         const double *x,      const double *v)
{
  int i;
  const double *xP, *vP;

  _p = p;
  _q = q;
  _n = n;
  if (C < 1) throw returnR("Error in PredictorCL.cpp: PredictorCL::PredictorCL(...), C < 1", 1);
  _C = C;
  _nTheta = _p + _C*_q;
  _work = MatrixRect<double>(1, 5*_C);

  if (_nTheta){
    _Theta = MatrixRect<double>(1, _nTheta, theta);
    _PropTheta = MatrixRect<double>(1, _nTheta, theta);
    _PriorMean = MatrixRect<double>(1, _nTheta, priorp);
    _PriorInvVar = MatrixRect<double>(1, _nTheta, priorp+_nTheta);
    _PriorWMean = MatrixRect<double>(1, _nTheta, priorp+_nTheta, priorp, 1);

    _PropMean = MatrixRect<double>(1, _nTheta);
    _U = MatrixRect<double>(1, _nTheta);
    _I = MatrixLT<double>(_nTheta);

    if (_n){
      _x = MatrixRect<double>(_p, _n, x);
      _v = MatrixRect<double>(_q, _n, v);

      _etaX = MatrixRect<double>(1, _n);
      _PropetaX = MatrixRect<double>(1, _n);
      if (_p){
        _etaX.bA(&_Theta, &_x, _C*_q);
        _PropetaX = _etaX;
      }

      _etaV = MatrixRect<double>(_C, _n);
      _PropetaV = MatrixRect<double>(_C, _n);
      if (_q){
        _etaV.bA2(&_Theta, &_v, 0);
        _PropetaV = _etaV;
      }

      _eta = MatrixRect<double>(_C, _n);
      _Propeta = MatrixRect<double>(_C, _n);
      _eta.b_plus_rowsA(_etaX.aconst(), &_etaV);
      _Propeta = _eta;

      _xx = new MatrixLT<double>[_n];
      _vv = new MatrixLT<double>[_n];
      _xv = new MatrixRect<double>[_n];
      if (!_xx || !_vv || !_xv) throw returnR("Out of memory in PredictorCL.cpp: PredictorCL::PredictorCL(...)", 99);

      xP = x;
      vP = v;           
      for (i = 0; i < _n; i++){
        _xx[i] = MatrixLT<double>(_p, xP);
        _vv[i] = MatrixLT<double>(_q, vP);
        _xv[i] = MatrixRect<double>(_p, _q, xP, vP);

        xP += _p;
        vP += _q;
      }
    }
  }
  else{    /*** else from if (_nTheta) ***/
    if (_n) _eta = MatrixRect<double>(_C, _n);
  }

  if (!_nTheta || !_n){
    _xx = NULL;
    _vv = NULL;
    _xv = NULL;
  }
}


/***** Copy constructor *****/
PredictorCL::PredictorCL(const PredictorCL &P)
{
  int i;

  _p = P._p;
  _q = P._q;
  _n = P._n;
  _C = P._C;
  _nTheta = P._nTheta;
  _work = P._work;

  _Theta = P._Theta;
  _PropTheta = P._PropTheta;
  _PriorMean = P._PriorMean;
  _PriorInvVar = P._PriorInvVar;
  _PriorWMean = P._PriorWMean;

  _PropMean = P._PropMean;
  _U = P._U;
  _I = P._I;

  _x = P._x;
  _v = P._v;

  _etaX = P._etaX;
  _etaV = P._etaV;
  _eta = P._eta;
  _PropetaX = P._PropetaX;
  _PropetaV = P._PropetaV;
  _Propeta = P._Propeta;

  if (_nTheta && _n){
    _xx = new MatrixLT<double>[_n];
    _vv = new MatrixLT<double>[_n];
    _xv = new MatrixRect<double>[_n];
    if (!_xx || !_vv || !_xv) throw returnR("Out of memory in PredictorCL.cpp: PredictorCL::PredictorCL(P)", 99);
    for (i = 0; i < _n; i++){
      _xx[i] = P._xx[i];
      _vv[i] = P._vv[i];
      _xv[i] = P._xv[i];
    }
  }
  else{
    _xx = NULL;
    _vv = NULL;
    _xv = NULL;
  }
}


/***** Assignment operator *****/
PredictorCL&
PredictorCL::operator=(const PredictorCL &P)
{
  int i;

  if (_n && _nTheta){
    delete [] _xx;
    delete [] _vv;
    delete [] _xv;
  }

  _p = P._p;
  _q = P._q;
  _n = P._n;
  _C = P._C;
  _nTheta = P._nTheta;
  _work = P._work;

  _Theta = P._Theta;
  _PropTheta = P._PropTheta;
  _PriorMean = P._PriorMean;
  _PriorInvVar = P._PriorInvVar;
  _PriorWMean = P._PriorWMean;

  _PropMean = P._PropMean;
  _U = P._U;
  _I = P._I;

  _x = P._x;
  _v = P._v;

  _etaX = P._etaX;
  _etaV = P._etaV;
  _eta = P._eta;
  _PropetaX = P._PropetaX;
  _PropetaV = P._PropetaV;
  _Propeta = P._Propeta;

  if (_nTheta && _n){
    _xx = new MatrixLT<double>[_n];
    _vv = new MatrixLT<double>[_n];
    _xv = new MatrixRect<double>[_n];
    if (!_xx || !_vv || !_xv) throw returnR("Out of memory in PredictorCL.cpp: PredictorCL::operator=", 99);
    for (i = 0; i < _n; i++){
      _xx[i] = P._xx[i];
      _vv[i] = P._vv[i];
      _xv[i] = P._xv[i];
    }
  }
  else{
    _xx = NULL;
    _vv = NULL;
    _xv = NULL;
  }

  return *this;
}


/***** Destructor *****/
PredictorCL::~PredictorCL()
{
  if (_n && _nTheta){
    delete [] _xx;
    delete [] _vv;
    delete [] _xv;
  }
}


/* ************************************************************************************************************************ */
/* Useful utilities                                                                                                         */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* change_Theta:  Change _Theta and update _etaX, _etaV, _eta                                                               */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
//
//  newTheta[_nTheta]:   New value of theta
//
void
PredictorCL::change_Theta(const double *newTheta)
{
  static int i;
  static double *dP;
  static const double *cdP;  

  if (_nTheta > 0){

    /*** Copy newTheta to _Theta ***/
    dP  = _Theta.a();
    cdP = newTheta;
    for (i = 0; i < _nTheta; i++){
      *dP = *cdP;
      dP++;
      cdP++;
    }

    /*** Update _etaX, _etaV, _eta  ***/
    if (_n){
      if (_q) _etaV.bA2(&_Theta, &_v, 0);
      if (_p) _etaX.bA(&_Theta, &_x, _C*_q);
      _eta.b_plus_rowsA(_etaX.aconst(), &_etaV);     
    }
  }

  return;
}


/* ************************************************************************************************************************ */
/* Model related functions                                                                                                  */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* ll_cumlogit2:  Apply either ll_cumlogitFS2 or ll_cumlogitNR2 function to the objects of PredictorCL                      */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
void
PredictorCL::ll_cumlogit2(double *ll,           double *prob,   int *anyZero,  
			  const double *offset,  const int *y,   const int &update_etas,  const int &order,
                          void (*ll_2)(double*, MatrixRect<double>*, MatrixLT<double>*, 
                                       double*, double*, double*, double*, double*, int*,
                                       const double*, const double*, const int*, const double*, const double*,
                                       const MatrixLT<double>*, const MatrixLT<double>*, const MatrixRect<double>*,
                                       const int&, const int&, const int&, const int&, const int&, const int&)
                          )
{
  ll_2(ll, &_U, &_I, 
       _etaX.a(),  _etaV.a(), _eta.a(), prob, _work.a(), anyZero,
       offset, _Theta.aconst(), y, _x.aconst(), _v.aconst(), 
       _xx, _vv, _xv,
       _n, _C, _p, _q, update_etas, order);
  return;
}


/* ************************************************************************************************************************ */
/* MCMC related functions                                                                                                   */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* update1: Update of regression coefficients                                                                               */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/*     accept:  OUTPUT: indicator of acceptance (1 if accepted, 0 if not accepted)                                          */
/*         ll:   INPUT: current value of the log-likelihood                                                                 */
/*              OUTPUT: new value of the log-likelihood                                                                     */
/*       prob:   INPUT: current values of category probabilities                                                            */
/*              OUTPUT: new values of category probabilities                                                                */
/*              * matrix (C+1) x n                                                                                          */
/*   Propprob:  working space                                                                                               */
/*              * matrix (C+1) x n                                                                                          */
/*     offset:  offset vector                                                                                               */
/*              * vector of length n                                                                                        */
/*          y:  vector of observations                                                                                      */
/*              * vector of length n                                                                                        */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
void
PredictorCL::update1(int *accept,                       double *ll, 
                     MatrixRect<double> *prob,          MatrixRect<double> *Propprob,
                     const MatrixRect<double> *offset,  const MatrixRect<int> *y)
{
  mcmc_common::update_reg_gamermanCL(accept, ll, &_U, &_I, 
                  _etaX.a(), _PropetaX.a(), _etaV.a(), _PropetaV.a(), _eta.a(), _Propeta.a(), prob->a(), Propprob->a(), _work.a(),
                  offset->aconst(), _Theta.a(), _PropTheta.a(), y->aconst(), _x.aconst(), _v.aconst(), _xx, _vv, _xv,
                  _n, _C, _p, _q, _PriorMean.aconst(), _PriorWMean.aconst(), _PriorInvVar.aconst(), NULL, true, &_PropMean, 
                  ll_cumlogitFS2, ll_cumlogitNR2, "PredictorCL::update1");

  return;
}


/* ************************************************************************************************************************ */
/* Utilities                                                                                                                */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
void
PredictorCL::print(const int &colmax) const
{
  int cmax = (colmax < 0 ? _n : (colmax <= _n ? colmax : _n));

  Rprintf("\nObject of class PredictorCL:\n");
  Rprintf("============================\n");
  Rprintf("q=%d,  C=%d,  p=%d,  nTheta=%d\n", _q, _C, _p, _nTheta);
  Rprintf("n=%d\n", _n);

  Rprintf("\nTHETA:\n");
  _Theta.print(0);
  
  Rprintf("\nPRIOR FOR THETA:\n");
  Rprintf("* PriorMean: ");
  _PriorMean.print(0);
  Rprintf("* PriorInvVar: ");
  _PriorInvVar.print(0);
  Rprintf("* PriorWMean: ");
  _PriorWMean.print(0);

  Rprintf("\nLINEAR PREDICTORS:\n");
  Rprintf("* eta:\n");
  _eta.print(0, cmax);
  Rprintf("* etaV:\n");
  _etaV.print(0, cmax);
  Rprintf("* etaX: ");
  _etaX.print(0, cmax);

  Rprintf("\nDESIGN MATRICES:\n");
  Rprintf("* V:\n");
  _v.print(0, cmax);
  Rprintf("* X:\n");
  _x.print(0, cmax);

  Rprintf("\n");
  return;
}
