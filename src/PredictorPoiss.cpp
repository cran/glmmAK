/*** PredictorPoiss.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  13/02/2007
//                            
//                  BASICS:  13/02/2007
//            change_Theta:  13/02/2007
//              ll_poisson:  13/02/2007
//                 update1:  13/02/2007
//                   print:  13/02/2007
//
// PURPOSE: Class to hold parameters related to the (linear) predictor in the Poisson log-linear model
//
//          * it also holds design matrices and current values of the predictor itself  
//          * it implements also MCMC updates of these parameters
//
/* ********************************************************************************* */

#include "PredictorPoiss.h"

/* ********************************************************************************* */
/* Constructors and destructors                                                      */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */

/***** Nonparametric constructor *****/
PredictorPoiss::PredictorPoiss()
  : _p(0), _nTheta(0), _n(0), _xx(NULL)
{
}


/***** Parametric constructor 1 *****/
/*                                                                                              */
/* priorp[2*p]: parameters of the prior of theta                                                */
/*               priorp[0,...,p-1]: prior means                                                 */
/*             priorp[p,...,2*p-1]: prior inverse variances                                     */
/*                                                                                              */
/*                                                                                              */
PredictorPoiss::PredictorPoiss(const int &p,  const int &n,  const double *theta,  const double *priorp,  const double *x)
{
  int i, LTp;
  const double *xP;

  _p = p;
  _n = n;
  _nTheta = _p;

  LTp = (_p*(_p+1))/2;
  if (LTp) _work = MatrixRect<double>(1, LTp);
  else     _work = MatrixRect<double>(1, 1);  

  if (_nTheta){
    _Theta       = MatrixRect<double>(1, _nTheta, theta);
    _PropTheta   = MatrixRect<double>(1, _nTheta, theta);
    _PriorMean   = MatrixRect<double>(1, _nTheta, priorp);
    _PriorInvVar = MatrixRect<double>(1, _nTheta, priorp+_nTheta);
    _PriorWMean  = MatrixRect<double>(1, _nTheta, priorp+_nTheta, priorp, 1);

    _PropMean = MatrixRect<double>(1, _nTheta);
    _U        = MatrixRect<double>(1, _nTheta);
    _I        = MatrixLT<double>(_nTheta);

    if (_n){
      _x = MatrixRect<double>(_p, _n, x);

      _eta     = MatrixRect<double>(1, _n);
      _Propeta = MatrixRect<double>(1, _n);
      if (_p){
        _eta.bA(&_Theta, &_x, 0);
        _Propeta = _eta;
      }

      _xx = new MatrixLT<double>[_n];
      if (!_xx) throw returnR("Out of memory in PredictorPoiss.cpp: PredictorPoiss::PredictorPoiss(...)", 99);

      xP = x;
      for (i = 0; i < _n; i++){
        _xx[i] = MatrixLT<double>(_p, xP);

        xP += _p;
      }
    }
  }
  else{    /*** else from if (_nTheta) ***/
    if (_n){
      _eta    = MatrixRect<double>(1, _n);
    }
  }

  if (!_nTheta || !_n){
    _xx = NULL;
  }
}


/***** Copy constructor *****/
PredictorPoiss::PredictorPoiss(const PredictorPoiss &P)
{
  int i;

  _p = P._p;
  _n = P._n;
  _nTheta = P._nTheta;
  _work = P._work;

  _Theta       = P._Theta;
  _PropTheta   = P._PropTheta;
  _PriorMean   = P._PriorMean;
  _PriorInvVar = P._PriorInvVar;
  _PriorWMean  = P._PriorWMean;

  _PropMean = P._PropMean;
  _U = P._U;
  _I = P._I;

  _x = P._x;

  _eta     = P._eta;
  _Propeta = P._Propeta;

  if (_nTheta && _n){
    _xx = new MatrixLT<double>[_n];
    if (!_xx) throw returnR("Out of memory in PredictorPoiss.cpp: PredictorPoiss::PredictorPoiss(P)", 99);
    for (i = 0; i < _n; i++){
      _xx[i] = P._xx[i];
    }
  }
  else{
    _xx = NULL;
  }
}


/***** Assignment operator *****/
PredictorPoiss&
PredictorPoiss::operator=(const PredictorPoiss &P)
{
  int i;

  if (_n && _nTheta){
    delete [] _xx;
  }

  _p = P._p;
  _n = P._n;
  _nTheta = P._nTheta;
  _work = P._work;

  _Theta       = P._Theta;
  _PropTheta   = P._PropTheta;
  _PriorMean   = P._PriorMean;
  _PriorInvVar = P._PriorInvVar;
  _PriorWMean  = P._PriorWMean;

  _PropMean = P._PropMean;
  _U = P._U;
  _I = P._I;

  _x = P._x;

  _eta     = P._eta;
  _Propeta = P._Propeta;

  if (_nTheta && _n){
    _xx = new MatrixLT<double>[_n];
    if (!_xx) throw returnR("Out of memory in PredictorPoiss.cpp: PredictorPoiss::operator=", 99);
    for (i = 0; i < _n; i++){
      _xx[i] = P._xx[i];
    }
  }
  else{
    _xx = NULL;
  }

  return *this;
}


/***** Destructor *****/
PredictorPoiss::~PredictorPoiss()
{
  if (_n && _nTheta){
    delete [] _xx;
  }
}


/* ************************************************************************************************************************ */
/* Useful utilities                                                                                                         */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* change_Theta:  Change _Theta and update _eta                                                                             */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
//
//  newTheta[_nTheta]:   New value of theta
//
void
PredictorPoiss::change_Theta(const double *newTheta)
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

    /*** Update _eta  ***/
    if (_n){
      if (_p){ 
        _eta.bA(&_Theta, &_x, 0);
      }
    }
  }

  return;
}


/* ************************************************************************************************************************ */
/* Model related functions                                                                                                  */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* ll_poisson:  Apply either ll_poisson function to the objects of PredictorPoiss                                           */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
void
PredictorPoiss::ll_poisson(double *ll,            double *mu,
			   const double *offset,  const int *y,  const double *log_y_factor,  const int &order)
{
  Ll_poisson::ll_poisson(ll, _U.a(), _I.a(), _eta.a(), mu, offset, _Theta.aconst(), y,  log_y_factor, 
                         _x.aconst(), _xx, _n, _p, order);
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
/*            accept:    OUTPUT: indicator of acceptance (1 if accepted, 0 if not accepted)                                 */
/*                ll:     INPUT: current value of the log-likelihood                                                        */
/*                     OUTPUT: new value of the log-likelihood                                                              */
/*             mu[n]:   INPUT: current values of expected counts                                                            */
/*                     OUTPUT: new values of expected counts                                                                */
/*         Propmu[n]:  working space                                                                                        */
/*         offset[n]:  offset vector                                                                                        */
/*              y[n]:  observed counts                                                                                      */
/*   log_y_factor[n]:  vector of logarithms of factorials of the observed counts                                            */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
void
PredictorPoiss::update1(int *accept,           double *ll, 
                        double *mu,            double *Propmu,
                        const double *offset,  const int *y,     const double *log_y_factor)
{
  mcmc_common::update_reg_gamermanPoiss(accept, ll, _U.a(), _I.a(), 
				   _eta.a(), _Propeta.a(), mu, Propmu, _work.a(), 
                                   offset, _Theta.a(), _PropTheta.a(), y, log_y_factor, _x.aconst(), _xx, 
                                   _n, _p, _PriorMean.aconst(), _PriorWMean.aconst(), _PriorInvVar.aconst(), NULL, true, _PropMean.a(),
                                   "PredictorPoiss::update1");

  return;
}


/* ************************************************************************************************************************ */
/* Utilities                                                                                                                */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
void
PredictorPoiss::print(const int &colmax) const
{
  int cmax = (colmax < 0 ? _n : (colmax <= _n ? colmax : _n));

  Rprintf("\nObject of class PredictorPoiss:\n");
  Rprintf("=================================\n");
  Rprintf("p=%d,  nTheta=%d\n", _p, _nTheta);
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

  Rprintf("\nDESIGN MATRICES:\n");
  Rprintf("* X:\n");
  _x.print(0, cmax);

  Rprintf("\n");
  return;
}
