/*** RandomPoiss.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//               CREATED:  13/02/2007
//
//                                            BASICS:  13/02/2007
//                                        ll_poisson:  13/02/2007
//
//                                         updateRE1:  14/02/2007
//
//                                         updateRE2:  14/02/2007
//
//                                       updateRE2Bi:  29/03/2007
//
//                                   invert_REInvVar:  14/02/2007
//                              invert_REInvVarFromL:  14/02/2007
//                                             print:  14/02/2007
//
// PURPOSE: Class to hold parameters related to the random effects in the poisson log-linear GALMM
//
//          * it also holds design matrices and current values of the predictor itself  
//          * it implements also MCMC updates of these parameters
//
/* ********************************************************************************* */

#include "RandomPoiss.h"


/* ********************************************************************************* */
/* Constructors and destructors                                                      */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */

/***** Nonparametric constructor *****/
RandomPoiss::RandomPoiss()
  : _p(0), _nTheta(0), _N(0), _max_ni(0), _n(0), _REdist(mcmc_Random::_None), _prior_for_REMean(0), _prior_for_REInvVar(0),
    _xx(NULL)
{
}


/***** Parametric constructor 1 *****/
/*                                                                                                             */
/*                                                                                                             */
/* mean_ivar[p+LT(p)]:  initial values of the mean and inverse variance of random effects                      */
/*         mean_ivar[p] = initial values of the mean of random effects                                         */
/*    mean_ivar[+LT(p)] = initial value of the inverse variance of random effects (lower triangle)             */
/*                                                                                                             */
/* mean_ivarPrior[2]:                                                                                          */
/*                     mean_ivarPrior[0] = _prior_for_mean                                                     */
/*                     mean_ivarPrior[1] = _prior_for_ivar                                                     */
/*                                                                                                             */
/* meanPrior[2*p]: parameters of the prior of the means of random effects                                      */
/*               meanPrior[p] = prior means for RE means                                                       */
/*            meanPrior[+p] = prior inverse variances for RE means                                             */
/*                                                                                                             */
/* ivarPrior: parameters of the Wishart prior of the inverse variance matrix of random effects                 */
/*     _prior_for_ivar = _Wishart:                                                                             */
/*                  ivarPrior[0] = degrees of freedom                                                          */
/*             ivarPrior[+LT(p)] = inverse scale matrix (lower triangle)                                       */
/*                                                                                                             */
/*     _prior_for_ivar = _SDUnif:                                                                              */
/*                  ivarPrior[p]: upper limits of the uniform distribution for each standard deviation         */
/*                                                                                                             */
/*     _prior_for_ivar = _GammaIndep:                                                                          */
/*                  ivarPrior[p]: shape parameters                                                             */
/*                 ivarPrior[+p]: rate (inverse scale) parameters                                              */
/*                                                                                                             */
/*                                                                                                             */
RandomPoiss::RandomPoiss(const int &p,           const int &N,             const int *ni,
                         const double *theta,    const double *mean_ivar,
                         const int &REdist,      const int *mean_ivarPrior,  const double *meanPrior,  const double *ivarPrior, 
                         const double *x)
{
  int i, j;
  const double *xP, *cdP;
  double *dP;

  int LTp;

  /*** _p, _N, _ni, _max_ni, _n ***/
  /*** ======================== ***/
  _p = p;
  _N = N;

  if (_N){
    _ni = MatrixRect<int>(1, _N, ni);
    if (_ni.anyNonNeg()) throw returnR("Error in RandomPoiss.cpp: RandomPoiss::RandomPoiss(...), N > 0 & any ni <= 0", 1);
    _max_ni = _ni.max();
    _n = _ni.sum();    
  }
  else{
    _max_ni = 0;
    _n = 0;
  }

  /*** _nTheta, _work ***/
  /*** ================== ***/
  _nTheta = _p;
  LTp = (_p*(1+_p))/2;
  if (LTp) _work = MatrixRect<double>(1, LTp);
  else     _work = MatrixRect<double>(1, 1);

  /*** _REdist ***/
  /*** ======= ***/
  switch (REdist){
  case mcmc_Random::_None:
  case mcmc_Random::_Normal:
  case mcmc_Random::_Gspline:
    _REdist = REdist;
    break;
  default:
    throw returnR("Error in RandomPoiss.cpp: RandomPoiss::RandomPoiss(...), incorrect REdist argument", 1);
  }

  if (_nTheta){
    /*** _ThetaBar, _PropTheta, _PropMean ***/
    /*** ================================ ***/
    _ThetaBar = MatrixRect<double>(1, _nTheta);
    _PropTheta = MatrixRect<double>(1, _nTheta);
    _PropMean = MatrixRect<double>(1, _nTheta);

    /*** _REMean, _REInvVar, _REInvVarL, _REVar, _REWMean, _U, _I ***/
    /*** ======================================================== ***/
    _REMean = MatrixRect<double>(1, _nTheta, mean_ivar);
    _REInvVar = MatrixLT<double>(_nTheta, mean_ivar+_nTheta, 0);
    if (mean_ivarPrior[1] == mcmc_Random::_SDUnif || mean_ivarPrior[1] == mcmc_Random::_GammaIndep){   /*** _REInvVar should be diagonal ***/
      dP = _REInvVar.a();
      for (j = 0; j < _nTheta; j++){
        dP++;
        for (i = j+1; i < _nTheta; i++){
          *dP = 0.0;
          dP++;
        }
      }
    }
    _REInvVarL = _REInvVar;
    i = _REInvVarL.cholesky(0);
    if (i < _nTheta){
      Rprintf("WARNING: RandomPoiss.cpp: RandomPoiss::RandomPoiss(...), supplied _REInvVar is not of full rank\n");
      Rprintf("_REInvVar:\n");
      _REInvVar.print(0);
    }
    _REVar = _REInvVarL;
    _REVar.chinv(0);
    _REStdDev = MatrixRect<double>(1, _nTheta);
    _REVar.sqrtDiag(_REStdDev.a());
    
    _REWMean = MatrixRect<double>(1, _nTheta);
    Ab2(_REWMean.a(), &_REInvVar, _REMean.aconst());

    _U = MatrixRect<double>(1, _nTheta);
    _I = MatrixLT<double>(_nTheta);  

    /*** _prior_for_REMean, _REMeanPriorMean, _REMeanPriorInvVar, _REMeanPriorWMean ***/
    /*** ========================================================================== ***/
    switch (mean_ivarPrior[0]){
    case mcmc_Random::_Fixed_:
      _prior_for_REMean = mcmc_Random::_Fixed_;
      _REMeanPriorMean = MatrixRect<double>(1, _nTheta);
      _REMeanPriorInvVar = MatrixRect<double>(1, _nTheta);
      _REMeanPriorWMean = MatrixRect<double>(1, _nTheta);
      break;
    case mcmc_Random::_Normal_:
      _prior_for_REMean = mcmc_Random::_Normal_;
      _REMeanPriorMean = MatrixRect<double>(1, _nTheta, meanPrior);
      _REMeanPriorInvVar = MatrixRect<double>(1, _nTheta, meanPrior+_nTheta);
      _REMeanPriorWMean = MatrixRect<double>(1, _nTheta, meanPrior+_nTheta, meanPrior, 1);
      break;
    default:
      throw returnR("Error in RandomPoiss.cpp: RandomPoiss::RandomPoiss(...), incorrect mean_ivarPrior argument", 1);
    }

    /*** _prior_for_REInvVar, _REInvVarPriorDF, _REInvVarPriorInvScale, _work_rwishart, _work_invVarSlice ***/
    /*** ================================================================================================ ***/
    switch (mean_ivarPrior[1]){
    case mcmc_Random::_Fixed:
      _prior_for_REInvVar = mcmc_Random::_Fixed;
      _REInvVarPriorDF = MatrixRect<double>(1, 1);
      _REInvVarPriorInvScale = MatrixLT<double>(_nTheta);
      break;
    case mcmc_Random::_Wishart:
      _prior_for_REInvVar = mcmc_Random::_Wishart;
      _REInvVarPriorDF = MatrixRect<double>(1, 1, ivarPrior);
      _REInvVarPriorInvScale = MatrixLT<double>(_nTheta, ivarPrior+1, 0);
      break;
    case  mcmc_Random::_SDUnif:
      _prior_for_REInvVar = mcmc_Random::_SDUnif;
      _REInvVarPriorDF = MatrixRect<double>(1, 1);
      _REInvVarPriorInvScale = MatrixLT<double>(_nTheta);
      cdP = ivarPrior;
      dP = _REInvVarPriorInvScale.a();
      for (i = 0; i < _nTheta; i++){
        *dP = 1/((*cdP)*(*cdP));
        dP++;
        cdP++;
      }
      break;
    case mcmc_Random::_GammaIndep:
      _prior_for_REInvVar = mcmc_Random::_GammaIndep;
      _REInvVarPriorDF = MatrixRect<double>(1, _nTheta, ivarPrior);
      _REInvVarPriorInvScale = MatrixLT<double>(_nTheta);
      cdP = ivarPrior + _nTheta;
      dP = _REInvVarPriorInvScale.a();
      for (i = 0; i < _nTheta; i++){
        *dP = *cdP;
        dP++;
        cdP++;
      }
      break;
    default:
      throw returnR("Error in RandomPoiss.cpp: RandomPoiss::RandomPoiss(...), incorrect mean_ivarPrior argument", 1);
    }
 
     _work_rwishart = MatrixRect<double>(1, 2*_nTheta*_nTheta); 
     //_work_rwishart = MatrixRect<double>(1, 2*_I.length() + _nTheta*_nTheta);    /* changed on 12/01/2007 */
    _work_invVarSlice = MatrixRect<double>(1, 4*_nTheta);

    if (_N){
      /*** _Theta, _Theta_REMean ***/
      /*** ===================== ***/
      _Theta = MatrixRect<double>(_nTheta, _N, theta);
      _Theta.mean(_ThetaBar.a(), 1);
      
      _Theta_REMean = MatrixRect<double>(_nTheta, _N, theta);
      NSampleVar(&_I, _Theta_REMean.a(), _Theta.aconst(), _REMean.aconst(), _N);

      /*** _x,_eta, _Propeta ***/  
      /*** ================= ***/
      _x = MatrixRect<double>(_p, _n, x);

      _eta     = MatrixRect<double>(1, _n);
      _Propeta = MatrixRect<double>(1, _max_ni);
      if (_p){
	_eta.BAcolProd(&_Theta, &_ni, &_x, 0);
      }

      /*** _xx ***/
      /*** ====***/
      _xx = new MatrixLT<double>[_n];
      if (!_xx) throw returnR("Out of memory in RandomPoiss.cpp: RandomPoiss::RandomPoiss(...)", 99);

      xP = x;
      for (i = 0; i < _n; i++){
        _xx[i] = MatrixLT<double>(_p, xP);

        xP += _p;
      }
    }
  }
  else{    /*** else from if (_nTheta) ***/
    if (_N) _eta = MatrixRect<double>(1, _n);
  }

  if (!_nTheta || !_n){
    _xx = NULL;
  }
}


/***** Copy constructor *****/
RandomPoiss::RandomPoiss(const RandomPoiss &P)
{
  int i;

  _p = P._p;
  _nTheta = P._nTheta;
  _N = P._N;
  _ni = P._ni;
  _max_ni = P._max_ni;
  _n = P._n;

  _REdist = P._REdist;
  _prior_for_REMean = P._prior_for_REMean;
  _prior_for_REInvVar = P._prior_for_REInvVar;

  _Theta = P._Theta;
  _ThetaBar = P._ThetaBar;
  _Theta_REMean = P._Theta_REMean;
  _PropTheta = P._PropTheta;

  _REMean = P._REMean;
  _REInvVar = P._REInvVar;
  _REInvVarL = P._REInvVarL;
  _REVar = P._REVar;
  _REStdDev = P._REStdDev;
  _REWMean = P._REWMean;

  _PropMean = P._PropMean;
  _U = P._U;
  _I = P._I;

  _REMeanPriorMean = P._REMeanPriorMean;
  _REMeanPriorInvVar = P._REMeanPriorInvVar;
  _REMeanPriorWMean = P._REMeanPriorWMean;

  _REInvVarPriorDF = P._REInvVarPriorDF;
  _REInvVarPriorInvScale = P._REInvVarPriorInvScale;

  _eta = P._eta;
  _Propeta = P._Propeta;

  _x = P._x;

  _work = P._work;
  _work_rwishart = P._work_rwishart;
  _work_invVarSlice = P._work_invVarSlice;

  if (_nTheta && _n){
    _xx = new MatrixLT<double>[_n];
    if (!_xx) throw returnR("Out of memory in RandomPoiss.cpp: RandomPoiss::RandomPoiss(P)", 99);
    for (i = 0; i < _n; i++){
      _xx[i] = P._xx[i];
    }
  }
  else{
    _xx = NULL;
  }
}


/***** Assignment operator *****/
RandomPoiss&
RandomPoiss::operator=(const RandomPoiss &P)
{
  int i;

  if (_n && _nTheta){
    delete [] _xx;
  }

  _p = P._p;
  _nTheta = P._nTheta;
  _N = P._N;
  _ni = P._ni;
  _max_ni = P._max_ni;
  _n = P._n;

  _REdist = P._REdist;
  _prior_for_REMean = P._prior_for_REMean;
  _prior_for_REInvVar = P._prior_for_REInvVar;

  _Theta = P._Theta;
  _ThetaBar = P._ThetaBar;
  _Theta_REMean = P._Theta_REMean;
  _PropTheta = P._PropTheta;

  _REMean = P._REMean;
  _REInvVar = P._REInvVar;
  _REVar = P._REVar;
  _REStdDev = P._REStdDev;
  _REInvVarL = P._REInvVarL;
  _REWMean = P._REWMean;

  _PropMean = P._PropMean;
  _U = P._U;
  _I = P._I;

  _REMeanPriorMean = P._REMeanPriorMean;
  _REMeanPriorInvVar = P._REMeanPriorInvVar;
  _REMeanPriorWMean = P._REMeanPriorWMean;

  _REInvVarPriorDF = P._REInvVarPriorDF;
  _REInvVarPriorInvScale = P._REInvVarPriorInvScale;

  _eta = P._eta;
  _Propeta = P._Propeta;

  _x = P._x;

  _work = P._work;
  _work_rwishart = P._work_rwishart;
  _work_invVarSlice = P._work_invVarSlice;

  if (_nTheta && _n){
    _xx = new MatrixLT<double>[_n];
    if (!_xx) throw returnR("Out of memory in RandomPoiss.cpp: RandomPoiss::operator=", 99);
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
RandomPoiss::~RandomPoiss()
{
  if (_n && _nTheta){
    delete [] _xx;
  }
}


/* ************************************************************************************************************************ */
/* Model related functions                                                                                                  */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* ll_poisson:  Apply ll_poisson to the objects of RandomPoiss                                                              */
/*              (loop over clusters)                                                                                        */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/*         mu:  vector of length _n                                                                                         */
/*     offset:  vector of length _n                                                                                         */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
void
RandomPoiss::ll_poisson(double *ll,            double *mu,
			const double *offset,  const int *y,  const double *log_y_factor,  const int &order)
{
  static int i, j;
  static double llCl;
  const int *yP, *niP;
  const double *thetaP, *xP, *offsetP, *log_y_factorP;
  const MatrixLT<double> *xxP;
  double *muP, *etaP;

  *ll           = 0;
  yP            = y;
  log_y_factorP = log_y_factor;
  muP           = mu;
  thetaP        = _Theta.aconst();
  etaP          = _eta.a();
  offsetP       = offset;
  xP            = _x.aconst();
  xxP           = _xx;
  niP           = _ni.aconst(); 

  for (i = 0; i < _N; i++){
    Ll_poisson::ll_poisson(&llCl, _U.a(), _I.a(), etaP, muP, offsetP, thetaP, yP,  log_y_factorP, 
                           xP, xxP, *niP, _p, order);
    *ll += llCl; 

    yP            += (*niP);
    log_y_factorP += (*niP);
    muP           += (*niP);
    thetaP        += _nTheta;
    etaP          += (*niP);
    offsetP       += (*niP);
    xP            += _p * (*niP);
    xxP           += (*niP);
    niP++;
  }

  return;
}


/* ************************************************************************************************************************ */
/* MCMC related functions                                                                                                   */
/*   for model with NORMAL random effects                                                                                   */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* updateRE1: Update of normally distributed random effects                                                                 */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/*        accept[_N]:  OUTPUT: indicators of acceptance (1 if accepted, 0 if not accepted)                                  */
/*                ll:   INPUT: current value of the log-likelihood                                                          */
/*                     OUTPUT: new value of the log-likelihood                                                              */
/*            mu[_n]:   INPUT: current values of expected counts                                                            */
/*                     OUTPUT: new values of expected counts                                                                */
/*        Propmu[_n]:  working space                                                                                        */
/*                     * it is sufficient if its length is max(_ni)                                                         */
/*        offset[_n]:  offset vector                                                                                        */
/*             y[_n]:  vector of observations                                                                               */
/*  log_y_factor[_n]:  vector with log(y!)                                                                                  */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
void
RandomPoiss::updateRE1(int *accept,           double *ll, 
                       double *mu,            double *Propmu,
                       const double *offset,  const int *y,   const double *log_y_factor)
{
  static int i;
  static double llCl;
  const int *yP, *niP;
  const double *xP, *offsetP, *log_y_factorP;
  const MatrixLT<double> *xxP;
  int *acceptP;
  double *muP, *thetaP, *etaP;

  /*** Update _REWMean = D^{-1}*Mean ***/
  Ab2(_REWMean.a(), &_REInvVar, _REMean.aconst());

  /*** Update _Theta: Loop over clusters ***/
  *ll           = 0;
  acceptP       = accept;
  yP            = y;
  log_y_factorP = log_y_factor;
  muP           = mu;
  thetaP        = _Theta.a();
  etaP          = _eta.a();
  offsetP       = offset;
  xP            = _x.aconst();
  xxP           = _xx;
  niP           = _ni.aconst(); 

  //Rprintf("\nREInvVar=%g,  REInvVarL=%g,  REMean=%g,  REWMean=%g", _REInvVar.aconst()[0], _REInvVarL.aconst()[0], _REMean.aconst()[0], _REWMean.aconst()[0]);
  for (i = 0; i < _N; i++){
    //if (i < 20){Rprintf("\nb(%d)=%g", i, *thetaP); }
    mcmc_common::update_reg_gamermanPoiss(acceptP, &llCl, _U.a(), _I.a(), 
		                          etaP, _Propeta.a(), muP, Propmu, _work.a(),
 		                          offsetP, thetaP, _PropTheta.a(), yP, log_y_factorP, xP, xxP,
      	                                  *niP, _p, _REMean.aconst(), _REWMean.aconst(), 
                                          _REInvVar.aconst(), _REInvVarL.aconst(), false, _PropMean.a(),
                                          "RandomPoiss::updateRE1");
    //if (i < 20){Rprintf(",  new b(%d)=%g", i, *thetaP); }
    *ll += llCl; 
    
    acceptP++;
    yP            += (*niP);
    log_y_factorP += (*niP);
    muP           += (*niP);
    thetaP        += _nTheta;
    etaP          += (*niP);
    offsetP       += (*niP);
    xP            += _p * (*niP);
    xxP           += (*niP);
    niP++;
  }

  return;
}


/* ************************************************************************************************************************ */
/* MCMC related functions                                                                                                   */
/*   for model with UNIVARIATE G-spline random effects                                                                      */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* updateRE2: Update of G-spline distributed random effects                                                                 */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
void
RandomPoiss::updateRE2(int *accept,           double *ll, 
                       double *mu,            double *Propmu,
                       const Gspline2 *Gspl,
                       const double *offset,  const int *y,    const double *log_y_factor,
                       const MatrixRect<int> *alloc)
{
  static int i;
  static double llCl;
  const int *yP, *niP, *allocP, *KP;
  const double *log_y_factorP, *xP, *offsetP, *REMeanP, *inv_sigma2_d2P, *d_knotsP;
  const MatrixLT<double> *xxP;
  const MatrixRect<double> *d_knotsM;
  int *acceptP;
  double *muP, *thetaP, *etaP, *meanP, *REWMeanP;

  *ll           = 0;
  acceptP       = accept;
  yP            = y;
  log_y_factorP = log_y_factor;
  muP           = mu;
  thetaP        = _Theta.a();
  etaP          = _eta.a();
  offsetP       = offset;
  xP            = _x.aconst();
  xxP           = _xx;
  niP           = _ni.aconst(); 
  allocP        = alloc->aconst();

  REMeanP  = _REMean.aconst();                 /*** constant over cycles of the loop over observations ***/
  meanP    = _ThetaBar.a();                    /*** working space at each cycle of the loop            ***/
  REWMeanP = _REWMean.a();                     /*** working space at each cycle of the loop            ***/

  switch(_nTheta){
  case 0:
    return;

  /*** UNIVARIATE random effects ***/
  /*** ========================= ***/
  case 1:
    KP             = Gspl->KAconst();
    d_knotsM       = Gspl->d_knotsconst();
    d_knotsP       = d_knotsM->aconst();
    inv_sigma2_d2P = Gspl->inv_sigma2_d2Aconst();

    //Rprintf("d_knotsP: "); d_knotsM->print(0);
    //Rprintf("REMeanP: %g\n\n", *REMeanP);
    for (i = 0; i < _N; i++){

      /*** Mean of the random effect given allocation, store it in _ThetaBar. ***/
      *meanP = *REMeanP + d_knotsP[*allocP + (*KP)];

      /*** Update _REWMean = (d*sigma)^{-2}*(intcpt + d*knot[alloc]).  ***/
      *REWMeanP = *inv_sigma2_d2P * (*meanP);

      /*** Sample new value of the random effect ***/
      mcmc_common::update_reg_gamermanPoiss(acceptP, &llCl, _U.a(), _I.a(), 
                                            etaP, _Propeta.a(), muP, Propmu, _work.a(),
		                            offsetP, thetaP, _PropTheta.a(), yP, log_y_factorP, xP, xxP,
     	                                    *niP, _p, meanP, REWMeanP, inv_sigma2_d2P, NULL, true, _PropMean.a(),
                                            "RandomPoiss::updateRE2");
      *ll += llCl; 
     
      /*** Increase pointers ***/
      acceptP++;
      yP            += (*niP);
      log_y_factorP += (*niP);
      muP           += (*niP);
      thetaP        += _nTheta;
      etaP          += (*niP);
      offsetP       += (*niP);
      xP            += _p * (*niP);
      xxP           += (*niP);
      niP++;
      allocP        += _nTheta;
    }

    return;

  /*** BIVARIATE random effects (copula???) ***/
  /*** ==================================== ***/
  case 2:
    throw returnR("Error in RandomPoiss.cpp: RandomPoiss::updateRE2. Not implemented for _nTheta = 2", 1);
    return;

  /*** MULTI(>2)VARIATE random effects ***/
  /*** =============================== ***/
  default:
    throw returnR("Error in RandomPoiss.cpp: RandomPoiss::updateRE2. Not implemented for _nTheta > 2", 1);
    return;
  }
}


/* ************************************************************************************************************************ */
/* MCMC related functions                                                                                                   */
/*   for model with BIVARIATE G-spline random effects                                                                       */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* updateRE2Bi: Update of BIVARIATE G-spline distributed random effects                                                     */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
//
// alloc[this->_N]:       vector of single index (0,...,Gspl->total_length-1) allocations
//
void
RandomPoiss::updateRE2Bi(int *accept,             double *ll, 
                         double *mu,              double *Propmu,
                         const BiGspline2 *Gspl,
                         const double *offset,    const int *y,    const double *log_y_factor,
                         const MatrixRect<int> *alloc)
{
  static int i, length0, i0, i1;
  static double llCl;
  const int *yP, *niP, *allocP;
  const double *log_y_factorP, *xP, *offsetP, *REMeanP, *inv_sigma2_d2P, *d_knotsP0, *d_knotsP1;
  const MatrixLT<double> *xxP;
  int *acceptP;
  double *muP, *thetaP, *etaP, *meanP, *REWMeanP;

  if (_nTheta != BiGspline2A::_dim){
    REprintf("_nTheta=%d,  BiGspline2A::_dim=%d\n", _nTheta, BiGspline2A::_dim);
    throw returnR("Error in RandomPoiss.cpp: updateRE2Bi(). Not implemented for this dimension.", 1);
  }

  *ll           = 0;
  acceptP       = accept;
  yP            = y;
  log_y_factorP = log_y_factor;
  muP           = mu;
  thetaP        = _Theta.a();
  etaP          = _eta.a();
  offsetP       = offset;
  xP            = _x.aconst();
  xxP           = _xx;
  niP           = _ni.aconst(); 
  allocP        = alloc->aconst();

  REMeanP  = _REMean.aconst();                 /*** constant over cycles of the loop over observations ***/
  meanP    = _ThetaBar.a();                    /*** working space at each cycle of the loop            ***/
  REWMeanP = _REWMean.a();                     /*** working space at each cycle of the loop            ***/

  length0        = Gspl->lengthconst()[0];
  d_knotsP0      = Gspl->d_knotsconst()[0].aconst();
  d_knotsP1      = Gspl->d_knotsconst()[1].aconst();
  inv_sigma2_d2P = Gspl->inv_sigma2_d2const();

  for (i = 0; i < _N; i++){

    /*** Mean of the random effect given allocation, store it in _ThetaBar. ***/
    i0 = *allocP % length0;
    i1 = *allocP / length0;
    meanP[0] = REMeanP[0] + d_knotsP0[i0];
    meanP[1] = REMeanP[1] + d_knotsP1[i1];

    /*** Update _REWMean = (d*sigma)^{-2}*(intcpt + d*knot[alloc]).  ***/
    REWMeanP[0] = inv_sigma2_d2P[0] * meanP[0];
    REWMeanP[1] = inv_sigma2_d2P[1] * meanP[1];

    /*** Sample new value of the random effect ***/
    mcmc_common::update_reg_gamermanPoiss(acceptP, &llCl, _U.a(), _I.a(), 
                                          etaP, _Propeta.a(), muP, Propmu, _work.a(),
       	                                  offsetP, thetaP, _PropTheta.a(), yP, log_y_factorP, xP, xxP,
   	                                  *niP, _p, meanP, REWMeanP, inv_sigma2_d2P, NULL, true, _PropMean.a(),
                                          "RandomPoiss::updateRE2Bi");
    *ll += llCl; 
   
    /*** Increase pointers ***/
    acceptP++;
    yP            += (*niP);
    log_y_factorP += (*niP);
    muP           += (*niP);
    thetaP        += _nTheta;
    etaP          += (*niP);
    offsetP       += (*niP);
    xP            += _p * (*niP);
    xxP           += (*niP);
    niP++;
    allocP++;
  }

  return;
}


/* ************************************************************************************************************************ */
/* Utilities                                                                                                                */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* invert_REInvVar:  Invert _REInvVar and store it in _REVar                                                                */
/*                                                                                                                          */
/*     RETURN: rank of _REInvVar (as described in chol_inv function)                                                        */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
int
RandomPoiss::invert_REInvVar()
{
  int Rank;

  _REVar.array2mat(_REInvVar.aconst());
  Rank = _REVar.chol_inv(0, NULL);

  return Rank;
}

void
RandomPoiss::invert_REInvVarFromL()
{
  _REVar.array2mat(_REInvVarL.aconst());
  _REVar.chinv(0);

  return;
}


/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* print                                                                                                                    */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
void
RandomPoiss::print(const int &colmax) const
{
  int cmax = (colmax < 0 ? _N : (colmax <= _N ? colmax : _N));

  Rprintf("\nObject of class RandomPoiss:\n");
  Rprintf("==============================\n");
  Rprintf("p=%d,  nTheta=%d\n", _p, _nTheta);
  Rprintf("N=%d,  n=%d,  max(ni)=%d\n",  _N, _n, _max_ni);
  Rprintf("ni: ");
  _ni.printI(0, cmax);
  Rprintf("RE distribution: %s\n", (_REdist == mcmc_Random::_None ? "None" : (_REdist==mcmc_Random::_Normal ? "Normal" : (_REdist==mcmc_Random::_Gspline ? "G-spline" : "Unknown"))));


  Rprintf("\nTHETA:\n");
  _Theta.print(0, cmax);
  
  Rprintf("\nMOMENTS OF THETA:\n");
  Rprintf("* REMean: ");
  _REMean.print(0);
  Rprintf("* REInvVar:\n");
  _REInvVar.print(0);
  Rprintf("* REInvVarL:\n");
  _REInvVarL.print(0);

  Rprintf("\nPRIOR FOR MEANS:\n");
  Rprintf("* REMeanPriorMean: ");
  _REMeanPriorMean.print(0);
  Rprintf("* REMeanPriorInvVar: ");
  _REMeanPriorInvVar.print(0);
  Rprintf("* REMeanPriorWMean: ");
  _REMeanPriorWMean.print(0);

  Rprintf("\nPRIOR FOR INVERSE VARIANCE MATRIX:\n");
  Rprintf("* REInvVarPriorDF:\n");
  _REInvVarPriorDF.print(0);
  Rprintf("* REInvVarPriorInvScale:\n");
  _REInvVarPriorInvScale.print(0);

  Rprintf("\nLINEAR PREDICTORS:\n");
  Rprintf("* eta:\n");
  _eta.print(0, cmax);

  Rprintf("\nDESIGN MATRICES:\n");
  Rprintf("* X:\n");
  _x.print(0, cmax);

  Rprintf("\n");
  return;
}
