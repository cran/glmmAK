/*** RandomCL.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//               CREATED:  17/09/2006
//   MAJOR MODIFICATIONS:  13/10/2006  (order of 'p' and 'q' parameters in _Theta changed)
//
//                                            BASICS:  17/09/2006
//                                      ll_cumlogit2:  27/09/2006
//
//                                         updateRE1:  19/09/2006
//
//                                         updateRE2:  29/10/2006
//
//                                       updateRE2Bi:  29/03/2007
//
//                                   invert_REInvVar:  27/09/2006 
//                              invert_REInvVarFromL:  02/10/2006
//                                             print:  27/09/2006
//
// PURPOSE: Class to hold parameters related to the random effects in the cumulative logits GALMM
//
//          * it also holds design matrices and current values of the predictor itself  
//          * it implements also MCMC updates of these parameters
//
/* ********************************************************************************* */

#include "RandomCL.h"


/* ********************************************************************************* */
/* Constructors and destructors                                                      */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */

/***** Nonparametric constructor *****/
RandomCL::RandomCL()
  : _p(0), _q(0), _nTheta(0), _N(0), _max_ni(0), _n(0), _C(1), _REdist(mcmc_Random::_None), _prior_for_REMean(0), _prior_for_REInvVar(0),
    _xx(NULL), _vv(NULL), _xv(NULL)
{
}


/***** Parametric constructor 1 *****/
/*                                                                                                             */
/*                                                                                                             */
/* mean_ivar[(C*q+p)+LT(C*q+p)]:  initial values of the mean and inverse variance of random effects            */
/*         mean_ivar[C*q+p] = initial values of the mean of random effects                                     */
/*    mean_ivar[+LT(C*q+p)] = initial value of the inverse variance of random effects (lower triangle)         */
/*                                                                                                             */
/* mean_ivarPrior[2]:                                                                                          */
/*                     mean_ivarPrior[0] = _prior_for_mean                                                     */
/*                     mean_ivarPrior[1] = _prior_for_ivar                                                     */
/*                                                                                                             */
/* meanPrior[2*(p+C*q)]: parameters of the prior of the means of random effects                                */
/*               meanPrior[C*q+p] = prior means for RE means                                                   */
/*            meanPrior[+(C*q+p)] = prior inverse variances for RE means                                       */
/*                                                                                                             */
/* ivarPrior[]: parameters of the Wishart prior of the inverse variance matrix of random effects               */
/*     _prior_for_ivar = _Wishart:                                                                             */
/*                  ivarPrior[0] = degrees of freedom                                                          */
/*         ivarPrior[+LT(C*q+p)] = inverse scale matrix (lower triangle)                                       */
/*                                                                                                             */
/*     _prior_for_ivar = _SDUnif:                                                                              */
/*                  ivarPrior[C*q+p]: upper limits of the uniform distribution for each standard deviation     */
/*                                                                                                             */
/*     _prior_for_ivar = _GammaIndep:                                                                          */
/*                  ivarPrior[p]: shape parameters                                                             */
/*                 ivarPrior[+p]: rate (inverse scale) parameters                                              */
/*                                                                                                             */
RandomCL::RandomCL(const int &p,           const int &q,               const int &N,             const int *ni,            const int &C,
                   const double *theta,    const double *mean_ivar,
                   const int &REdist,      const int *mean_ivarPrior,  const double *meanPrior,  const double *ivarPrior, 
                   const double *x,        const double *v)
{
  int i, j;
  const double *xP, *vP, *cdP;
  double *dP;

  /*** _p, _q, _N, _ni, _max_ni, _n ***/
  /*** ============================ ***/
  _p = p;
  _q = q;
  _N = N;

  if (_N){
    _ni = MatrixRect<int>(1, _N, ni);
    if (_ni.anyNonNeg()) throw returnR("Error in RandomCL.cpp: RandomCL::RandomCL(...), N > 0 & any ni <= 0", 1);
    _max_ni = _ni.max();
    _n = _ni.sum();    
  }
  else{
    _max_ni = 0;
    _n = 0;
  }

  /*** _C, _nTheta, _work ***/
  /*** ================== ***/
  if (C < 1) throw returnR("Error in RandomCL.cpp: RandomCL::RandomCL(...), C < 1", 1);
  _C = C;
  _nTheta = _C*_q + _p;
  _work = MatrixRect<double>(1, 5*_C);

  /*** _REdist ***/
  /*** ======= ***/
  switch (REdist){
  case mcmc_Random::_None:
  case mcmc_Random::_Normal:
  case mcmc_Random::_Gspline:
    _REdist = REdist;
    break;
  default:
    throw returnR("Error in RandomCL.cpp: RandomCL::RandomCL(...), incorrect REdist argument", 1);
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
      Rprintf("WARNING: RandomCL.cpp: RandomCL::RandomCL(...), supplied _REInvVar is not of full rank\n");
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
      throw returnR("Error in RandomCL.cpp: RandomCL::RandomCL(...), incorrect mean_ivarPrior argument", 1);
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
      throw returnR("Error in RandomCL.cpp: RandomCL::RandomCL(...), incorrect mean_ivarPrior argument", 1);
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

      /*** _x, _v, _etaX, _PropetaX, _etaV, _PropetaV, _eta, _Propeta ***/  
      /*** ========================================================== ***/
      _x = MatrixRect<double>(_p, _n, x);
      _v = MatrixRect<double>(_q, _n, v);

      _etaX = MatrixRect<double>(1, _n);
      _PropetaX = MatrixRect<double>(1, _max_ni);
      if (_p){
	_etaX.BAcolProd(&_Theta, &_ni, &_x, _C*_q);
      }

      _etaV = MatrixRect<double>(_C, _n);
      _PropetaV = MatrixRect<double>(_C, _max_ni);
      if (_q){
	_etaV.BAcolProd2(&_Theta, &_ni, &_v, 0);
      }

      _eta = MatrixRect<double>(_C, _n);
      _Propeta = MatrixRect<double>(_C, _max_ni);
      _eta.b_plus_rowsA(_etaX.aconst(), &_etaV);

      /*** _xx, _vv, _xv ***/
      /*** ============= ***/
      _xx = new MatrixLT<double>[_n];
      _vv = new MatrixLT<double>[_n];
      _xv = new MatrixRect<double>[_n];
      if (!_xx || !_vv || !_xv) throw returnR("Out of memory in RandomCL.cpp: RandomCL::RandomCL(...)", 99);

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
    if (_N) _eta = MatrixRect<double>(_C, _n);
  }

  if (!_nTheta || !_n){
    _xx = NULL;
    _vv = NULL;
    _xv = NULL;
  }
}


/***** Copy constructor *****/
RandomCL::RandomCL(const RandomCL &P)
{
  int i;

  _p = P._p;
  _q = P._q;
  _nTheta = P._nTheta;
  _N = P._N;
  _ni = P._ni;
  _max_ni = P._max_ni;
  _n = P._n;
  _C = P._C;

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

  _etaX = P._etaX;
  _etaV = P._etaV;
  _eta = P._eta;
  _PropetaX = P._PropetaX;
  _PropetaV = P._PropetaV;
  _Propeta = P._Propeta;

  _x = P._x;
  _v = P._v;

  _work = P._work;
  _work_rwishart = P._work_rwishart;
  _work_invVarSlice = P._work_invVarSlice;

  if (_nTheta && _n){
    _xx = new MatrixLT<double>[_n];
    _vv = new MatrixLT<double>[_n];
    _xv = new MatrixRect<double>[_n];
    if (!_xx || !_vv || !_xv) throw returnR("Out of memory in RandomCL.cpp: RandomCL::RandomCL(P)", 99);
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
RandomCL&
RandomCL::operator=(const RandomCL &P)
{
  int i;

  if (_n && _nTheta){
    delete [] _xx;
    delete [] _vv;
    delete [] _xv;
  }

  _p = P._p;
  _q = P._q;
  _nTheta = P._nTheta;
  _N = P._N;
  _ni = P._ni;
  _max_ni = P._max_ni;
  _n = P._n;
  _C = P._C;

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

  _etaX = P._etaX;
  _etaV = P._etaV;
  _eta = P._eta;
  _PropetaX = P._PropetaX;
  _PropetaV = P._PropetaV;
  _Propeta = P._Propeta;

  _x = P._x;
  _v = P._v;

  _work = P._work;
  _work_rwishart = P._work_rwishart;
  _work_invVarSlice = P._work_invVarSlice;

  if (_nTheta && _n){
    _xx = new MatrixLT<double>[_n];
    _vv = new MatrixLT<double>[_n];
    _xv = new MatrixRect<double>[_n];
    if (!_xx || !_vv || !_xv) throw returnR("Out of memory in RandomCL.cpp: RandomCL::operator=", 99);
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
RandomCL::~RandomCL()
{
  if (_n && _nTheta){
    delete [] _xx;
    delete [] _vv;
    delete [] _xv;
  }
}


/* ************************************************************************************************************************ */
/* Model related functions                                                                                                  */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* ll_cumlogit2:  Apply either ll_cumlogitFS2 or ll_cumlogitNR2 function to the objects of RandomC   L                      */
/*                (loop over clusters)                                                                                      */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/*       prob:  matrix (C+1) x _n                                                                                           */
/*     offset:  matrix C x _n                                                                                               */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
void
RandomCL::ll_cumlogit2(double *ll,            double *prob,   int *anyZero,  
	               const double *offset,  const int *y,   const int &update_etas,  const int &order,
                       void (*ll_2)(double*, MatrixRect<double>*, MatrixLT<double>*, 
                                    double*, double*, double*, double*, double*, int*,
                                    const double*, const double*, const int*, const double*, const double*,
                                    const MatrixLT<double>*, const MatrixLT<double>*, const MatrixRect<double>*,
                                    const int&, const int&, const int&, const int&, const int&, const int&)
                       )
{
  static int i, j;
  static double llCl;
  const int *yP, *niP;
  const double *xP, *vP, *offsetP;
  const MatrixLT<double> *xxP, *vvP;
  const MatrixRect<double> *xvP;
  double *probP, *thetaP, *etaXP, *etaVP, *etaP;
  int anyZeroCl = 0;

  *ll = 0;
  *anyZero = 0;
  yP = y;
  probP = prob;
  thetaP = _Theta.a();
  etaXP = _etaX.a();
  etaVP = _etaV.a();
  etaP = _eta.a();
  offsetP = offset;
  xP = _x.aconst();
  vP = _v.aconst();
  xxP = _xx;
  vvP = _vv;
  xvP = _xv;
  niP = _ni.aconst(); 

  for (i = 0; i < _N; i++){
    ll_2(&llCl, &_U, &_I, 
         etaXP,  etaVP, etaP, probP, _work.a(), &anyZeroCl,
         offsetP, thetaP, yP, xP, vP, 
         xxP, vvP, xvP,
         *niP, _C, _p, _q, update_etas, order);

    if (anyZeroCl) *anyZero = 1;
    *ll += llCl; 

    yP += (*niP);
    j = _C * (*niP);    
    probP += j + (*niP);
    thetaP += _nTheta;
    etaXP += (*niP);
    etaVP += j;
    etaP += j;
    offsetP += j;
    xP += _p * (*niP);
    vP += _q * (*niP);
    xxP += (*niP);
    vvP += (*niP);
    xvP += (*niP);
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
/*     accept:  OUTPUT: indicator of acceptance (1 if accepted, 0 if not accepted)                                          */
/*              * array of length _N (number of clusters)                                                                   */
/*         ll:   INPUT: current value of the log-likelihood                                                                 */
/*              OUTPUT: new value of the log-likelihood                                                                     */
/*       prob:   INPUT: current values of category probabilities                                                            */
/*              OUTPUT: new values of category probabilities                                                                */
/*              * matrix (C+1) x _n                                                                                         */
/*   Propprob:  working space                                                                                               */
/*              * matrix (C+1) x _n                                                                                         */
/*              * or, it is sufficient if its dimension is (C+1) x max(_ni)                                                 */
/*     offset:  offset matrix                                                                                               */
/*              * matrix C x _n                                                                                             */
/*          y:  vector of observations                                                                                      */
/*              * vector of length _n                                                                                       */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
void
RandomCL::updateRE1(int *accept,                       double *ll, 
                    MatrixRect<double> *prob,          MatrixRect<double> *Propprob,
                    const MatrixRect<double> *offset,  const MatrixRect<int> *y)
{
  static int i, j;
  static double llCl;
  const int *yP, *niP;
  const double *xP, *vP, *offsetP;
  const MatrixLT<double> *xxP, *vvP;
  const MatrixRect<double> *xvP;
  int *acceptP;
  double *probP, *thetaP, *etaXP, *etaVP, *etaP;

  /*** Update _REWMean = D^{-1}*Mean ***/
  Ab2(_REWMean.a(), &_REInvVar, _REMean.aconst());

  /*** Update _Theta: Loop over clusters ***/
  *ll = 0;
  acceptP = accept;
  yP = y->aconst();
  probP = prob->a();
  thetaP = _Theta.a();
  etaXP = _etaX.a();
  etaVP = _etaV.a();
  etaP = _eta.a();
  offsetP = offset->aconst();
  xP = _x.aconst();
  vP = _v.aconst();
  xxP = _xx;
  vvP = _vv;
  xvP = _xv;
  niP = _ni.aconst(); 

  for (i = 0; i < _N; i++){
    mcmc_common::update_reg_gamermanCL(acceptP, &llCl, &_U, &_I, 
                                       etaXP, _PropetaX.a(), etaVP, _PropetaV.a(), etaP, _Propeta.a(), probP, Propprob->a(), _work.a(),
                                       offsetP, thetaP, _PropTheta.a(), yP, xP, vP, xxP, vvP, xvP, 
                                       *niP, _C, _p, _q, _REMean.aconst(), _REWMean.aconst(), _REInvVar.aconst(), &_REInvVarL, false, 
                                       &_PropMean,
                                       ll_cumlogitFS2, ll_cumlogitNR2, "RandomCL::updateRE1");
    *ll += llCl; 
    
    acceptP++;
    yP += (*niP);
    j = _C * (*niP);
    probP += j + (*niP);
    thetaP += _nTheta;
    etaXP += (*niP);
    etaVP += j;
    etaP += j;
    offsetP += j;
    xP += _p * (*niP);
    vP += _q * (*niP);
    xxP += (*niP);
    vvP += (*niP);
    xvP += (*niP);
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
RandomCL::updateRE2(int *accept,                       double *ll, 
                    MatrixRect<double> *prob,          MatrixRect<double> *Propprob,
                    const Gspline2 *Gspl,
                    const MatrixRect<double> *offset,  const MatrixRect<int> *y,
                    const MatrixRect<int> *alloc)
{
  static int i, j;
  static double llCl;
  const int *yP, *niP, *allocP, *KP;
  const double *xP, *vP, *offsetP, *REMeanP, *inv_sigma2_d2P, *d_knotsP;
  const MatrixLT<double> *xxP, *vvP;
  const MatrixRect<double> *xvP, *d_knotsM;
  int *acceptP;
  double *probP, *thetaP, *etaXP, *etaVP, *etaP, *meanP, *REWMeanP;

  *ll     = 0;
  acceptP = accept;
  yP      = y->aconst();
  probP   = prob->a();
  thetaP  = _Theta.a();
  etaXP   = _etaX.a();
  etaVP   = _etaV.a();
  etaP    = _eta.a();
  offsetP = offset->aconst();
  xP      = _x.aconst();
  vP      = _v.aconst();
  xxP     = _xx;
  vvP     = _vv;
  xvP     = _xv;
  niP     = _ni.aconst(); 
  allocP  = alloc->aconst();

  REMeanP  = _REMean.aconst();              /*** constant over cycles of the loop over observations ***/
  meanP    = _ThetaBar.a();                 /*** working space at each cycle of the loop            ***/
  REWMeanP = _REWMean.a();                  /*** working space at each cycle of the loop            ***/

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
      mcmc_common::update_reg_gamermanCL(acceptP, &llCl, &_U, &_I, 
                                         etaXP, _PropetaX.a(), etaVP, _PropetaV.a(), etaP, _Propeta.a(), probP, Propprob->a(), _work.a(),
                                         offsetP, thetaP, _PropTheta.a(), yP, xP, vP, xxP, vvP, xvP, 
                                         *niP, _C, _p, _q, 
                                         meanP, REWMeanP, inv_sigma2_d2P, NULL, true, &_PropMean,
                                         ll_cumlogitFS2, ll_cumlogitNR2, "RandomCL::updateRE2");
      *ll += llCl; 
     
      /*** Increase pointers ***/
      acceptP++;
      yP      += (*niP);
      j       = _C * (*niP);
      probP   += j + (*niP);
      thetaP  += _nTheta;
      etaXP   += (*niP);
      etaVP   += j;
      etaP    += j;
      offsetP += j;
      xP      += _p * (*niP);
      vP      += _q * (*niP);
      xxP     += (*niP);
      vvP     += (*niP);
      xvP     += (*niP);
      niP++;
      allocP  += _nTheta;
    }

    return;

  /*** BIVARIATE random effects with univariately modelled  margins (copula???) ***/
  /*** ======================================================================== ***/
  case 2:
    throw returnR("Error in RandomCL.cpp: RandomCL::updateRE2. Not implemented for _nTheta = 2", 1);
    return;

  /*** MULTI(>2)VARIATE random effects ***/
  /*** =============================== ***/
  default:
    throw returnR("Error in RandomCL.cpp: RandomCL::updateRE2. Not implemented for _nTheta > 2", 1);
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
RandomCL::updateRE2Bi(int *accept,                       double *ll, 
                      MatrixRect<double> *prob,          MatrixRect<double> *Propprob,
                      const BiGspline2 *Gspl,
                      const MatrixRect<double> *offset,  const MatrixRect<int> *y,
                      const MatrixRect<int> *alloc)
{
  static int i, j, length0, i0, i1;
  static double llCl;
  const int *yP, *niP, *allocP;
  const double *xP, *vP, *offsetP, *REMeanP, *inv_sigma2_d2P, *d_knotsP0, *d_knotsP1;
  const MatrixLT<double> *xxP, *vvP;
  const MatrixRect<double> *xvP;
  int *acceptP;
  double *probP, *thetaP, *etaXP, *etaVP, *etaP, *meanP, *REWMeanP;

  if (_nTheta != BiGspline2A::_dim){
    REprintf("_nTheta=%d,  BiGspline2A::_dim=%d\n", _nTheta, BiGspline2A::_dim);
    throw returnR("Error in RandomCL.cpp: updateRE2Bi(). Not implemented for this dimension.", 1);
  }

  *ll     = 0;
  acceptP = accept;
  yP      = y->aconst();
  probP   = prob->a();
  thetaP  = _Theta.a();
  etaXP   = _etaX.a();
  etaVP   = _etaV.a();
  etaP    = _eta.a();
  offsetP = offset->aconst();
  xP      = _x.aconst();
  vP      = _v.aconst();
  xxP     = _xx;
  vvP     = _vv;
  xvP     = _xv;
  niP     = _ni.aconst(); 
  allocP  = alloc->aconst();

  REMeanP  = _REMean.aconst();              /*** constant over cycles of the loop over observations ***/
  meanP    = _ThetaBar.a();                 /*** working space at each cycle of the loop            ***/
  REWMeanP = _REWMean.a();                  /*** working space at each cycle of the loop            ***/

  length0        = Gspl->lengthconst()[0];
  d_knotsP0      = Gspl->d_knotsconst()[0].aconst();
  d_knotsP1      = Gspl->d_knotsconst()[1].aconst();
  inv_sigma2_d2P = Gspl->inv_sigma2_d2const();

  for (i = 0; i < _N; i++){  /** loop over clusters            **/
    
    /*** Mean of the random effect given allocation, store it in _ThetaBar. ***/
    i0 = *allocP % length0;
    i1 = *allocP / length0;
    meanP[0] = REMeanP[0] + d_knotsP0[i0];
    meanP[1] = REMeanP[1] + d_knotsP1[i1];

    /*** Update _REWMean = (d*sigma)^{-2}*(intcpt + d*knot[alloc]).  ***/
    REWMeanP[0] = inv_sigma2_d2P[0] * meanP[0];
    REWMeanP[1] = inv_sigma2_d2P[1] * meanP[1];

    /*** Sample new value of the random effect ***/
    mcmc_common::update_reg_gamermanCL(acceptP, &llCl, &_U, &_I, 
                                       etaXP, _PropetaX.a(), etaVP, _PropetaV.a(), etaP, _Propeta.a(), probP, Propprob->a(), _work.a(),
                                       offsetP, thetaP, _PropTheta.a(), yP, xP, vP, xxP, vvP, xvP, 
                                       *niP, _C, _p, _q, 
                                       meanP, REWMeanP, inv_sigma2_d2P, NULL, true, &_PropMean,
                                       ll_cumlogitFS2, ll_cumlogitNR2, "RandomCL::updateRE2Bi");
    *ll += llCl; 

    /*** Increase pointers ***/
    acceptP++;
    yP      += (*niP);
    j       = _C * (*niP);
    probP   += j + (*niP);
    thetaP  += _nTheta;
    etaXP   += (*niP);
    etaVP   += j;
    etaP    += j;
    offsetP += j;
    xP      += _p * (*niP);
    vP      += _q * (*niP);
    xxP     += (*niP);
    vvP     += (*niP);
    xvP     += (*niP);
    niP++;
    allocP++;
  }                          /** end of the loop over clusters **/

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
RandomCL::invert_REInvVar()
{
  int Rank;

  _REVar.array2mat(_REInvVar.aconst());
  Rank = _REVar.chol_inv(0, NULL);

  return Rank;
}

void
RandomCL::invert_REInvVarFromL()
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
RandomCL::print(const int &colmax) const
{
  int cmax = (colmax < 0 ? _N : (colmax <= _N ? colmax : _N));

  Rprintf("\nObject of class RandomCL:\n");
  Rprintf("=========================\n");
  Rprintf("q=%d,  C=%d,  p=%d,  nTheta=%d\n", _q, _C, _p, _nTheta);
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
