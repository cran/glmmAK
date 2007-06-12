/*** Gspline2.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  17/10/2006
//                            
//                                       BASICS:  18/10/2006
//
// PURPOSE: Class to hold parameters related to the G-spline
//          with weights that are composed as a product of marginal ones
//
//   * Methods for class Gspline2 are implemented in Gspline2_Util.cpp
//                                                   Gspline2_updateLambda.cpp
//                                                   Gspline2_updateWeights.cpp
//
/* ********************************************************************************* */

#include "Gspline2.h"

/* ********************************************************************************* */
/* Constructors and destructors                                                      */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */

/***** Nonparametric constructor *****/
Gspline2::Gspline2()
  : _dim(0), _sum_length(0), _max_length(0), _knots(NULL), _knots2(NULL), _d_knots(NULL),
    _a(NULL), _expa(NULL), _aIdent(0), _aContrast(NULL), _type_update_a(0), _abscis(NULL), _minw(0),
    _equal_lambda(true), _prior_for_lambda(0)
{
  for (int i=0; i < 6; i++) _par_rscale[i] = 0.0;
}


/***** Parametric constructor *****/
//
//  dimPar[_dim+1]:    
//                     dimPar[1] = _dim
//                     dimPar[+_dim] = _K[0,...,_dim-1]
//
//  penMix_iPar[2+_dim]:
//                     penMix_iPar[1]     = _type_update_a                          // added on 16/01/2007
//                     penMix_iPar[+1]    = _equal_lambda
//                     penMix_iPar[+1]    = _prior_for_lambda
//                     penMix_iPar[+_dim] = order of the penalty for each margin
//                     penMix_iPar[+_dim] = neighbor system (IGNORED HERE)
//
//  penMix_dPar[_dim+2*_sum_length+_dim]:    
//                     penMix_dPar[_dim]         = basis standard deviations
//                     penMix_dPar[+_sum_length] = knots for each margin
//
//  aPar[1+_dim]:
//              aPar[1] = type of identifiability constraint (0=_Mean_, 1=_Reference_, 2=_General_)
//              aPar[+1] = indeces of reference a coefficients
//                         BOTH INPUT and STORED : indeces from -K,...,K
//                         some time ago I have STORED: indeces from 0,...,2*K
//
//  aContrast[_sum_length]:
//              identifiability constraints for each margin
//              typically: _aIdent = _Mean_     => _aContrast[i] = (1/_length[i],...,1/_length[i])
//                         _aIdent = _Reference => _aContrast[i] = (0,...,0,1,0,...,0)
//
//  lambdaPrior[2*_dim]:
//                     if (_prior_for_lambda == Gamma_) lambdaPrior[_dim] = prior shape parameters for gamma prior
//                                                      lambdaPrior[+_dim] = prior inverse scale (rate) parameters for gamma prior
//                     if (_prior_for_lambda == SDUnif_) lambdaPrior[_dim] = upper limits for the uniform distributions
//
//  lambda_a[_dim+_sum_length]:
//                     lambda_a[_dim] = lambda for each margin
//                     lambda_a[+_sum_length] = a coefficients for each margin 
//           
//  invVar[LT(dim)]:   lower triangle of the overall covariance matrix D,
//                     it is assumed that it is diagonal, only diagonal elements are used to compute _d_knots, _sigma_d, _inv_sigma2_d2
//
Gspline2::Gspline2(const int *dimPar,  
                   const int *penMix_iPar,     const double *penMix_dPar,  
                   const int *aPar,            const double *aContrast,
                   const double *lambdaPrior,  const double *lambda_a,
                   const double *invVar)
{
  int i, j;
  double dd, dd2;
  int *iP;
  double *dP, *dP2, *dP3, *dP4, *dP5;
  MatrixRect<double> *MP;
  const MatrixRect<double> *cMP;
  const int *ciP, *ciP2; 
  const double *cdP, *cdP2, *cdP3;  

  int nc, nworkGMRF2, nworksGMRF2[5];
  double FF;

  _dim = dimPar[0];
  if (_dim < 0) throw returnR("Error in Gspline2.cpp: Gspline2::Gspline2(...), _dim < 0", 1);
  if (_dim == 0){
    _sum_length       = 0;
    _max_length       = 0;
    _knots            = NULL;
    _knots2           = NULL;
    _d_knots          = NULL;
    _a                = NULL;
    _expa             = NULL;
    _aIdent           = 0;
    _aContrast        = NULL;
    _type_update_a    = 0;
    _abscis           = NULL;
    _minw             = 0;
    _equal_lambda     = true;
    _prior_for_lambda = 0; 
    for (i=0; i < 6; i++) _par_rscale[i] = 0.0;
  }
  else{    
    /*** _K, _length, _sum_length, _max_length, _workD, _workI ***/
    /*** ===================================================== ***/
    _K = MatrixRect<int>(1, _dim, dimPar+1);
    _length = MatrixRect<int>(1, _dim);
    iP = _length.a();
    ciP = _K.aconst();
    for (i = 0; i < _dim; i++){
      *iP = 2*(*ciP) + 1;
      iP++;
      ciP++;      
    }
    _sum_length = _length.sum();
    _max_length = _length.max();
    _workD = MatrixRect<double>(1, _max_length);
    _workI = MatrixRect<int>(1, _max_length);

    /*** _log_null_w ***/
    /*** =========== ***/
    _log_null_w = MatrixRect<double>(1, _dim);
    dP = _log_null_w.a();
    ciP = _length.aconst();
    for (i = 0; i < _dim; i++){
      *dP = log(mcmc_Gspline2::_null_mass) - log(1 - mcmc_Gspline2::_null_mass) - log(double(*ciP));
      dP++;
      ciP++;
    }

    /*** _equal_lambda, _prior_for_lambda, _order ***/
    /*** ======================================== ***/
    _equal_lambda = (penMix_iPar[1] == 0 ? false : true);
    switch (penMix_iPar[2]){
    case GMRF_Gspline_Util::_Fixed_:
    case GMRF_Gspline_Util::_Gamma_:
    case GMRF_Gspline_Util::_SDUnif_:
      _prior_for_lambda = penMix_iPar[2];
      break;
    default:
      throw returnR("Error in Gspline2.cpp: Gspline2::Gspline2(...), unknown _prior_for_lambda", 1);
    }

    _order = MatrixRect<int>(1, _dim, penMix_iPar + 3);

    /*** _sigma, _invsigma2, _knots, _knots2, _a, _expa, _sumexpa, _penalty, _lambda ***/
    /*** =========================================================================== ***/
    _sigma     = MatrixRect<double>(1, _dim, penMix_dPar);
    _invsigma2 = MatrixRect<double>(1, _dim);
    _knots     = new MatrixRect<double>[_dim];
    _knots2    = new MatrixRect<double>[_dim];
    _a         = new MatrixRect<double>[_dim];
    _expa      = new MatrixRect<double>[_dim];
    if (!_knots || !_knots2 || !_a || !_expa) throw returnR("Out of memory in Gspline2.cpp: Gspline2::Gspline2(...)", 99);       
    _sumexpa = MatrixRect<double>(1, _dim);   
    _penalty = MatrixRect<double>(1, _dim);
        
    cdP = penMix_dPar;                    // *sigma
    cdP2 = cdP + _dim;                    // *knots
    _lambda = MatrixRect<double>(1, _dim, lambda_a);
    cdP3 = lambda_a + _dim;               // *a

    if (_equal_lambda){
      dP = _lambda.a();
      dd = *dP;
      for (i = 1; i < _dim; i++){
        dP++;
        *dP = dd;
      }
    }

    ciP = _length.aconst();
    ciP2 = _order.aconst();
    dP = _invsigma2.a();
    dP2 = _sumexpa.a();
    dP4 = _penalty.a();
    for (i = 0; i < _dim; i++){
      *dP = 1/((*cdP) * (*cdP));         // _invsigma2[i] = 1/(_sigma[i]*_sigma[i]);
      dP++;                              // _invsigma2++;
      cdP++;                             // _sigma++;

      _knots[i] = MatrixRect<double>(1, *ciP, cdP2);
      _knots2[i] = MatrixRect<double>(1, *ciP);
      dP5 = _knots2[i].a();

      _a[i] = MatrixRect<double>(1, *ciP, cdP3);

      this->penalty(dP4, cdP3, ciP2, ciP); 

      _expa[i] = MatrixRect<double>(1, *ciP);
      dP3 = _expa[i].a();
      *dP2 = 0;                       // _sumexpa[i] = 0;   
      for (j = 0; j < *ciP; j++){
        *dP3 = exp_AK(*cdP3);         // compute _expa[i];
        *dP2 += *dP3;                 // _sumexpa[i] += exp(_a[i]);
        dP3++;                        // _expa[i]++;
        *dP5 = (*cdP2) * (*cdP2);     // _knots2[i][j] = _knots[i][j]^2;
        dP5++;                        // _knots2[i]++;
        cdP2++;                       // knots++;
        cdP3++;                       // _a++;
      }
      dP2++;                          // _sumexpa++;
      dP4++;                          // _penalty++;
      ciP++;                          // _length++;
      ciP2++;                         // _order++;
    }

    /*** _gMean, _gVar ***/
    /*** ============= ***/
    _gMean = MatrixRect<double>(1, _dim);
    _gVar  = MatrixRect<double>(1, _dim);
    this->gMoments();


    /*** _d_knots, _sigma_d, _inv_sigma2_d2 ***/
    /*** ================================== ***/
    _d_knots = new MatrixRect<double>[_dim];
    if (!_d_knots) throw returnR("Out of memory in Gspline2.cpp: Gspline2::Gspline2(...)", 99);
    _sigma_d = MatrixRect<double>(1, _dim);
    _inv_sigma2_d2 = MatrixRect<double>(1, _dim);

    cdP = invVar;
    cdP2 = _sigma.aconst();
    ciP = _length.aconst();
    cMP = _knots;

    dP2 = _sigma_d.a();
    dP3 = _inv_sigma2_d2.a();
    MP = _d_knots;
    for (j = _dim; j > 0; j--){
      dd = sqrt(*cdP);               // dd = sqrt(invVar[j,j]);
      *dP2 = (*cdP2)/dd;             // _sigma_d[j] = _sigma[j]/dd;
      *dP3 = 1/((*dP2) * (*dP2));    // _inv_sigma2_d2[j] = 1/_sigma_d[j]^2;  
           
      *MP = MatrixRect<double>(1, *ciP);
      dP4 = MP->a();               // _d_knots[j];
      cdP3 = cMP->aconst();        // _knots[j];
      for (i = 0; i < *ciP; i++){
        *dP4 = (*cdP3)/dd;         // _d_knots[j][i] = _knots[j][i]/dd;
        dP4++;                     // _d_knots[j]++;
        cdP3++;                    // _knots[j]++;
      }

      cdP += j;                     // _invVar += j (to the next diagonal element);
      cdP2++;                       // _sigma++;
      ciP++;                        // _length++;
      cMP++;                        // _knots++;

      dP2++;                        // _sigma_d++;
      dP3++;                        // _inv_sigma2_d2++;
      MP++;                         // _d_knots++;
    }

    /*** _lambdaPriorShape, _lambdaPriorInvScale, _lambdaPriorPar ***/
    /*** ======================================================== ***/
    _lambdaPriorPar = MatrixRect<double>(2, _dim);

    switch (_prior_for_lambda){
    case GMRF_Gspline_Util::_Fixed_:
      _lambdaPriorShape = MatrixRect<double>(1, _dim);
      _lambdaPriorInvScale = MatrixRect<double>(1, _dim);
      break;

    case GMRF_Gspline_Util::_Gamma_:
      _lambdaPriorShape = MatrixRect<double>(1, _dim, lambdaPrior);
      _lambdaPriorInvScale = MatrixRect<double>(1, _dim, lambdaPrior+_dim);
      if (_equal_lambda){
        dP  = _lambdaPriorShape.a();
        dP2 = _lambdaPriorInvScale.a();
        dd  = *dP;
        dd2 = *dP2;
        for (i = 1; i < _dim; i++){
          dP++;
          dP2++;
          *dP = dd;
          *dP2 = dd2;
        }
      }
      cdP  = _lambdaPriorShape.aconst();
      cdP2 = _lambdaPriorInvScale.aconst();
      dP   = _lambdaPriorPar.a();
      for (i = 0; i < _dim; i++){
        *dP = *cdP;
        dP++;
        cdP++;
        *dP = *cdP2;
        dP++;
        cdP2++;
      }
      break;

    case GMRF_Gspline_Util::_SDUnif_:
      _lambdaPriorShape = MatrixRect<double>(1, _dim);
      _lambdaPriorInvScale = MatrixRect<double>(1, _dim);
      if (_equal_lambda){
        dP = _lambdaPriorShape.a();
        dd = lambdaPrior[0];
        for (i = 0; i < _dim; i++){
          *dP = 1/(dd*dd);
          dP++;
        }
      }
      else{
        dP = _lambdaPriorShape.a();
        cdP = lambdaPrior;
        for (i = 0; i < _dim; i++){
          *dP = 1/((*cdP)*(*cdP));
          dP++;
          cdP++;
        }
      }
      cdP  = _lambdaPriorShape.aconst();
      dP   = _lambdaPriorPar.a() + 1;
      for (i = 0; i < _dim; i++){
        *dP = *cdP;
        dP += 2;
        cdP++;
      }
      break;
    }

    /*** _aIdent ***/
    /*** ======= ***/
    switch (*aPar){
    case GMRF_Gspline_Util::_Mean_:
    case GMRF_Gspline_Util::_Reference_:
    case GMRF_Gspline_Util::_General_:
      _aIdent = *aPar;
      break;
    default:
      REprintf("_aIdent = %d\n", *aPar);
      throw returnR("Error in Gspline2.cpp: Gspline2::Gspline2(...), incorrect _aIdent argument", 1);
    }

    /*** _aReference (do not shift it from scale -K,...,K to scale 0,...,2*K), _aContrast ***/
    /*** ================================================================================ ***/
    _aContrast = new MatrixRect<double>[_dim];
    if (!_aContrast) throw returnR("Out of memory in Gspline2.cpp: Gspline2::Gspline2(...)", 99);       

    _aReference = MatrixRect<int>(1, _dim, aPar+1);
    iP = _aReference.a();
    ciP = _length.aconst();
    ciP2 = _K.aconst();
    cdP = aContrast;
    for (i = 0; i < _dim; i++){
      _aContrast[i] = MatrixRect<double>(1, *ciP, cdP);
      cdP += *ciP;

      /***  *iP += *ciP2;                 _aReference[i] += K[i]                                                             ***/
      /***  if (*iP >= *ciP) throw returnR("Error in Gspline2.cpp: Gspline2::Gspline2(...), _aReference out of range", 1);   ***/
      if (*iP < -(*ciP2) || *iP > *ciP2){
        REprintf("_dim=%d,  margin=%d, _aReference=%d,  _K=%d\n", _dim, i, *iP, *ciP2);
        throw returnR("Error in Gspline2.cpp: Gspline2::Gspline2(...), _aReference out of range", 1);
      }
      iP++;
      ciP++;
      ciP2++;    
    }

    /*** _type_update_a ***/ 
    /*** ============== ***/
    switch (*penMix_iPar){
    case mcmc_Gspline2::Slice:
    case mcmc_Gspline2::ARS_quantile:
    case mcmc_Gspline2::ARS_mode:
    case mcmc_Gspline2::Block:
      _type_update_a = *penMix_iPar;
      break;
    default:
      REprintf("_type_update_a = %d\n", *penMix_iPar);
      throw returnR("Error in Gspline2.cpp: Gspline2::Gspline2(...), incorrect _type_update_a argument", 1);
    }

    /*** _abscis, _iwv, _rwv, _hx, _hpx ***/    
    /*** ============================== ***/
    _abscis = new MatrixRect<double>[_dim];
    if (!_abscis) throw returnR("Out of memory in Gspline2.cpp: Gspline2::Gspline2(...)", 99);
    for (i = 0; i < _dim; i++){
      _abscis[i] = MatrixRect<double>(mcmc_Gspline2::_nabscis, _length.aconst()[i]);      
    }
    if (_dim == 1) this->find_start_abscis1();
    else{
      REprintf("_dim=%d\n", _dim);
      throw returnR("Error in Gspline2.cpp: Gspline2::Gspline2(...). Starting abscissae for ARS not implemented for _dim != 1", 1);
    }

    _iwv  = MatrixRect<int>(1, mcmc_Gspline2::_liwv);
    _rwv  = MatrixRect<double>(1, mcmc_Gspline2::_lrwv);
    _hx   = MatrixRect<double>(1, mcmc_Gspline2::_nabscis);
    _hpx  = MatrixRect<double>(1, mcmc_Gspline2::_nabscis);
   
    /*** _w, _minw, _Q, _Da, _Qa, _diffOper, _par_rscale, _LTna, _LTna_1, _nworkML, _nworka, _nworkGMRF, _workML, _worka, _workGMRF ***/
    /*** ========================================================================================================================== ***/
    if (_dim == 1){
      _LTna = (_sum_length * (1 + _sum_length))/2;
      _LTna_1 = ((_sum_length - 1) * _sum_length)/2;

      _nworkML = _sum_length + 1 + _sum_length + 1 + 2*_sum_length + 1 + (_sum_length-1) + _LTna_1 + _sum_length + _LTna;
      _nworka  = _sum_length + _LTna_1 + (_sum_length-1) + 2*_sum_length + 1 + _sum_length + 1 + 2*_sum_length + 1;
      nc = 0;
      nworksGMRF2[0] = nc*nc;                                            /** log_density_Ax_x **/
      nworksGMRF2[1] = (nc*(1+nc))/2 + _sum_length*nc + nc;              /** rGMRF_inputArgs  **/
      nworksGMRF2[2] = (_sum_length > nc ? _sum_length : nc);            /** rGMRF            **/
      nworksGMRF2[3] = nc;                                               /** dGMRF_inputArgs  **/ 
      nworksGMRF2[4] = _sum_length;                                      /** dGMRF            **/
      nworkGMRF2     = maxArray(nworksGMRF2, 5);
      _nworkGMRF     = 1 + 4 + _sum_length*nc + (nc*(1+nc))/2 + nworkGMRF2;

      _w        = MatrixRect<double>(1, _sum_length);
      _Q        = MatrixLT<double>(_sum_length);
      _Da       = MatrixRect<double>(1, _sum_length);
      _Qa       = MatrixRect<double>(1, _sum_length);
      _diffOper = MatrixRect<int>(1, _order.aconst()[0]+1);

      _workML   = MatrixRect<double>(1, _nworkML);
      _worka    = MatrixRect<double>(1, _nworka);
      _workGMRF = MatrixRect<double>(1, _nworkGMRF);

      if (_type_update_a == mcmc_Gspline2::Block){
        GMRF::diff_operator(_diffOper.a(), _order.aconst());
        GMRF::Q_matrix(_Q.a(), _order.aconst(), &_sum_length);      
        GMRF_Gspline_Util::update4_ll12(_expa->a(), _sumexpa.a(), _Da.a(), _penalty.a(), _Qa.a(), 
                                        _w.a(), &_minw, _a->aconst(), _order.aconst(), _diffOper.aconst(), &_sum_length);

        FF = 1.1;
        GMRF::dscale_norm_const(&FF, _par_rscale);
      }
    }
    else{
      REprintf("_dim=%d\n", _dim);
      throw returnR("Error in Gspline2.cpp: Gspline2::Gspline2(...). Stuff for block update of a not implemented for _dim != 1", 1);
    }
  }
}


/***** Copy constructor *****/
Gspline2::Gspline2(const Gspline2 &Gs)
{
  int i;

  _dim = Gs._dim;
  if (_dim > 0){
    _knots     = new MatrixRect<double>[_dim];
    _knots2    = new MatrixRect<double>[_dim];
    _d_knots   = new MatrixRect<double>[_dim];
    _a         = new MatrixRect<double>[_dim];
    _expa      = new MatrixRect<double>[_dim];
    _aContrast = new MatrixRect<double>[_dim];
    _abscis    = new MatrixRect<double>[_dim];
    if (!_knots || !_knots2 || !_d_knots || !_a || !_expa || !_aContrast || !_abscis) 
      throw returnR("Out of memory in Gspline2.cpp: Gspline2::Gspline2(const Gspline2&)", 99);   
  }
  else{
    _knots     = NULL;
    _knots2    = NULL;
    _d_knots   = NULL;
    _a         = NULL;
    _expa      = NULL;
    _aContrast = NULL;
    _abscis    = NULL;   
  }

  _gMean = Gs._gMean;
  _gVar  = Gs._gVar;

  _K          = Gs._K;
  _length     = Gs._length;
  _sum_length = Gs._sum_length;  
  _max_length = Gs._max_length;  
  _workD      = Gs._workD;
  _workI      = Gs._workI;

  _log_null_w = Gs._log_null_w;

  _aIdent        = Gs._aIdent;
  _aReference    = Gs._aReference;
  _type_update_a = Gs._type_update_a;
  for (i = 0; i < _dim; i++){
    _knots[i]     = Gs._knots[i];
    _knots2[i]    = Gs._knots2[i];
    _d_knots[i]   = Gs._d_knots[i];
    _a[i]         = Gs._a[i];
    _expa[i]      = Gs._expa[i];
    _aContrast[i] = Gs._aContrast[i];
    _abscis[i]    = Gs._abscis[i];
  }

  _invsigma2     = Gs._invsigma2;
  _sigma         = Gs._sigma;
  _sigma_d       = Gs._sigma_d;
  _inv_sigma2_d2 = Gs._inv_sigma2_d2;
  _sumexpa       = Gs._sumexpa;

  _lambda  = Gs._lambda;
  _penalty = Gs._penalty;

  _order               = Gs._order;
  _equal_lambda        = Gs._equal_lambda;  
  _prior_for_lambda    = Gs._prior_for_lambda;
  _lambdaPriorShape    = Gs._lambdaPriorShape;
  _lambdaPriorInvScale = Gs._lambdaPriorInvScale;
  _lambdaPriorPar      = Gs._lambdaPriorPar;  

  _iwv  = Gs._iwv;
  _rwv  = Gs._rwv;  
  _hx   = Gs._hx;
  _hpx  = Gs._hpx;

  _w        = Gs._w;
  _minw     = Gs._minw;
  _Q        = Gs._Q;  
  _Da       = Gs._Da;  
  _Qa       = Gs._Qa;  
  _diffOper = Gs._diffOper;  

  for (i = 0; i < 6; i++) _par_rscale[i] = Gs._par_rscale[i];

  _LTna      = Gs._LTna;  
  _LTna_1    = Gs._LTna_1;  
  _nworkML   = Gs._nworkML;  
  _nworka    = Gs._nworka;  
  _nworkGMRF = Gs._nworkGMRF;  
  _workML    = Gs._workML;  
  _worka     = Gs._worka;  
  _workGMRF  = Gs._workGMRF;  
}


/***** Assignment operator *****/
Gspline2&
Gspline2::operator=(const Gspline2 &Gs)
{
  int i;

  if (_dim > 0){
    delete [] _knots;
    delete [] _knots2;
    delete [] _d_knots;
    delete [] _a;
    delete [] _expa;
    delete [] _aContrast;
  }

  _dim = Gs._dim;
  if (_dim > 0){
    _knots     = new MatrixRect<double>[_dim];
    _knots2    = new MatrixRect<double>[_dim];
    _d_knots   = new MatrixRect<double>[_dim];
    _a         = new MatrixRect<double>[_dim];
    _expa      = new MatrixRect<double>[_dim];
    _aContrast = new MatrixRect<double>[_dim];
    _abscis    = new MatrixRect<double>[_dim];
    if (!_knots || !_knots2 || !_a || !_expa || !_aContrast || !_abscis) 
      throw returnR("Out of memory in Gspline2.cpp: Gspline2::operator=(const Gspline2&)", 99);   
  }
  else{
    _knots     = NULL;
    _knots2    = NULL;
    _d_knots   = NULL;
    _a         = NULL;
    _expa      = NULL;
    _aContrast = NULL;
    _abscis    = NULL;   
  }

  _gMean = Gs._gMean;
  _gVar  = Gs._gVar;

  _K          = Gs._K;
  _length     = Gs._length;
  _sum_length = Gs._sum_length;
  _max_length = Gs._max_length;  
  _workD      = Gs._workD;
  _workI      = Gs._workI;

  _log_null_w = Gs._log_null_w;

  _aIdent        = Gs._aIdent;
  _aReference    = Gs._aReference;
  _type_update_a = Gs._type_update_a;
  for (i = 0; i < _dim; i++){
    _knots[i]     = Gs._knots[i];
    _knots2[i]    = Gs._knots2[i];
    _d_knots[i]   = Gs._d_knots[i];
    _a[i]         = Gs._a[i];
    _expa[i]      = Gs._expa[i];
    _aContrast[i] = Gs._aContrast[i];
    _abscis[i]    = Gs._abscis[i];
  }

  _invsigma2     = Gs._invsigma2;
  _sigma         = Gs._sigma;
  _sigma_d       = Gs._sigma_d;
  _inv_sigma2_d2 = Gs._inv_sigma2_d2;
  _sumexpa       = Gs._sumexpa;

  _lambda  = Gs._lambda;
  _penalty = Gs._penalty;

  _order               = Gs._order;
  _equal_lambda        = Gs._equal_lambda;  
  _prior_for_lambda    = Gs._prior_for_lambda;
  _lambdaPriorShape    = Gs._lambdaPriorShape;
  _lambdaPriorInvScale = Gs._lambdaPriorInvScale;
  _lambdaPriorPar      = Gs._lambdaPriorPar;  

  _iwv  = Gs._iwv;
  _rwv  = Gs._rwv;  
  _hx   = Gs._hx;
  _hpx  = Gs._hpx;

  _w        = Gs._w;
  _minw     = Gs._minw;
  _Q        = Gs._Q;  
  _Da       = Gs._Da;  
  _Qa       = Gs._Qa;  
  _diffOper = Gs._diffOper;  

  for (i = 0; i < 6; i++) _par_rscale[i] = Gs._par_rscale[i];

  _LTna      = Gs._LTna;  
  _LTna_1    = Gs._LTna_1;  
  _nworkML   = Gs._nworkML;  
  _nworka    = Gs._nworka;  
  _nworkGMRF = Gs._nworkGMRF;  
  _workML    = Gs._workML;  
  _worka     = Gs._worka;  
  _workGMRF  = Gs._workGMRF;  

  return *this;
}


/***** Destructor *****/
Gspline2::~Gspline2()
{
  if (_dim > 0){
    delete [] _knots;
    delete [] _knots2;
    delete [] _d_knots;
    delete [] _a;
    delete [] _expa;
    delete [] _aContrast;
    delete [] _abscis;
  }
}

