/*** BiGspline2.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  21/03/2007
//                            
//                                       BASICS:  ??/??/2007
//
// PURPOSE: Class to hold parameters related to the bivariate G-spline
//
//   * Methods for class BiGspline2 are implemented in BiGspline2_Util.cpp
//                                                     BiGspline2_updateLambda.cpp
//                                                     BiGspline2_updateWeights.cpp
//
/* ********************************************************************************* */

#include "BiGspline2.h"

/* ********************************************************************************* */
/* Constructors and destructors                                                      */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */

/***** Nonparametric constructor *****/
BiGspline2::BiGspline2()
{
  int j;
  _total_length = 0;
  _max_length   = 0;

  _a          = MatrixRect<double>();
  _expa       = MatrixRect<double>();
  _sumexpa    = 0;
  _a_max      = 0;
  _log_null_w = 0;

  _k_effect     = 0;
  _ind_w_effect = MatrixRect<int>();

  _neighbor_system  = 0;
  _order            = 0;
  _type_update_a    = 0;
  _prior_for_lambda = 0;
  _equal_lambda     = 1;

  _workD = MatrixRect<double>();
  _workI = MatrixRect<int>();

  for (j = 0; j < BiGspline2A::_dim; j++){
    _K[j]      = 0;
    _length[j] = 0;

    _knots[j]     = MatrixRect<double>();
    _knots2[j]    = MatrixRect<double>();
    _invsigma2[j] = 0;
    _sigma[j]     = 0;

    _d_knots[j]       = MatrixRect<double>();
    _sigma_d[j]       = 0;
    _inv_sigma2_d2[j] = 0;

    _sumexpa_margin[j] = MatrixRect<double>();

    _aReference[j] = 0;

    _lambda[j]  = 0;
    _penalty[j] = 0;
    _lambdaPriorShape[j]    = 0;
    _lambdaPriorInvScale[j] = 0;

    _gMean[j] = 0;
  }

  for (j = 0; j < BiGspline2A::_LTdim; j++) _gVar[j]  = 0;

  _abscis = MatrixRect<double>();
  for (j = 0; j < mcmc_Gspline2::_liwv; j++) _iwv[j] = 0;
  for (j = 0; j < mcmc_Gspline2::_lrwv; j++) _rwv[j] = 0;
  for (j = 0; j < mcmc_Gspline2::_nabscis; j++){
    _hx[j]  = 0;
    _hpx[j] = 0;
  }
}


/***** Parametric constructor *****/
//
//  dimPar[_dim+1]:    
//                     dimPar[1] = _dim                                             // only used to check whether it is equal to 2
//                     dimPar[+_dim] = _K[0,...,_dim-1]
//
//  penMix_iPar[2+_dim]:
//                     penMix_iPar[1]     = _type_update_a
//                     penMix_iPar[+1]    = _equal_lambda
//                     penMix_iPar[+1]    = _prior_for_lambda
//                     penMix_iPar[+_dim] = order of the penalty for each margin   // ignored if neighbor system is not uniCAR
//                                          only the first value is used, the rest is ignored (all margins must have the same _order)
//                     penMix_iPar[+1]    = neighbor system
//
//  penMix_dPar[_dim+2*_sum_length+_dim]:    
//                     penMix_dPar[_dim]                   = basis standard deviations
//                     penMix_dPar[+_length[0]+_length[1]] = knots for each margin
//
//  aPar[1+_dim]:
//              aPar[1] = type of identifiability constraint (0=_Mean_, 1=_Reference_, 2=_General_)   // only checked whether equal to _Reference_
//              aPar[+1] = indeces of reference a coefficients
//                         BOTH INPUT and STORED : indeces from -K,...,K
//                         some time ago I have STORED: indeces from 0,...,2*K
//
//  lambdaPrior[2*_dim]:
//                     if (_prior_for_lambda == Gamma_) lambdaPrior[_dim] = prior shape parameters for gamma prior
//                                                      lambdaPrior[+_dim] = prior inverse scale (rate) parameters for gamma prior
//                     if (_prior_for_lambda == SDUnif_) lambdaPrior[_dim] = upper limits for the uniform distributions
//
//  lambda_a[_dim+_sum_length]:
//                     lambda_a[_dim]           = lambda for each margin
//                     lambda_a[+_total_length] = a coefficients (matrix in column major order)
//           
//  invVar[LT(dim)]:   lower triangle of the inverse of the overall covariance matrix D,
//                     it is assumed that it is diagonal, only diagonal elements are used to compute _d_knots, _sigma_d, _inv_sigma2_d2
//
BiGspline2::BiGspline2(const int *dimPar,  
                       const int *penMix_iPar,     const double *penMix_dPar,  
                       const int *aPar,
                       const double *lambdaPrior,  const double *lambda_a,
                       const double *invVar)
{
  int j, k0, k1;  
  bool a_inf;
  const double *cdP, *cdP2;
  const int *ciP;
  double *dP, *dP0, *dP1;
  int *iP;
  double dd;

  if (dimPar[0] != BiGspline2A::_dim){
    REprintf("dimPar[0]=%d. Dimension should be equal to %d\n", dimPar[0], BiGspline2A::_dim);
    throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2() constructor", 1);
  }

  /*** dimPar --> _K, _length, _max_length, _total_length ***/
  /*** ================================================== ***/
  _max_length   = 0;
  _total_length = 1;
  for (j = 0; j < BiGspline2A::_dim; j++){
    _K[j]       = dimPar[1+j];
    _length[j]  = 2*_K[j] + 1;
    if (_length[j] <= 0){
      REprintf("_K[%d]=%d leads to non-positive _length[%d]=%d\n", j, _K[j], j, _length[j]);
      throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2() constructor", 1);
    }
    if (_length[j] > _max_length) _max_length = _length[j];
    _total_length *= _length[j];
  }


  /*** _workD, _workI ***/
  /*** ============== ***/
  _workD = MatrixRect<double>(1, _total_length);
  _workI = MatrixRect<int>(1, _total_length);


  /*** _log_null_w ***/
  /*** =========== ***/
  _log_null_w = log(mcmc_Gspline2::_null_mass) - log(1 - mcmc_Gspline2::_null_mass) - log(double(_total_length));
  

  /*** penMix_dPar, invVar --> _knots, _knots2, _sigma, _invsigma2, _d_knots, _sigma_d, _inv_sigma2_d2 ***/
  /*** =============================================================================================== ***/  
  cdP  = penMix_dPar + BiGspline2A::_dim;     /** beginning of knots **/
  cdP2 = invVar;
  for (j = 0; j < BiGspline2A::_dim; j++){
    _sigma[j] = penMix_dPar[j];
    if (_sigma[j] <= 0){
      REprintf("_sigma[%d]=%g is not positive\n", j, _sigma[j]);
      throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2() constructor", 1);
    }
    _invsigma2[j] = 1/(_sigma[j]*_sigma[j]);

    _knots[j]   = MatrixRect<double>(1, _length[j], cdP);
    _knots2[j]  = MatrixRect<double>(1, _length[j], cdP);
    _d_knots[j] = MatrixRect<double>(1, _length[j], cdP);
  
    _knots2[j].square();

    dd  = sqrt(*cdP2);                     /** inverse of the overall standard deviation tau for j-th margin **/
    _sigma_d[j]       = _sigma[j]/dd;
    _inv_sigma2_d2[j] = 1/(_sigma_d[j]*_sigma_d[j]); 
    _d_knots[j].divide(dd);

    cdP  += _length[j];
    cdP2 += (BiGspline2A::_dim - j);
  }

  /*** penMix_iPar --> _type_update_a, _equal_lambda, _prior_for_lambda, _order, _neighbor_system ***/
  /*** ========================================================================================== ***/
  /** _neighbor_system **/
  ciP = penMix_iPar + 1 + 1 + 1 + BiGspline2A::_dim;
  switch (*ciP){
  case BiGspline2A::uniCAR:
  case BiGspline2A::eight_neighbors:
  case BiGspline2A::twelve_neighbors:
    _neighbor_system = *ciP;
    break;
  default:
    REprintf("_neighbor_system=%d\n", *ciP);
    throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...), incorrect _neighbor_system argument", 1);    
  }

  ciP = penMix_iPar;

  /** _type_update_a **/
  switch (*ciP){
  case mcmc_Gspline2::Slice:
  case mcmc_Gspline2::ARS_quantile:
  case mcmc_Gspline2::ARS_mode:
    _type_update_a = *ciP;
    break;
  case mcmc_Gspline2::Block:
    REprintf("Block update of a coefficients not implemented for multivariate G-splines\n");
  default:
    REprintf("_type_update_a=%d\n", *ciP);
    throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...), incorrect _type_update_a argument", 1);
  }
  ciP++;

  /** _equal_lambda **/
  _equal_lambda = (*ciP != 0) ? true : false;
  if (_neighbor_system == BiGspline2A::eight_neighbors || _neighbor_system == BiGspline2A::twelve_neighbors) _equal_lambda = true;
  ciP++;

  /** _prior_for_lambda **/
  switch (*ciP){
  case GMRF_Gspline_Util::_Fixed_:
  case GMRF_Gspline_Util::_Gamma_:
  case GMRF_Gspline_Util::_SDUnif_:
    _prior_for_lambda = *ciP;
    break;
  default:
    REprintf("_prior_for_lambda=%d\n", *ciP);
    throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...), unknown _prior_for_lambda", 1);
  }
  ciP++;

  /** _order **/
  switch (_neighbor_system){
  case BiGspline2A::uniCAR:
    if (*ciP < 0){
      REprintf("_order[%d]=%d\n", 0, *ciP);
      throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...), order may not be negative", 1);
    }
    _order = *ciP;
    ciP++;
    for (j = 1; j < BiGspline2A::_dim; j++){
      if (*ciP < 0){
        REprintf("_order[%d]=%d\n", j, *ciP);
        throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...), order may not be negative", 1);
      }
      if (_order != *ciP){
        REprintf("_order[0]=%d,  _order[1]=%d\n", _order, *ciP);
        throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...). Both _order values must be the same", 1);
      }
      ciP++;
    }
    break;
  case BiGspline2A::eight_neighbors:
    _order = 2;                                  /** 2 here to update 'a' with as small as possible autocorrelation **/
    for (j = 0; j < BiGspline2A::_dim; j++) ciP++;
    break;
  case BiGspline2A::twelve_neighbors:
    _order = 3;                                 /** 3 here to update 'a' with as small as possible autocorrelation **/
    for (j = 0; j < BiGspline2A::_dim; j++) ciP++;
    break;
  default:
    REprintf("_neighbor_system=%d\n", *ciP);
    throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...), incorrect _neighbor_system argument", 1);
  }


  /*** aPar --> _aReference ***/
  /*** ==================== ***/
  /** Check '_aIdent' **/
  ciP = aPar;
  switch (*ciP){
  case GMRF_Gspline_Util::_Mean_:
    throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...). _aIdent = _Mean_ not implemented", 1);    
  case GMRF_Gspline_Util::_Reference_:
    break;
  case GMRF_Gspline_Util::_General_:
    throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...). _aIdent = _General_ not implemented", 1);
  default:
    REprintf("_aIdent = %d\n", *ciP);
    throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...). Incorrect _aIdent argument", 1);    
  }
  ciP++;

  /** _aReference **/
  for (j = 0; j < BiGspline2A::_dim; j++){
    if (*ciP < -_K[j] || *ciP > _K[j]){
      REprintf("_K[%d]=%d,  _aReference[%d]=%d\n", j, _K[j], j, *ciP);
      throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...). Incorrect _aReference argument", 1);
      _aReference[j] = *ciP;
      ciP++;
    }
  }


  /*** lambdaPrior --> _lambdaPriorShape, _lambdaPriorInvScale ***/
  /*** ======================================================= ***/
  switch (_prior_for_lambda){
  case GMRF_Gspline_Util::_Fixed_:
    for (j = 0; j < BiGspline2A::_dim; j++) _lambdaPriorShape[j] = _lambdaPriorInvScale[j] = 0;    
    break;
  case GMRF_Gspline_Util::_Gamma_:
    cdP  = lambdaPrior;
    cdP2 = lambdaPrior + BiGspline2A::_dim;
    if (_equal_lambda){
      if (*cdP <= 0 || *cdP2 <= 0){
        REprintf("_lambdaPriorShape=%g,  _lambdaPriorInvScale=%g\n", *cdP, *cdP2);
        throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...). Incorrect parameters for the lambda prior", 1);
      }
      for (j = 0; j < BiGspline2A::_dim; j++){
        _lambdaPriorShape[j]    = *cdP;
        _lambdaPriorInvScale[j] = *cdP2;
      }
    }
    else{
      for (j = 0; j < BiGspline2A::_dim; j++){
        if (*cdP <= 0 || *cdP2 <= 0){
          REprintf("_lambdaPriorShape[%d]=%g,  _lambdaPriorInvScale[%d]=%g\n", j, *cdP, j, *cdP2);
          throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...). Incorrect parameters for the lambda prior", 1);
        }
        _lambdaPriorShape[j]    = *cdP;
        _lambdaPriorInvScale[j] = *cdP2;
        cdP++;
        cdP2++;
      }
    }
    break;
  case GMRF_Gspline_Util::_SDUnif_:
    cdP  = lambdaPrior;
    if (_equal_lambda){
      if (*cdP <= 0){
        REprintf("Upper limit for _SDUnif_ lambda prior=%g\n", *cdP);
        throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...). Incorrect parameters for the lambda prior", 1);
      }
      for (j = 0; j < BiGspline2A::_dim; j++){
        _lambdaPriorShape[j]    = 1/((*cdP)*(*cdP));
        _lambdaPriorInvScale[j] = 0;
      }
    }
    else{
      for (j = 0; j < BiGspline2A::_dim; j++){
        if (*cdP <= 0){
          REprintf("%d-th upper limit for _SDUnif_ lambda prior=%g\n", j, *cdP);
          throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...). Incorrect parameters for the lambda prior", 1);
        }
        _lambdaPriorShape[j]    = 1/((*cdP)*(*cdP));
        _lambdaPriorInvScale[j] = 0;
        cdP++;
      }
    }
    break;
  default:
    REprintf("_prior_for_lambda=%d\n", *ciP);
    throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...), unknown _prior_for_lambda", 1);
  }


  /*** lambda_a --> _lambda, _a, _expa, _sumexpa, _a_max, _sumexpa_margin, _k_effect, _ind_w_effect ***/  
  /*** ============================================================================================ ***/
  /** _lambda **/
  cdP = lambda_a;
  if (_equal_lambda){
    if (*cdP <= 0){
      REprintf("_lambda=%g\n", *cdP);
      throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...). Incorrect initial value of _lambda", 1);
    }
    _lambda[0] = *cdP;    
    cdP++;
    for (j = 1; j < BiGspline2A::_dim; j++){
      _lambda[j] = _lambda[0];
      cdP++;
    }
  }
  else{
    for (j = 0; j < BiGspline2A::_dim; j++){
      if (*cdP <= 0){
        REprintf("_lambda[%d]=%g\n", j, *cdP);
        throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...). Incorrect initial value of _lambda", 1);
      }
      _lambda[j] = *cdP;
      cdP++;
    }
  }

  /** _a, _expa, _sumexpa, _sumexpa_margin, _a_max       **/
  /** _k_effect, _ind_w_effect (initially all indeces)   **/
  if (BiGspline2A::_dim != 2){
    REprintf("BiGspline2A::_dim=%d\n", BiGspline2A::_dim);
    throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...). The following is implemented only for _dim=2", 1);
  }

  for (j = 0; j < BiGspline2A::_dim; j++) _sumexpa_margin[j] = MatrixRect<double>(1, _length[j]);
  _a            = MatrixRect<double>(_length[0], _length[1], cdP);
  _expa         = MatrixRect<double>(_length[0], _length[1]);  
  _ind_w_effect = MatrixRect<int>(_length[0], _length[1]);

  cdP      = _a.aconst();
  dP       = _expa.a();
  dP1      = _sumexpa_margin[1].a();
  _a_max   = *cdP;
  _sumexpa = 0;
  a_inf    = false;
  
  iP        = _ind_w_effect.a();
  _k_effect = 0;  

  for (k1 = -_K[1]; k1 <= _K[1]; k1++){
    dP0 = _sumexpa_margin[0].a();
    for (k0 = -_K[0]; k0 <= _K[0]; k0++){
      if (*cdP > _AK_EMAX){
        REprintf("a[%d,%d]=%g\n", k0, k1, *cdP);
        a_inf = true;
      }
      if (*cdP > _a_max) _a_max = *cdP;
      *dP = exp_AK(*cdP);                    // _expa[j] = exp(a[j]);
      _sumexpa += *dP;
      *dP0     += *dP;
      *dP1     += *dP;
      cdP++;
      dP++;
      dP0++;

      *iP = 1;
      iP++;
      _k_effect++;
    }
    dP1++;
  }
  if (a_inf){
    _sumexpa = R_PosInf;
    REprintf("_sumexpa = Infinity\n");
    throw returnR("Error in BiGspline2.cpp: BiGspline2::BiGspline2(...). Too high values of 'a' coefficients", 1);
  }

  
  /*** _penalty ***/
  /*** ======== ***/
  this->penalty();


  /*** _gMean, _gVar ***/
  /*** ============= ***/
  this->gMoments();
  

  /*** ARS stuff (initialize _abscis, _iwv, _rwv, _hx, _hpx) ***/
  /*** ===================================================== ***/
  for (j = 0; j < mcmc_Gspline2::_liwv; j++) _iwv[j] = 0;
  for (j = 0; j < mcmc_Gspline2::_lrwv; j++) _rwv[j] = 0;
  for (j = 0; j < mcmc_Gspline2::_nabscis; j++) _hx[j] = _hpx[j] = 0;
  _abscis = MatrixRect<double>(mcmc_Gspline2::_nabscis, _total_length);
  this->find_start_abscisBi();

}


/***** Copy constructor *****/
BiGspline2::BiGspline2(const BiGspline2 &Gs)
{
  throw returnR("CONTACT THE AUTHOR: Copy constructor for class BiGspline2 has not been written yet!", 99);
}


/***** Assignment operator *****/
BiGspline2&
BiGspline2::operator=(const BiGspline2 &Gs)
{
  int j;
  _total_length = Gs._total_length;
  _max_length   = Gs._max_length;

  _a          = Gs._a;
  _expa       = Gs._expa;
  _sumexpa    = Gs._sumexpa;
  _a_max      = Gs._a_max;
  _log_null_w = Gs._log_null_w;

  _k_effect     = Gs._k_effect;
  _ind_w_effect = Gs._ind_w_effect;

  _neighbor_system  = Gs._neighbor_system;
  _order            = Gs._order;
  _type_update_a    = Gs._type_update_a;
  _prior_for_lambda = Gs._prior_for_lambda;
  _equal_lambda     = Gs._equal_lambda;

  _workD = Gs._workD;
  _workI = Gs._workI;

  for (j = 0; j < BiGspline2A::_dim; j++){
    _K[j]      = Gs._K[j];
    _length[j] = Gs._length[j];

    _knots[j]     = Gs._knots[j];
    _knots2[j]    = Gs._knots2[j];
    _invsigma2[j] = Gs._invsigma2[j];
    _sigma[j]     = Gs._sigma[j];

    _d_knots[j]       = Gs._d_knots[j];
    _sigma_d[j]       = Gs._sigma_d[j];
    _inv_sigma2_d2[j] = Gs._inv_sigma2_d2[j];

    _sumexpa_margin[j] = Gs._sumexpa_margin[j];

    _aReference[j] = Gs._aReference[j];

    _lambda[j]  = Gs._lambda[j];
    _penalty[j] = Gs._penalty[j];
    _lambdaPriorShape[j]    = Gs._lambdaPriorShape[j];
    _lambdaPriorInvScale[j] = Gs._lambdaPriorInvScale[j];

    _gMean[j] = Gs._gMean[j];
  }

  for (j = 0; j < BiGspline2A::_LTdim; j++) _gVar[j]  = Gs._gVar[j];

  _abscis = Gs._abscis;
  for (j = 0; j < mcmc_Gspline2::_liwv; j++) _iwv[j] = Gs._iwv[j];
  for (j = 0; j < mcmc_Gspline2::_lrwv; j++) _rwv[j] = Gs._rwv[j];
  for (j = 0; j < mcmc_Gspline2::_nabscis; j++){
    _hx[j]  = Gs._hx[j];
    _hpx[j] = Gs._hpx[j];
  }

  return *this;
}


/***** Destructor *****/
BiGspline2::~BiGspline2()
{
}
