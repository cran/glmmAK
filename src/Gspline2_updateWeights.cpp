/*** Gspline2_updateWeights.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  15/01/2007
//                       (many parts are only slightly adjusted code from Gspline_update_a.cpp from bayesSurv package)
//                            
//                       update_a1:  15/01/2007
//             update_alla_lambda1:  15/01/2007
//                    full_a_pars1:  15/01/2007
//              find_start_abscis1:  16/01/2007
//
// PURPOSE: Update of the G-spline weights
//
/* ********************************************************************************* */

#include "Gspline2.h"


/***** ============================================================================ *****/
/***** update_a1:  Update one particular 'a' coefficient of the UNIVARIATE G-spline *****/
/*****                                                                              *****/
/***** ============================================================================ *****/
//
// aa:                           Pointer to _a that is updated
// expaa:                        Pointer to _expa corresponding to _a that is updated
// Abscis[mcmc_Gspline2::_nabscis]:  Starting abscisae for ARS or working space for the slice sampler
// ia:                           Index of the a that is updated (on the scale -K,...,K)
// a_ipars[2]:                   a_ipars[0] = number of all observations
//                               a_ipars[1] = number of observations currently belonging to the component of the a
// overrelax:                    1/0 indicating whether overrelaxation is to be used (used only when the slice sampler is used,
//                               ignored otherwise)
//
void
Gspline2::update_a1(double *aa,     double *expaa,       double *Abscis,
                    const int *ia,  const int *a_ipars,  const int *overrelax)
{
  static double a_pars[4];
  static double newa;
  static double *sumexpa, *hx, *hpx, *abscis;
  static int i;
  static int _ONE_INT = 1;

  sumexpa = _sumexpa.a();

  this->full_a_pars1(a_pars + 0, a_pars + 1, ia, aa);          /* compute mean and inv. variance of [a[ia] | a[-ia], lambda] */
  a_pars[2] = *expaa;
  a_pars[3] = *sumexpa;

  /*** Find the mode of the full conditional (if necessary)                                              ***/
  /*** Store this mode in _abscis[ia][1]                                                                 ***/
  /*** Compute either starting abscissae for ARS or initial guesses for interval defining the slice      ***/
  switch (_type_update_a){
  case mcmc_Gspline2::Slice:
  case mcmc_Gspline2::ARS_mode:
    mcmc_Gspline2::find_eval_abscis(Abscis, _hx.a(), _hpx.a(), ia, &_ONE_INT, aa, a_pars, a_ipars);
    break;

  case mcmc_Gspline2::ARS_quantile:    /** Starting abscissae are taken as quantiles of an upper hull from the previous iteration **/
                                       /** Evaluate the function to sample from in starting abscissae                             **/
    abscis = Abscis;
    hx     = _hx.a();
    hpx    = _hpx.a();
    for (i = 0; i < mcmc_Gspline2::_nabscis; i++){
      mcmc_Gspline2::full_a_logdens(abscis, hx, hpx, a_pars, a_ipars);            
      abscis++;
      hx++;
      hpx++;
    }
    break;
  default: 
    throw returnR("Error in Gspline2_updateWeights.cpp: Gspline2::update_a1. Unimplemented _type_update_a", 1);
  }

  /** Check whether starting abscissae/initial guesses for the interval defining the slice lie on correct size of the mode   **/
  mcmc_Gspline2::check_abscis(Abscis, _hx.a(), _hpx.a(), a_pars, a_ipars);

  /** Sample new a **/
  switch (_type_update_a){  
  case mcmc_Gspline2::Slice:
    mcmc_Gspline2::sample_a_by_slice(&newa, Abscis, _hx.a(), _hpx.a(), ia, &_ONE_INT, aa, a_pars, a_ipars, overrelax);
    break;

  case mcmc_Gspline2::ARS_quantile:
  case mcmc_Gspline2::ARS_mode:
    mcmc_Gspline2::sample_a_by_ARS(&newa, Abscis, _hx.a(), _hpx.a(), _rwv.a(), _iwv.a(), ia, &_ONE_INT, aa, a_pars, a_ipars, &_type_update_a);
    break;    

  default:
    throw returnR("Error in Gspline2_updateWeights.cpp: Gspline2::update_a1. Unimplemented _type_update_a", 1);
  }

  /** Update exp(a) and sum(exp(a)) **/
  *aa       = newa;
  *sumexpa -= *expaa;
  if (*aa >= mcmc_Gspline2::_log_inf){
    *aa      = mcmc_Gspline2::_log_inf;
    *expaa   = mcmc_Gspline2::_exp_emax;
    *sumexpa = mcmc_Gspline2::_exp_emax;
  }
  else{
    *expaa    = exp(*aa);
    *sumexpa += *expaa;
  }
 
  return;
}


/***** ============================================================================ *****/
/***** update_alla_lambda1                                                          *****/
/*****                                                                              *****/
/***** ============================================================================ *****/
//
// Update all G-spline coefficients and also smoothing hyperparameter(s) lambda
//  * function for UNIVARIATE G-spline
//
//  Sequence of updating (if not in 1 block): try to minimize autocorrelation
//    * e.g. _dim = 1, _order = 3, update a[0], a[4], a[8], ...
//                                        a[1], a[5], a[9], ...
//                                        a[2], a[6], a[10], ...
//                                        a[3], a[7], a[11], ...
//
// a_ipars[2]:          a_ipars[0] = number of all observations
//                      a_ipars[1] = working place to store current mixtureNM
// mixtureNM:           numbers of observations belonging to each component (long vector)  [_total_length]
// overrelax:           1/0 indicating whether overrelaxation is to be used (used only when the slice sampler is used,
//                      ignored otherwise)
//
void
Gspline2::update_alla_lambda1(int *a_ipars,  const int *mixtureNM,  const int *overrelax)
{
  static int accept;
  static const int *order;
  static int order_1;
  static int K, k, ia, imax;
  static double *aa, *expaa, *Abscis, *aaS, *expaaS, *AbscisS, *sumexpa;
  static const int *mNM, *mNMS;
  static double amax;
  static double tmp, aRef;

  order   = _order.aconst();
  order_1 = *order + 1;
  K       = _K.aconst()[0];

  if (*order == 0) return;                    /*** a's are fixed, and there is no lambda ***/

  if (_type_update_a < mcmc_Gspline2::Block){     /*** Slice, ARS_quantile, ARS_mode ***/

    /*** Update lambda ***/
    /*** ------------- ***/
    this->updateLambda();

    /*** Update a coefficients ***/
    /*** --------------------- ***/
    aaS     = _a->a();
    expaaS  = _expa->a();
    AbscisS = _abscis->a();
    mNMS    = mixtureNM;
    for (k = -K; k <= -K+(*order); k++){
      aa     = aaS;
      expaa  = expaaS;
      Abscis = AbscisS;
      mNM   = mNMS;
      for (ia = k; ia <= K; ia += order_1){
        a_ipars[1] = *mNM;
        this->update_a1(aa, expaa, Abscis, &ia, a_ipars, overrelax);

        aa     += order_1;
        expaa  += order_1; 
        Abscis += order_1*mcmc_Gspline2::_nabscis;
        mNM    += order_1;
      }
      aaS++;
      expaaS++;
      AbscisS += mcmc_Gspline2::_nabscis;
      mNMS++;
    }

    /*** Adjust a's to satisfy constraints ***/
    /*** --------------------------------- ***/
    switch (_aIdent){
    case GMRF_Gspline_Util::_Mean_:
      /* Compute the mean of a's */
      aa  = _a->a();
      tmp = *aa;
      for (k = -K+1; k <= K; k++){
        aa++;
        tmp += *aa;
      }
      tmp /= _sum_length;

      /* Center a and compute exp(a), sum(exp(a)) */
      aa       = _a->a();
      expaa    = _expa->a();
      sumexpa  = _sumexpa.a();
      *sumexpa = 0.0;
      for (k = -K; k <= K; k++){
        *aa      -= tmp;
        *expaa    = exp(*aa);
        *sumexpa += *expaa;
        aa++;
        expaa++;
      }
      break;

    case GMRF_Gspline_Util::_Reference_:
      /* Find maximal a */
      imax = -K;
      aa   = _a->a();
      amax = *aa;
      for (k = -K+1; k <= K; k++){
        aa++;
        if (*aa >= amax){
          amax = *aa;
          imax = k;
        }
      }

      /* Adjust a's in the case they are too high and find the reference a */
      if (amax > mcmc_Gspline2::_a_ceil){
        aa   = _a->a();
        _aReference.a()[0] = imax;
        for (k = -K; k <= K; k++){
          *aa -= amax;
          aa++;
        }
        aRef = 0.0;
      }
      else{
        aa  = _a->a();
        for (k = -K; k < _aReference.a()[0]; k++) aa++;
        aRef = *aa;
      }

      /* Compute exp(a), sum(exp(a)) */
      aa       = _a->a();
      expaa    = _expa->a();
      sumexpa  = _sumexpa.a();
      *sumexpa = 0.0;
      for (k = -K; k <= K; k++){
        *aa      -= aRef;
        *expaa    = exp(*aa);
        *sumexpa += *expaa;
        aa++;
        expaa++;
      }
      break;

    case GMRF_Gspline_Util::_General_:
      throw returnR("Error in Gspline2_updateWeights.cpp: Gspline2::update_alla_lambda1(). _aIdent=_General_ not implemented", 1);

    default:
      throw returnR("Error in Gspline2_updateWeights.cpp: Gspline2::update_alla_lambda1(). Unimplemented _aIdent", 1);      
    }

    /*** Compute the value of the penalty ***/
    /*** -------------------------------- ***/
    this->penalty(_penalty.a(), _a->aconst(), _order.aconst(), _length.aconst());
  }
  else{                                /* _type_update_a == Block */
    //Rprintf("\n_a (in margin 0): ");
    //AK_BLAS_LAPACK::printArray(_a->aconst(), _length.aconst()[0]);
    //Rprintf("_expa (sum(exp(a)=%g), in margin 0): ", _sumexpa.aconst()[0]);
    //AK_BLAS_LAPACK::printArray(_expa->aconst(), _length.aconst()[0]);
    //Rprintf("_w (_minw = %g): ", _minw);
    //AK_BLAS_LAPACK::printArray(_w.aconst(), _sum_length);
    //Rprintf("_Da: ");
    //AK_BLAS_LAPACK::printArray(_Da.aconst(), _sum_length);
    //Rprintf("_Qa (penalty=%g): ", _penalty.aconst()[0]);
    //AK_BLAS_LAPACK::printArray(_Qa.aconst(), _sum_length);

    /*** Update a coefficients ***/
    /*** --------------------- ***/
    GMRF_Gspline::update(&accept, _a->a(), _lambda.a(), _expa->a(), _sumexpa.a(), _w.a(), &_minw, 
                         _Da.a(), _Qa.a(), _penalty.a(), _workML.a(), _worka.a(), _workGMRF.a(),
                         mixtureNM, &_prior_for_lambda, _lambdaPriorPar.aconst(), _par_rscale, 
                         _Q.aconst(), _order.aconst(), _diffOper.aconst(), 
                         &GMRF_Gspline::_null_mass, &_aIdent, _aReference.aconst(), &_sum_length, a_ipars, &_AK_ZERO_INT);
    if (accept){
      if (_aIdent == GMRF_Gspline_Util::_Reference_){
        /*** Find maximal a ***/
        /*** -------------- ***/
        imax = -K;
        aa   = _a->a();
        amax = *aa;
        for (k = -K+1; k <= K; k++){
          aa++;
          if (*aa >= amax){
            amax = *aa;
            imax = k;
          }
        }

        /*** Adjust a's in the case they are too high and find the reference a ***/
        /*** ----------------------------------------------------------------- ***/
        if (amax > mcmc_Gspline2::_a_ceil){
          aa   = _a->a();
          _aReference.a()[0] = imax;
          for (k = -K; k <= K; k++){
            *aa -= amax;
            aa++;
          } 
          GMRF_Gspline_Util::update4_ll12(_expa->a(), _sumexpa.a(), _Da.a(), _penalty.a(), _Qa.a(), 
                                          _w.a(), &_minw, _a->a(), _order.aconst(), _diffOper.aconst(), &_sum_length);
        }
      }
    }
  }

  return;
}  /** end of the function update_alla_lambda1 **/


/***** ============================================================================================= *****/
/***** full_a_pars1:  Parameters of the full conditional distribution of 1 a in the UNIVARIATE case  *****/
/*****                                                                                               *****/
/***** ============================================================================================= *****/
//
// mean:       Computed mean of the full conditional distribution
// invvar:     Computed inverse variance of the full conditional distribution
//
// ia:         Index of the a that is updated (on the scale -K,...,K)
// aa:         Pointer to _a that is updated
//
void
Gspline2::full_a_pars1(double* mean, double* invvar,  const int *ia, const double *aa) const
{
  int K = _K.aconst()[0];
  if (*ia < -K || *ia > K){
    REprintf("K=%d, ia=%d\n", K, *ia);
    throw returnR("Error in Gspline2_updateWeights: Gspline2::full_a_pars1(). Argument ia out of the range", 1);
  }

  switch (_order.aconst()[0]){
  case 1:
    if (*ia > -K && *ia < K){
      *mean = (aa[-1] + aa[1])/2;
      *invvar = 2*_lambda.aconst()[0];
    }      
    else{     // ia = -K or K
      if (*ia == -K) *mean = aa[1];
      else           *mean = aa[-1];
      *invvar = _lambda.aconst()[0];
    }
    return;

  case 2:
    if (*ia >= -K+2 && *ia <= K-2){
      *mean = (-aa[-2] + 4*aa[-1] + 4*aa[1] - aa[2])/6;
      *invvar = 6*_lambda.aconst()[0];
    }
    else{
      if (*ia == -K+1 || *ia == K-1){
        if (*ia == -K+1) *mean = (2*aa[-1] + 4*aa[1] - aa[2])/5;
        else             *mean = (-aa[-2] + 4*aa[-1] + 2*aa[1])/5;
        *invvar = 5*_lambda.aconst()[0];
      }
      else{     // ia = -K or K
        if (*ia == -K) *mean = 2*aa[1] - aa[2];
        else           *mean = -aa[-2] + 2*aa[-1];
        *invvar = _lambda.aconst()[0];
      }
    }
    return;

  case 3:
    if (*ia >= -K+3 && *ia <= K-3){
      *mean = (aa[-3] - 6*aa[-2] + 15*aa[-1] + 15*aa[1] - 6*aa[2] + aa[3])/20;
      *invvar = 20*_lambda.aconst()[0];
    }
    else{
      if (*ia == -K+2 || *ia == K-2){
        if (*ia == -K+2) *mean = (-3*aa[-2] + 12*aa[-1] + 15*aa[1] - 6*aa[2] + aa[3])/19;
        else             *mean = (aa[-3] - 6*aa[-2] + 15*aa[-1] + 12*aa[1] - 3*aa[2])/19;
        *invvar = 19*_lambda.aconst()[0];
      }
      else{
        if (*ia == -K+1 || *ia == K-1){
          if (*ia == -K+1) *mean = (3*aa[-1] + 12*aa[1] - 6*aa[2] + aa[3])/10; 
          else             *mean = (aa[-3] - 6*aa[-2] + 12*aa[-1] + 3*aa[1])/10;
          *invvar = 10*_lambda.aconst()[0];
        }
        else{     // ia = -K or K
          if (*ia == -K) *mean = 3*aa[1] - 3*aa[2] + aa[3];
          else           *mean = aa[-3] - 3*aa[-2] + 3*aa[-1];
          *invvar = _lambda.aconst()[0];
        }  
      }  
    }  
    return;

  default:
    REprintf("_order=%d\n", _order.aconst()[0]);
    throw returnR("Error in Gspline2_updateWeights: Gspline2::full_a_pars1(). Unimplemented _order.", 1);
  }
}


/***** ============================================================================================================ *****/
/***** find_start_abscis1:  Find starting abscissae for a particular a in the UNIVARIATE case                       *****/
/*****                                                                                                              *****/
/***** ============================================================================================================ *****/
//
//  * use mean and mean +- 3*sd of the distribution [a[ia] | a[-ia], lambda]
//    as starting abscissae
//  * do not check whether they lie on both sides of the mode,
//    this will be done always before sampling is done
//
void
Gspline2::find_start_abscis1()
{
  if (_dim != 1){
    throw returnR("Error in Gspline2_updateWeights.cpp: Gspline2::find_start_abscis1. Implemented only for UNIVARIATE G-splines", 1);
  }

  if (mcmc_Gspline2::_nabscis != 3){
    throw returnR("Dear Arnost, please update Gspline2::find_start_abscis1() function after changing _nabscis ;-)", 1);
  }

  double mean, invvar, three_sd;
  int K            = _K.aconst()[0];  
  const double *aa = _a->aconst();
  double *Abscis   = _abscis->a();
  for (int ia=-K; ia <= K; ia++){
    this->full_a_pars1(&mean, &invvar, &ia, aa);
    three_sd = 3/sqrt(invvar);
    Abscis[0] = mean - three_sd;
    Abscis[1] = mean;
    Abscis[2] = mean + three_sd;
    
    aa++;
    Abscis += mcmc_Gspline2::_nabscis;
  }

  return;
}
