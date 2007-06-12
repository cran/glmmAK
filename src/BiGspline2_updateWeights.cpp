/*** BiGspline2_updateWeights.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  23/03/2007
//                       (many parts are only slightly adjusted code from Gspline_update_a.cpp from bayesSurv package
//                        and from Gspline2_updateWeights.cpp of this package)
//                            
//                      update_aBi:  27/03/2007
//            update_alla_lambdaBi:  29/03/2007
//                   full_a_parsBi:  23/03/2007  
//             find_start_abscisBi:  27/03/2007
//
//
// PURPOSE: Update of the G-spline weights
//
/* ********************************************************************************* */

#include "BiGspline2.h"


/***** ============================================================================= *****/
/***** update_aBi:  Update one particular 'a' coefficient of the BIVARIATE G-spline  *****/
/*****                                                                               *****/
/***** * It also updates _expa, _sumexpa_margin, _sumexpa                            *****/
/***** * It does not update _a_max!                                                  *****/
/*****                                                                               *****/
/***** ============================================================================= *****/
//
// aa:                               Pointer to _a that is updated
// a_log_inf:                        OUTPUT: if <> 0 then the sampled a is higher than _log_inf
//                                           and all a coefficients must be shifted down to avoid 
//                                           numerical problems in the update of subsequent a coefficients
//
// expaa:                            Pointer to _expa corresponding to _a that is updated
//                                   (updated only if a_log_inf = 0)
// sexpaam0:                         Pointer to _sumexpa_margin[0] corresponding to _a that is updated
//                                   (updated only if a_log_inf = 0)
// sexpaam1:                         Pointer to _sumexpa_margin[1] corresponding to _a that is updated
//                                   (updated only if a_log_inf = 0)
//
// Abscis[mcmc_Gspline2::_nabscis]:  Starting abscisae for ARS or working space for the slice sampler
//
// ia[2]:                            Double index of the a that is updated (on the scale -K,...,K)
// a_ipars[2]:                       a_ipars[0] = number of all observations
//                                   a_ipars[1] = number of observations currently belonging to the component of the a
// overrelax:                        1/0 indicating whether overrelaxation is to be used (used only when the slice sampler is used,
//                                   ignored otherwise)
//
void
BiGspline2::update_aBi(double *aa,     int *a_log_inf,      double *expaa,       double *sexpaam0,      double *sexpaam1,   double *Abscis,
                       const int *ia,  const int *a_ipars,  const int *overrelax)
{
  static double a_pars[4];
  static double newa;
  static double *hx, *hpx, *abscis;
  static int i;
  static int _TWO_INT = 2;

  this->full_a_parsBi(a_pars + 0, a_pars + 1, ia, aa);          /* compute mean and inv. variance of [a[ia] | a[-ia], lambda] */
  a_pars[2] = *expaa;
  a_pars[3] = _sumexpa;

  /*** Find the mode of the full conditional (if necessary)                                              ***/
  /*** Store this mode in _abscis[ia][1]                                                                 ***/
  /*** Compute either starting abscissae for ARS or initial guesses for interval defining the slice      ***/
  switch (_type_update_a){
  case mcmc_Gspline2::Slice:
  case mcmc_Gspline2::ARS_mode:
    mcmc_Gspline2::find_eval_abscis(Abscis, _hx, _hpx, ia, &_TWO_INT, aa, a_pars, a_ipars);
    break;

  case mcmc_Gspline2::ARS_quantile:    /** Starting abscissae are taken as quantiles of an upper hull from the previous iteration **/
                                       /** Evaluate the function to sample from in starting abscissae                             **/
    abscis = Abscis;
    hx     = _hx;
    hpx    = _hpx;
    for (i = 0; i < mcmc_Gspline2::_nabscis; i++){
      mcmc_Gspline2::full_a_logdens(abscis, hx, hpx, a_pars, a_ipars);            
      abscis++;
      hx++;
      hpx++;
    }
    break;
  default: 
    throw returnR("Error in BiGspline2_updateWeights.cpp: BiGspline2::update_aBi. Unimplemented _type_update_a", 1);
  }

  /** Check whether starting abscissae/initial guesses for the interval defining the slice lie on correct size of the mode   **/
  mcmc_Gspline2::check_abscis(Abscis, _hx, _hpx, a_pars, a_ipars);

  /** Sample new a **/
  switch (_type_update_a){  
  case mcmc_Gspline2::Slice:
    mcmc_Gspline2::sample_a_by_slice(&newa, Abscis, _hx, _hpx, ia, &_TWO_INT, aa, a_pars, a_ipars, overrelax);
    break;

  case mcmc_Gspline2::ARS_quantile:
  case mcmc_Gspline2::ARS_mode:
    mcmc_Gspline2::sample_a_by_ARS(&newa, Abscis, _hx, _hpx, _rwv, _iwv, ia, &_TWO_INT, aa, a_pars, a_ipars, &_type_update_a);
    break;    

  default:
    throw returnR("Error in BiGspline2_updateWeights.cpp: BiGspline2::update_aBi. Unimplemented _type_update_a", 1);
  }

  /** Update _expa, _sumexpa, _sumexpa_margin, _a_max **/
  //if (*aa >= mcmc_Gspline2::_a_ceil){ /** Alternative **/
  if (newa >= mcmc_Gspline2::_log_inf){ 
    //Rprintf("\nWARNING: a[%d, %d]: %g-->%g;  ipars=(%d, %d);  pars=(%g, %g, %g, %g)", ia[0], ia[1], *aa, newa, a_ipars[0], a_ipars[1], a_pars[0], a_pars[1], a_pars[2], a_pars[3]);
    *aa        = newa;
    *a_log_inf = 1;    
  }
  else{
    *aa        = newa;
    *a_log_inf = 0;

    _sumexpa  -= *expaa;
    *sexpaam0 -= *expaa;
    *sexpaam1 -= *expaa;

    *expaa     = exp(*aa);
    _sumexpa  += *expaa;
    *sexpaam0 += *expaa;
    *sexpaam1 += *expaa;
  }
 
  return;
}


/***** ============================================================================ *****/
/***** update_alla_lambdaBi                                                         *****/
/*****                                                                              *****/
/***** ============================================================================ *****/
//
// Update all G-spline coefficients and also smoothing hyperparameter(s) lambda
//  * function for BIVARIATE G-spline
//
//  Sequence of updating (if not in 1 block): try to minimize autocorrelation
//           _dim = 2, _order = 3, update a[0,0], a[4,0], a[8,0], ...
//                                        a[0,4], a[4,4], a[8,4], ...
//                                        a[0,8], a[4,8], a[8,8],...
//                                        .....
//                                        a[0,1], a[4,1], a[8,1], ...
//                                        a[0,5], a[4,5], a[8,5], ...
//                                        a[0,9], a[4,9], a[8,9],...
//                                        .....
//                                        a[0,2], a[4,2], a[8,2], ...
//                                        a[0,6], a[4,6], a[8,6], ...
//                                        a[0,10], a[4,10], a[8,10],...
//                                        .....
//                                        a[0,3], a[4,3], a[8,3], ...
//                                        a[0,7], a[4,7], a[8,7], ...
//                                        a[0,11], a[4,11], a[8,11],...
//                                        .....
//                                        a[1,0], a[5,0], a[9,0], ...
//                                        a[1,4], a[5,4], a[9,4], ...
//                                        a[1,8], a[5,8], a[9,8],...
//                                        .....
//
// a_ipars[2]:          a_ipars[0] = number of all observations
//                      a_ipars[1] = working place to store current mixtureNM
//
// mixtureNM:           numbers of observations belonging to each component 
// [_total_length]      (long vector corresponding to a matrix _length[0] x _length[1] in column major order)
//
// overrelax:           1/0 indicating whether overrelaxation is to be used (used only when the slice sampler is used,
//                      ignored otherwise)
//
void
BiGspline2::update_alla_lambdaBi(int *a_ipars,  const int *mixtureNM,  const int *overrelax)
{
  static int order_1, l0_order_1, nabs_order_1, l0_nabs, l0_nabs_order_1;
  static int a_log_inf;
  static int k0, k1, l0, l1;
  static double tmp;
  static double *aa, *aaS, *aaST, *aaSTT, *aaAdj;
  static double *expaa, *expaaS, *expaaST, *expaaSTT, *expaaAdj;
  static double *Abscis, *AbscisS, *AbscisST, *AbscisSTT;
  static double *sexpaam0, *sexpaam0STT, *sexpaam0Adj;
  static double *sexpaam1, *sexpaam1STT, *sexpaam1Adj;
  static const int *mNM, *mNMS, *mNMST, *mNMSTT;
  static int *iP;
  static int ia[2];

  order_1         = _order + 1;
  l0_order_1      = _length[0] * order_1;
  nabs_order_1    = mcmc_Gspline2::_nabscis * order_1;
  l0_nabs         = _length[0] * mcmc_Gspline2::_nabscis;
  l0_nabs_order_1 = _length[0] * nabs_order_1;

  if (_order == 0) return;                        /*** a's are fixed, and there is no lambda ***/

  if (_type_update_a < mcmc_Gspline2::Block){     /*** Slice, ARS_quantile, ARS_mode ***/

    /*** Update lambda ***/
    /*** ------------- ***/
    this->updateLambda();

    /*** Update a coefficients ***/
    /*** --------------------- ***/
    aaSTT     = _a.a();
    expaaSTT  = _expa.a();
    AbscisSTT = _abscis.a();
    mNMSTT    = mixtureNM;
    sexpaam0STT = _sumexpa_margin[0].a();
    sexpaam1STT = _sumexpa_margin[1].a();

    for (k0 = -_K[0]; k0 <= -_K[0] + _order; k0++){       /*** loop over a coefficients ***/
      aaST       = aaSTT;
      expaaST    = expaaSTT;
      AbscisST   = AbscisSTT;
      mNMST      = mNMSTT;

      for (k1 = -_K[1]; k1 <= -_K[1] + _order; k1++){
        aaS       = aaST;
        expaaS    = expaaST;
        AbscisS   = AbscisST;
        mNMS      = mNMST;

        sexpaam1 = sexpaam1STT + (k1 + _K[1]);        
        for (ia[1] = k1; ia[1] <= _K[1]; ia[1] += order_1){
          aa       = aaS;
          expaa    = expaaS;
          Abscis   = AbscisS;
          mNM      = mNMS;

          sexpaam0 = sexpaam0STT + (k0 + _K[0]);
          for (ia[0] = k0; ia[0] <= _K[0]; ia[0] += order_1){
            a_ipars[1] = *mNM;
            this->update_aBi(aa, &a_log_inf, expaa, sexpaam0, sexpaam1, Abscis, ia, a_ipars, overrelax);

            if (a_log_inf){   /** All a's must be shifted down, added on 06/04/2007 **/

              /* Reset _sumexpa_margin[0] */
              sexpaam0Adj = _sumexpa_margin[0].a();
              for (l0 = -_K[0]; l0 <= _K[0]; l0++){
              *sexpaam0Adj = 0;
              sexpaam0Adj++;
              }

              /* Decrease/increase all a's such that _a_max will be equal to 0           */      
              /* Compute _expa, _sumexpa_margin, _sumexpa                                */
              tmp = *aa;

              aaAdj        = _a.a();
              expaaAdj     = _expa.a();
              sexpaam1Adj  = _sumexpa_margin[1].a();
              _sumexpa     = 0; 
              for (l1 = -_K[1]; l1 <= _K[1]; l1++){
                *sexpaam1Adj = 0;                          /** Reset _sumexpa_margin[1] **/
                sexpaam0Adj  = _sumexpa_margin[0].a();
                for (l0 = -_K[0]; l0 <= _K[0]; l0++){
                  *aaAdj       -= tmp;
                  *expaaAdj     = exp_AK(*aaAdj);
                  _sumexpa     += *expaaAdj;
                  *sexpaam0Adj += *expaaAdj;
                  *sexpaam1Adj += *expaaAdj;

                  aaAdj++;
                  expaaAdj++;
                  sexpaam0Adj++;
                }
                sexpaam1Adj++;        
              }

              _a_max = 0;
            }

	    aa       += order_1;
            expaa    += order_1;
            Abscis   += nabs_order_1;
            mNM      += order_1;
            sexpaam0 += order_1;           

          }
          aaS      += l0_order_1;
          expaaS   += l0_order_1;
          AbscisS  += l0_nabs_order_1;
          mNMS     += l0_order_1;
          sexpaam1 += order_1;
        }
        aaST     += _length[0];
        expaaST  += _length[0];
        AbscisST += l0_nabs;
        mNMST    += _length[0];
      }
      aaSTT++;
      expaaSTT++;
      AbscisSTT += mcmc_Gspline2::_nabscis;
      mNMSTT++;
    }                                    /*** end of the loop over a coefficients ***/    

    /*** Determine _a_max and adjust a coefficients if _a_max is too high or too low ***/
    /*** --------------------------------------------------------------------------- ***/
    aa     = _a.a();
    _a_max = *aa;
    for (k0 = 0; k0 < _total_length; k0++){
      if (*aa > _a_max) _a_max = *aa;
      aa++;
    }

    /** Too high or too low _a_max                                                                            **/
    if (_a_max > mcmc_Gspline2::_a_ceil || _a_max < mcmc_Gspline2::_a_floor){

      /* Reset _sumexpa_margin[0] */
      sexpaam0 = _sumexpa_margin[0].a();
      for (k0 = -_K[0]; k0 <= _K[0]; k0++){
        *sexpaam0 = 0;
        sexpaam0++;
      }

      /* Decrease/increase all a's such that _a_max will be equal to 0           */      
      /* Compute _expa, _sumexpa_margin, _sumexpa                                */
      tmp = _a_max;

      aa        = _a.a();
      expaa     = _expa.a();
      sexpaam1  = _sumexpa_margin[1].a();
      _sumexpa  = 0; 
      for (k1 = -_K[1]; k1 <= _K[1]; k1++){
        *sexpaam1 = 0;                          /** Reset _sumexpa_margin[1] **/
        sexpaam0 = _sumexpa_margin[0].a();
        for (k0 = -_K[0]; k0 <= _K[0]; k0++){
          *aa -= tmp;
          *expaa = exp_AK(*aa);
          _sumexpa  += *expaa;
          *sexpaam0 += *expaa;
          *sexpaam1 += *expaa;

          aa++;
          expaa++;
          sexpaam0++;
        }
        sexpaam1++;        
      }

      _a_max -= tmp;
    }

    /*** Determine _k_effect, _ind_w_effect          ***/
    /*** ------------------------------------------- ***/
    aa        = _a.a();
    iP        = _ind_w_effect.a();       
    _k_effect = 0;

    for (k1 = -_K[1]; k1 <= _K[1]; k1++){
      for (k0 = -_K[0]; k0 <= _K[0]; k0++){
        if (*aa - _a_max > _log_null_w){
          *iP = 1;
          _k_effect++;
        }
        else{
          *iP = 0;
        }
        aa++;
        iP++;
      }
    }

    /*** Compute the value of the penalty ***/
    /*** -------------------------------- ***/
    this->penalty();
  }
  else{
    REprintf("_type_update_a=%d\n");
    throw returnR("Error in BiGspline2_updateWeights.cpp: update_alla_lambdaBi(). Not implemented _type_update_a", 1);
  }

  return;
}


/***** ============================================================================================== *****/
/***** full_a_parsBi:  Parameters of the full conditional distribution of 1 a in the BIVARIATE case   *****/
/*****                                                                                                *****/
/***** ============================================================================================== *****/
//
// mean:       Computed mean of the full conditional distribution
// invvar:     Computed inverse variance of the full conditional distribution
//
// ia[2]:      Double index of the a that is updated (on the scale -K,...,K)
// aa:         Pointer to _a that is updated
//
void
BiGspline2::full_a_parsBi(double* mean, double* invvar,  const int *ia, const double *aa) const
{
  static int nr, K0, K1;
  static double mean0_1, mean1_0, invvar0_1, invvar1_0;
  static const double *a1, *a2, *a3, *a4, *am1, *am2, *am3, *am4;
  nr = _length[0];

  if (ia[0] < -_K[0] || ia[0] > _K[0] || ia[1] < -_K[1] || ia[1] > _K[1]){
    REprintf("K=(%d, %d),  ia=(%d, %d)\n", _K[0], _K[1], ia[0], ia[1]);
    throw returnR("Error in BiGspline2_updateWeights: BiGspline2::full_a_parsBi(). Argument ia out of the range", 1);
  }

  switch (_neighbor_system){
  case BiGspline2A::uniCAR:
    switch (_order){
    case 1:
      /** first conditional (with fixed column) **/
      if (ia[0] > -_K[0] && ia[0] < _K[0]){
        mean0_1   = (aa[-1] + aa[1])/2;
        invvar0_1 = 2*_lambda[0];
      }      
      else{
        if (ia[0] == -_K[0]) mean0_1 = aa[1];
        else                 mean0_1 = aa[-1];
        invvar0_1 = _lambda[0];
      }
   
      /** second conditional (with fixed row) **/
      if (ia[1] > -_K[1] && ia[1] < _K[1]){
        mean1_0   = (aa[-nr] + aa[nr])/2;
        invvar1_0 = 2*_lambda[1];
      }      
      else{
        if (ia[1] == -_K[1]) mean1_0 = aa[nr];
        else                 mean1_0 = aa[-nr];
        invvar1_0 = _lambda[1];
      }
      break;
 
    case 2:
      /** first conditional (with fixed column) **/
      if (ia[0] >= -_K[0]+2 && ia[0] <= _K[0]-2){
        mean0_1   = (-aa[-2] + 4*aa[-1] + 4*aa[1] - aa[2])/6;
        invvar0_1 = 6*_lambda[0];
      }
      else{
        if (ia[0] == -_K[0]+1 || ia[0] == _K[0]-1){
          if (ia[0] == -_K[0]+1) mean0_1 = (2*aa[-1] + 4*aa[1] - aa[2])/5;
          else                   mean0_1 = (-aa[-2] + 4*aa[-1] + 2*aa[1])/5;
          invvar0_1 = 5*_lambda[0];
        }
        else{
          if (ia[0] == -_K[0]) mean0_1 = 2*aa[1] - aa[2];
          else                 mean0_1 = -aa[-2] + 2*aa[-1];
          invvar0_1 = _lambda[0];
        }
      }

      /** second conditional (with fixed row) **/
      if (ia[1] >= -_K[1]+2 && ia[1] <= _K[1]-2){
        a1  = aa + nr;
        a2  = a1 + nr;
        am1 = aa - nr;
        am2 = am1 - nr;
        mean1_0   = (-(*am2) + 4*(*am1) + 4*(*a1) - (*a2))/6;
        invvar1_0 = 6*_lambda[1];
      }
      else{
        if (ia[1] == -_K[1]+1 || ia[1] == _K[1]-1){
          if (ia[1] == -_K[1]+1){
            a1  = aa + nr;
            a2  = a1 + nr;
            am1 = aa - nr;            
            mean1_0 = (2*(*am1) + 4*(*a1) - (*a2))/5;
          }
          else{
            a1  = aa + nr;
            am1 = aa - nr;
            am2 = am1 - nr;
            mean1_0 = (-(*am2) + 4*(*am1) + 2*(*a1))/5;
          }
          invvar1_0 = 5*_lambda[1];
        }
        else{
          if (ia[1] == -_K[1]){
            a1  = aa + nr;
            a2  = a1 + nr;
            mean1_0 = 2*(*a1) - (*a2);
          }
          else{
            am1 = aa - nr;
            am2 = am1 - nr;
            mean1_0 = -(*am2) + 2*(*am1);
          }
          invvar1_0 = _lambda[1];
        }
      }
      break;

    case 3:
      /** first conditional (with fixed column) **/
      if (ia[0] >= -_K[0]+3 && ia[0] <= _K[0]-3){
        mean0_1   = (aa[-3] - 6*aa[-2] + 15*aa[-1] + 15*aa[1] - 6*aa[2] + aa[3])/20;
        invvar0_1 = 20*_lambda[0];
      }
      else{
        if (ia[0] == -_K[0]+2 || ia[0] == _K[0]-2){
          if (ia[0] == -_K[0]+2) mean0_1 = (-3*aa[-2] + 12*aa[-1] + 15*aa[1] - 6*aa[2] + aa[3])/19;
          else                   mean0_1 = (aa[-3] - 6*aa[-2] + 15*aa[-1] + 12*aa[1] - 3*aa[2])/19;
          invvar0_1 = 19*_lambda[0];
        }
        else{
          if (ia[0] == -_K[0]+1 || ia[0] == _K[0]-1){
            if (ia[0] == -_K[0]+1) mean0_1 = (3*aa[-1] + 12*aa[1] - 6*aa[2] + aa[3])/10; 
            else                   mean0_1 = (aa[-3] - 6*aa[-2] + 12*aa[-1] + 3*aa[1])/10;
            invvar0_1 = 10*_lambda[0];
          }
          else{
            if (ia[0] == -_K[0]) mean0_1 = 3*aa[1] - 3*aa[2] + aa[3];
            else                 mean0_1 = aa[-3] - 3*aa[-2] + 3*aa[-1];
            invvar0_1 = _lambda[0];
          }  
        }  
      }  

      /** second conditional (with fixed row) **/
      if (ia[1] >= -_K[1]+3 && ia[1] <= _K[1]-3){
        a1  = aa + nr;
        a2  = a1 + nr;
        a3  = a2 + nr;
        am1 = aa - nr;
        am2 = am1 - nr;
        am3 = am2 - nr;
        mean1_0   = ((*am3) - 6*(*am2) + 15*(*am1) + 15*(*a1) - 6*(*a2) + (*a3))/20;
        invvar1_0 = 20*_lambda[1];
      }
      else{
        if (ia[1] == -_K[1]+2 || ia[1] == _K[1]-2){
          if (ia[1] == -_K[1]+2){
            a1  = aa + nr;
            a2  = a1 + nr;
            a3  = a2 + nr;
            am1 = aa - nr;
            am2 = am1 - nr;
            mean1_0 = (-3*(*am2) + 12*(*am1) + 15*(*a1) - 6*(*a2) + (*a3))/19;
          }
          else{
            a1  = aa + nr;
            a2  = a1 + nr;
            am1 = aa - nr;
            am2 = am1 - nr;
            am3 = am2 - nr;
            mean1_0 = ((*am3) - 6*(*am2) + 15*(*am1) + 12*(*a1) - 3*(*a2))/19;
          }
          invvar1_0 = 19*_lambda[1];
        }
        else{
          if (ia[1] == -_K[1]+1 || ia[1] == _K[1]-1){
            if (ia[1] == -_K[1]+1){
              a1  = aa + nr;
              a2  = a1 + nr;
              a3  = a2 + nr;
              am1 = aa - nr;
              mean1_0 = (3*(*am1) + 12*(*a1) - 6*(*a2) + (*a3))/10; 
            }
            else{
              a1  = aa + nr;
              am1 = aa - nr;
              am2 = am1 - nr;
              am3 = am2 - nr;
              mean1_0 = ((*am3) - 6*(*am2) + 12*(*am1) + 3*(*a1))/10;
            }
            invvar1_0 = 10*_lambda[1];
          }
          else{
            if (ia[1] == -_K[1]){
              a1  = aa + nr;
              a2  = a1 + nr;
              a3  = a2 + nr;
              mean1_0 = 3*(*a1) - 3*(*a2) + (*a3);
            }
            else{
              am1 = aa - nr;
              am2 = am1 - nr;
              am3 = am2 - nr;
              mean1_0 = (*am3) - 3*(*am2) + 3*(*am1);
            }
            invvar1_0 = _lambda[1];
          }  
        }  
      }  
      break;

    default:
      REprintf("_order=%d\n", _order);
      throw returnR("Error in BiGspline2_updateWeights.cpp: BiGspline2::full_a_parsBi(...), not (yet) implemented for _order > 3", 1);
    }

    *invvar = invvar0_1 + invvar1_0;        
    *mean   = (invvar0_1*mean0_1 + invvar1_0*mean1_0) / (*invvar);    
    return;

  case BiGspline2A::eight_neighbors:           /***  am1 a4 am4  ***/
                  			       /***  a1  aa a3   ***/
					       /***  am2 a2 am3  ***/

    if (ia[0] > -_K[0] && ia[0] < _K[0] && ia[1] > -_K[1] && ia[1] < _K[1]){     /** not on the edges **/
      a1  = aa - nr;
      am2 = a1 + 1;
      am1 = a1 - 1;
      a3  = aa + nr;
      am3 = a3 + 1;
      am4 = a3 - 1;
      a2  = aa + 1;
      a4  = aa - 1;
      *mean   = (2*(*a1 + *a2 + *a3 + *a4) - (*am1 + *am2 + *am3 + *am4))/4;
      *invvar = 4*_lambda[0];
    }
    else{                                                                /* four edges                                            */
      if (ia[0] == -_K[0]){                                              /* upper edge                                            */
        if (ia[1] == -_K[1]){                                            /* (-K0, -K1) upper left corner                          */
          a3  = aa + nr;
          am3 = a3 + 1;
          a2  = aa + 1;
          *mean   = *a2 + *a3 - *am3;
          *invvar = _lambda[0];
        }
        else{
          if (ia[1] == _K[1]){                                           /* (-K0, K1) upper right corner                          */
            a1  = aa - nr;
            am2 = a1 + 1;
            a2  = aa + 1;
            *mean   = *a1 + *a2 - *am2;
            *invvar = _lambda[0];
          }
          else{                                                          /* (-K0; -K1+1, ..., K1-1) upper edge but not the corner */
            a1  = aa - nr;
            am2 = a1 + 1;
            a3  = aa + nr;
            am3 = a3 + 1;
            a2  = aa + 1;
            *mean = (2*(*a2) + (*a1 + *a3) - (*am2 + *am3))/2;
            *invvar = 2*_lambda[0];
          }
        }
      }                                                                  /* end of the upper edge                                 */
      else{
        if (ia[0] == _K[0]){                                             /* lower edge                                            */
          if (ia[1] == -_K[1]){                                          /* (K0, -K1) lower left corner                           */
            a3  = aa + nr;
            am4 = a3 - 1;
            a4  = aa - 1;
            *mean   = *a3 + *a4 - *am4;
            *invvar = _lambda[0];
          }
          else{
            if (ia[1] == _K[1]){                                         /* (K0, K1) lower right corner                           */
              a1  = aa - nr;
              am1 = a1 - 1;
              a4  = aa - 1;
              *mean   = *a1 + *a4 - *am1;         
              *invvar = _lambda[0];
            }
            else{                                                        /* (K0; -K1+1, ..., K1-1) lower edge but not the corner  */
              a1  = aa - nr;
              am1 = a1 - 1;
              a3  = aa + nr;
              am4 = a3 - 1;
              a4  = aa - 1;
              *mean = (2*(*a4) + (*a1 + *a3) - (*am1 + *am4))/2;
              *invvar = 2*_lambda[0];
            }
          }
        }  /* end of the lower edge */
        else{
          if (ia[1] == -_K[1]){                                          /* left edge but not the corner                          */
            a3  = aa + nr;
            am3 = a3 + 1;
            am4 = a3 - 1;
            a2  = aa + 1;
            a4  = aa - 1;
            *mean = (2*(*a3) + (*a2 + *a4) - (*am3 + *am4))/2;
            *invvar = 2*_lambda[0];
          }
          else{                                                          /* right edge but not the corner                         */
            a1  = aa - nr;
            am2 = a1 + 1;
            am1 = a1 - 1;
            a2  = aa + 1;
            a4  = aa - 1;
            *mean = (2*(*a1) + (*a2 + *a4) - (*am1 + *am2))/2;
            *invvar = 2*_lambda[0];
          }       
        }
      }
    }
    return;

  case BiGspline2A::twelve_neighbors:
    REprintf("_neighbor_system=%d (twelve_neighbors)\n", _neighbor_system);  
    throw returnR("Error in BiGspline2_updateWeights.cpp: BiGspline2::full_a_parsBi(...), not (yet) implemented _neighbor_system", 1);
    return;

  default:
    REprintf("_neighbor_system=%d\n", _neighbor_system);
    throw returnR("Error in BiGspline2_updateWeights.cpp: BiGspline2::full_a_parsBi(...), incorrect _neighbor_system argument", 1);        
  }
}


/***** ============================================================================================================ *****/
/***** find_start_abscisBi:  Find starting abscissae for a particular a in the BIVARIATE case                       *****/
/*****                                                                                                              *****/
/***** ============================================================================================================ *****/
//
//  * use mean and mean +- 3*sd of the distribution [a[ia] | a[-ia], lambda]
//    as starting abscissae
//  * do not check whether they lie on both sides of the mode,
//    this will be done always before sampling is done
//
void
BiGspline2::find_start_abscisBi()
{
  if (BiGspline2A::_dim != 2){
    throw returnR("Error in BiGspline2_updateWeights.cpp: BiGspline2::find_start_abscisBi. Implemented only for BIVARIATE G-splines", 1);
  }
  if (mcmc_Gspline2::_nabscis != 3){
    throw returnR("Dear Arnost, please update BiGspline2::find_start_abscisBi() function after changing _nabscis ;-)", 1);
  }

  double mean, invvar, three_sd;
  const double *aa = _a.aconst();
  double *Abscis   = _abscis.a();
  int ia[2];

  for (ia[1] = -_K[1]; ia[1] <= _K[1]; ia[1]++){
    for (ia[0] = -_K[0]; ia[0] <= _K[0]; ia[0]++){
      this->full_a_parsBi(&mean, &invvar, ia, aa);
      three_sd = 3/sqrt(invvar);
      Abscis[0] = mean - three_sd;
      Abscis[1] = mean;
      Abscis[2] = mean + three_sd;
    
      aa++;
      Abscis += mcmc_Gspline2::_nabscis;
    }
  }

  return;
}
