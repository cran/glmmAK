/*** mcmc_Gspline2.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  27/03/2007
//                            
// PURPOSE: Functions needed for the MCMC with the G-splines
//          * used mainly by objects of classes Gspline2 and BiGspline2
//

//                find_eval_abscic:  15/01/2007, 27/03/2007
//                    check_abscis:  15/01/2007, 27/03/2007
//
//               sample_a_by_slice:  15/01/2007, 27/03/2007
//                 sample_a_by_ARS:  15/01/2007, 27/03/2007
//
//                 full_a_logdens0:  15/01/2007, 27/03/2007
//                  full_a_logdens:  15/01/2007, 27/03/2007
//                 full_a_logdens2:  15/01/2007, 27/03/2007
//                 full_a_logdens3:  15/01/2007, 27/03/2007
// 
/* ********************************************************************************* */

#include "mcmc_Gspline2.h"

namespace mcmc_Gspline2{


/***** ============================================================================================================ *****/
/***** sample_a_by_slice: Sample new a using a slice sampler                                                        *****/
/*****                                                                                                              *****/
/***** ============================================================================================================ *****/
//
// ASSUMPTION: Abscis[1] = mode of the full conditional distribution
//             Abscis[0] = mode - 2*sd
//             Abscis[2] = mode + 2*sd
//
// _hx[mcmc_Gspline2::_nabscis]:  Working array
// _hpx[mcmc_Gspline2::_nabscis]: Working array
//
// ia[lia]:    Index of 'a' that will be updated using computed abscissae
// lia:        Length of ia (1 or 2)
//  aa:        Current value of a
//
void
sample_a_by_slice(double *newa,   double *Abscis,  double *_hx,       double *_hpx,
                  const int *ia,  const int *lia,  const double *aa,  const double *a_pars,  const int *a_ipars,  const int *overrelax)
{
  static double horiz;
  static int i, iter_nr, err_nr;
  static double *hx, *hpx;

  hx       = _hx;
  hpx      = _hpx;

  /*** Reshuffle Abscis such that the interval defining the slice will be at the beginnnig      ***/
  /***  and evaluate h in current point - store it in _hx[2]                                    ***/  
  Abscis[1] = Abscis[2];
  hx[1]     = hx[2];
  hpx[1]    = hpx[2];
  mcmc_Gspline2::full_a_logdens0(aa, hx+2, a_pars, a_ipars);    

  /*** Sample the horizontal level defining the slice (on log-scale), store it in horiz ***/
  horiz = hx[2] - rexp(1);   

  /*** Find the interval defining the slice ***/
  for (i = 0; i < 1; i++){
    try{
      AK_Optim::solver_newton_raphson(Abscis + i, hx + i, hpx + i, &horiz, a_pars, a_ipars, mcmc_Gspline2::full_a_logdens3,
                     &iter_nr, &mcmc_Gspline2::_maxiter_solver_nr, &mcmc_Gspline2::_toler_solver_nr, &mcmc_Gspline2::_epsilon, &err_nr);
    }
    catch(returnR){
      if (*lia == 1) REprintf("a[%d] = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
                              *ia, *aa, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
      if (*lia == 2) REprintf("a[%d,%d] = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
                              ia[0], ia[1], *aa, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
      REprintf("Trap in mcmc_Gspline2.cpp: sample_a_by_slice(). solver_newton_raphson failed\n");
      throw;
    }
    if (err_nr >= 3){
      REprintf("err_nr = %d\n", err_nr);
      if (*lia == 1) REprintf("a[%d] = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
                              *ia, *aa, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
      if (*lia == 2) REprintf("a[%d,%d] = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
                              ia[0], ia[1], *aa, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
      throw returnR("Trap in mcmc_Gspline2.cpp: sample_a_by_slice(). Unable to find an interval defining the slice", 1);
    }
  }

  /*** Sample the new point ***/
  if (*overrelax){
    Slice_sampler::ss_exact_overrelax(newa, Abscis, aa, &horiz, mcmc_Gspline2::full_a_logdens0, a_pars, a_ipars);
  }
  else{
    Slice_sampler::ss_exact_sample(newa, Abscis, hx, aa, &horiz, mcmc_Gspline2::full_a_logdens0, a_pars, a_ipars);
  }

  return;
}


/***** ============================================================================================================ *****/
/***** sample_a_by_ARS: Sample new 'a' using the ARS                                                                *****/
/*****                                                                                                              *****/
/***** ============================================================================================================ *****/
//
// _rwv[mcmc_Gspline2::_lrwv]:   Working array
// _iwv[mcmc_Gspline2::_liwv]:   Working array
//
void
sample_a_by_ARS(double *newa,   double *Abscis,  double *_hx,       double *_hpx,          double *_rwv,        int *_iwv,
                const int *ia,  const int *lia,  const double *aa,  const double *a_pars,  const int *a_ipars,  const int *_type_update_a)
{
  static int i, ifault, r_zero;
  static double hlb, hub;
  static double *hx, *hpx, *rwv;
  static int *iwv;
  hx   = _hx;
  hpx  = _hpx;
  rwv  = _rwv;
  iwv  = _iwv;


  /***  Initialize arrays for ARS, ifault can be either 0 or 5, other values (1, 2, 3, 4) are not possible  ***/
  ifault = 1;
  ARS::initial_(&mcmc_Gspline2::_ns, &mcmc_Gspline2::_nabscis, &mcmc_Gspline2::_emax, 
                Abscis, hx, hpx, &_AK_ZERO_INT, &hlb, &_AK_ZERO_INT, &hub, &ifault, iwv, rwv);
  if (ifault > 0){    /*  numerical non-log-concavity detected -> use slice sampler instead */
    mcmc_Gspline2::sample_a_by_slice(newa, Abscis, _hx, _hpx, ia, lia, aa, a_pars, a_ipars, &_AK_ZERO_INT);
    return;
//    if (*lia == 1) REprintf("a[%d] = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
//                            *aa, *ia, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
//    if (*lia == 2) REprintf("a[%d,%d] = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
//                            *aa, ia[0], ia[1], a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
//    throw returnR("Trap in ARS: Numerical non-log-concavity (initial_) detected when updating 'a' coefficients, use slice sampler instead", 1);
  }

  /***  Sample the new point ***/
  for (i = 0; i < mcmc_Gspline2::_n_ARS_step; i++){
    r_zero = 0;
    ifault = 6;    
    while (ifault == 6){    /* while random number generator generates zero */
      ARS::sample_(iwv, rwv, mcmc_Gspline2::full_a_logdens, a_pars, a_ipars, newa, &ifault);
      switch (ifault){
      case 5:    /*  numerical non-log-concavity detected -> use slice sampler instead */
        mcmc_Gspline2::sample_a_by_slice(newa, Abscis, _hx, _hpx, ia, lia, aa, a_pars, a_ipars, &_AK_ZERO_INT);
        return;
//        if (*lia == 1) REprintf("a[%d] = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
//                                *aa, *ia, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
//        if (*lia == 2) REprintf("a[%d,%d] = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
//                                *aa, ia[0], ia[1], a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
//        throw returnR("Trap in ARS: Numerical non-log-concavity (initial_) detected when updating 'a' coefficients, use slice sampler instead", 1);
      case 6:
        r_zero++;
        Rprintf("Warning: Random number generator generated zero during ARS.\n");
        if (r_zero >= 10) throw returnR("Trap in ARS: Too many zeros generated by the random number generator", 1);
        break;
      case 7:
        throw returnR("Trap in ARS: Numerical instability detected by sample_", 1);
      }
    }
  }

  /***  Compute starting abscissae for the next iteration of MCMC (if necessary)  ***/
  if (*_type_update_a == mcmc_Gspline2::ARS_quantile) 
    ARS::quantile_(iwv, rwv, &mcmc_Gspline2::_nabscis, mcmc_Gspline2::_prob + 0, Abscis, &_AK_ZERO_INT);
  return;
}


/***** ============================================================================================================ *****/
/***** find_eval_abscis: Find starting abscissae for the ARS or initial guess for the interval defining the slice   *****/
/*****                                                                                                              *****/
/***** ============================================================================================================ *****/
//
// it modifies: Abscis:
//    on OUTPUT: Abscis[0] = mode of full conditional - 2*sd
//               Abscis[1] = mode of full conditional
//               Abscis[2] = mode of full conditional + 2*sd
//
// _hx[mcmc_Gspline2::_nabscis]:  Working array
// _hpx[mcmc_Gspline2::_nabscis]: Working array
//
// ia[lia]:       index of 'a' that will be updated using computed abscissae
// lia:           length of ia (1 or 2)
// aa:            current value of 'a' that is going to be updated
// a_pars[4]:     double parameters for the full conditional distribution of a
// a_ipars[2]:    integer parameters for the full conditional distribution of a
//
void
find_eval_abscis(double *Abscis,  double *_hx,     double *_hpx,  
                 const int *ia,   const int *lia,  const double *aa,  const double* a_pars,  const int* a_ipars)
{
  static double hppx;
  static int err_nr, iter_nr;

  static double *hx, *hpx;
  hx  = _hx;
  hpx = _hpx;

  Abscis[1] = *aa;
  mcmc_Gspline2::full_a_logdens3(Abscis + 1, hx + 1, hpx + 1, &hppx, a_pars, a_ipars, 0);
  try{
    AK_Optim::optim_newton_raphson02(Abscis + 1, hx + 1, hpx + 1, &hppx, a_pars, a_ipars, mcmc_Gspline2::full_a_logdens3, 
        &iter_nr, &mcmc_Gspline2::_maxiter_nr, &mcmc_Gspline2::_max_stephalf, &mcmc_Gspline2::_toler_nr, &mcmc_Gspline2::_epsilon, &err_nr);
  }
  catch(returnR){
    if (*lia == 1) REprintf("a[%d] = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
                            *ia, *aa, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    if (*lia == 2) REprintf("a[%d,%d] = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
                            ia[0], ia[1], *aa, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    REprintf("Trap in mcmc_Gspline2.cpp: find_eval_abscis(). optim_newton_raphson02 failed\n");
    throw;  
  }
  if (err_nr >= 3){
    REprintf("err_nr = %d\n", err_nr);
    if (*lia == 1) REprintf("a[%d] = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
                            *ia, *aa, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    if (*lia == 2) REprintf("a[%d,%d] = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
                            ia[0], ia[1], *aa, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    throw returnR("Trap in mcmc_Gspline2.cpp. find_eval_abscis(). Unable to find a mode of the full conditional distribution", 1);
  }
  if (hppx <= mcmc_Gspline2::_epsilon) hppx = mcmc_Gspline2::_epsilon;
  hppx = 2/sqrt(hppx);                                            /* approx. 2*std. deviation of the full conditional */
  Abscis[0] = Abscis[1] - hppx;
  Abscis[2] = Abscis[1] + hppx;
  mcmc_Gspline2::full_a_logdens(Abscis + 0, hx + 0, hpx + 0, a_pars, a_ipars);
  mcmc_Gspline2::full_a_logdens(Abscis + 2, hx + 2, hpx + 2, a_pars, a_ipars);

  return;
}


/***** ============================================================================================================ *****/
/***** check_abscis: Check whether starting abscissae/initial guesses                                               *****/
/*****               for the interval defining the slice lie on correct size of the mode                            *****/
/*****                                                                                                              *****/
/***** ============================================================================================================ *****/
//
// Abscis:                        It is possibly shifted to satisfy the condition that it should lie on a correct size of the mode
// _hx[mcmc_Gspline2::_nabscis]:  Working array
// _hpx[mcmc_Gspline2::_nabscis]: Working array
// a_pars[4]:                     Double parameters for the full conditional distribution of a
// a_ipars[2]:                    Integer parameters for the full conditional distribution of a
//
void 
check_abscis(double *Abscis,        double *_hx,        double *_hpx,  
             const double* a_pars,  const int* a_ipars)
{
  static double *abscis_E;
  static double *hx, *hpx, *hx_E, *hpx_E;
  static double step_left, step_right;
  static bool left_bad, right_bad;

  hx       = _hx;
  hpx      = _hpx;
  abscis_E = Abscis + mcmc_Gspline2::_nabscis-1;
  hx_E     = hx  + mcmc_Gspline2::_nabscis - 1;
  hpx_E    = hpx + mcmc_Gspline2::_nabscis - 1;

  step_left  = Abscis[1] - Abscis[0];
  step_right = abscis_E[0] - abscis_E[-1];
  left_bad = right_bad = true;

  while (left_bad){
    if (hpx[0] < mcmc_Gspline2::_epsilon){
      Abscis[0] -= step_left;
      mcmc_Gspline2::full_a_logdens(Abscis + 0, hx + 0, hpx + 0, a_pars, a_ipars);
    }
    else
      left_bad = false;
  }

  while (right_bad){
    if (hpx_E[0] > -mcmc_Gspline2::_epsilon){
      abscis_E[0] += step_right;
      mcmc_Gspline2::full_a_logdens(abscis_E, hx_E, hpx_E, a_pars, a_ipars);
    }
    else
      right_bad = false;
  }
  return;
}


/***** ============================================================================================================ *****/
/***** Log-density of the full conditional distribution [a[ia] | a[-ia], lambda]                                    *****/
/***                                                                                                                  ***/ 
/*** full_a_logdens0: Compute only log-density                                                                        ***/
/*** full_a_logdens:  Compute log-density and its first derivative of the log of the full conditional of a            ***/
/*** full_a_logdens2: Compute also  minus second derivative of log-density                                            ***/
/*** full_a_logdens3: Compute selectively log-density or/and first and minus second derivative of log-density         ***/
/***                                                                                                                  ***/ 
/***** ============================================================================================================ *****/
//
// ai ......... new a, point where it should be evaluated
// yu ......... log-density
// ypu ........ derivative of the log-density
// yppu ....... minus second derivative of the log-density
// a_pars ..... a_pars[0] = mean of [a[ia] | a[-ia], lambda]
//              a_pars[1] = inverse variance of [a[ia] | a[-ia], lambda]
//              a_pars[2] = exp(old ai) 
//              a_pars[3] = sum(exp(old a)) (evaluated at current a vector)
// a_ipars .... a_ipars[0] = N     = number of all observations
//              a_ipars[1] = N[ia] = number of observations in a component ia
//
void
full_a_logdens0(const double *ai,  double *yu,  const double *a_pars,  const int *a_ipars)
{
  static double new_expai, new_sumexpa, a_min_A;
  if (*ai >= mcmc_Gspline2::_log_inf){
    //REprintf("\na = %e leads to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //         *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    //throw returnR("Trap in full_a_logdens0.", 1);
    //Rprintf("\na = %e may lead to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //        *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    new_expai   = mcmc_Gspline2::_exp_emax;
    new_sumexpa = mcmc_Gspline2::_exp_emax;
  }
  else{
    new_expai   = exp(*ai);
    new_sumexpa = a_pars[3] - a_pars[2] + new_expai;   
  }

  a_min_A = (*ai) - a_pars[0];

  *yu = a_ipars[1]*(*ai) - a_ipars[0]*log(new_sumexpa) - 0.5*a_pars[1]*a_min_A*a_min_A;

  if (!R_finite(*yu)){
    REprintf("\na = %e, yu = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
             *ai, *yu, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    throw returnR("Trap in mcmc_Gspline2.cpp: full_a_logdens0(). NaN is not allowed.", 1);
  }
  return;
}

void
full_a_logdens(const double* ai,  double* yu,  double* ypu,  const double* a_pars,  const int* a_ipars)
{
  static double new_expai, new_sumexpa, a_min_A, new_wi;
  if (*ai >= mcmc_Gspline2::_log_inf){
    //REprintf("\na = %e leads to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //         *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    //throw returnR("Trap in full_a_logdens.", 1);
    //Rprintf("\na = %e may lead to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //        *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    new_expai   = mcmc_Gspline2::_exp_emax;
    new_sumexpa = mcmc_Gspline2::_exp_emax;
  }
  else{
    new_expai   = exp(*ai);
    new_sumexpa = a_pars[3] - a_pars[2] + new_expai;   
  }

  a_min_A = (*ai) - a_pars[0];
  new_wi  = new_expai / new_sumexpa;

  *yu  = a_ipars[1]*(*ai) - a_ipars[0]*log(new_sumexpa) - 0.5*a_pars[1]*a_min_A*a_min_A;
  *ypu = a_ipars[1] - a_ipars[0]*new_wi - a_pars[1]*a_min_A;

  if (!R_finite(*yu)){
    REprintf("\na = %e, yu = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
             *ai, *yu, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    throw returnR("Trap in mcmc_Gspline2.cpp: full_a_logdens(). NaN is not allowed.", 1);
  }
  if (!R_finite(*ypu)){
    REprintf("\na = %e, yu = %e, ypu = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
             *ai, *yu, *ypu, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    throw returnR("Trap in mcmc_Gspline2.cpp: full_a_logdens(). NaN is not allowed.", 1);
  }
  return;
}

void
full_a_logdens2(const double* ai,  double* yu,  double* ypu,  double* yppu,  const double* a_pars,  const int* a_ipars)
{
  static double new_expai, new_sumexpa, a_min_A, new_wi;

  if (*ai >= mcmc_Gspline2::_log_inf){
    //REprintf("\na = %e leads to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //         *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    //throw returnR("Trap in full_a_logdens2.", 1);
    //Rprintf("\na = %e may lead to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //        *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    new_expai   = mcmc_Gspline2::_exp_emax;
    new_sumexpa = mcmc_Gspline2::_exp_emax;
  }
  else{
    new_expai   = exp(*ai);
    new_sumexpa = a_pars[3] - a_pars[2] + new_expai;   
  }

  a_min_A = (*ai) - a_pars[0];
  new_wi  = new_expai / new_sumexpa;

  *yu   = a_ipars[1]*(*ai) - a_ipars[0]*log(new_sumexpa) - 0.5*a_pars[1]*a_min_A*a_min_A;
  *ypu  = a_ipars[1] - a_ipars[0]*new_wi - a_pars[1]*a_min_A;
  *yppu = a_ipars[0]*new_wi*(1 - new_wi) + a_pars[1];

  if (!R_finite(*yu)){
    REprintf("\na = %e, yu = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
             *ai, *yu, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    throw returnR("Trap in mcmc_Gspline2.cpp: full_a_logdens2(). NaN is not allowed.", 1);
  }
  if (!R_finite(*ypu)){
    REprintf("\na = %e, yu = %e, ypu = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
             *ai, *yu, *ypu, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    throw returnR("Trap in mcmc_Gspline2.cpp: full_a_logdens2(). NaN is not allowed.", 1);
  }
  if (!R_finite(*yppu)){
    REprintf("\na = %e, yu = %e, ypu = %e, yppu = %e, pars[0] = %e, pars[1] = %e, pars[2] = %e, pars[3] = %e, ipars[0] = %d, ipars[1] = %d \n", 
             *ai, *yu, *ypu, *yppu, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    throw returnR("Trap in mcmc_Gspline2.cpp: full_a_logdens2(). NaN is not allowed.", 1);
  }
  return;
}

//
// what .... 0 = compute l(x), l'(x), -l''(x)
//           1 = compute only l(x)
//           2 = compute only l'(x), -l''(x)
//           3 = compute only l(x), l'(x)
//
void
full_a_logdens3(const double* ai,  double* yu,  double* ypu,  double* yppu,  const double* a_pars,  const int* a_ipars,  const int& what)
{
  static double new_expai, new_sumexpa, a_min_A, new_wi;

  if (*ai >= mcmc_Gspline2::_log_inf){
    //Rprintf("\na = %e may lead to exp(a) = Inf; a_pars = %e, %e, %e, %e;  a_ipars = %d, %d\n", 
    //        *ai, a_pars[0], a_pars[1], a_pars[2], a_pars[3], a_ipars[0], a_ipars[1]);
    new_expai   = mcmc_Gspline2::_exp_emax;
    new_sumexpa = mcmc_Gspline2::_exp_emax;
  }
  else{
    new_expai   = exp(*ai);
    new_sumexpa = a_pars[3] - a_pars[2] + new_expai;
  }
  a_min_A = (*ai) - a_pars[0];
  new_wi = new_expai / new_sumexpa;

  switch (what){
  case 0:
    *yu   = a_ipars[1]*(*ai) - a_ipars[0]*log(new_sumexpa) - 0.5*a_pars[1]*a_min_A*a_min_A;
    *ypu  = a_ipars[1] - a_ipars[0]*new_wi - a_pars[1]*a_min_A;
    *yppu = a_ipars[0]*new_wi*(1 - new_wi) + a_pars[1];
    return;
  case 1:
    *yu = a_ipars[1]*(*ai) - a_ipars[0]*log(new_sumexpa) - 0.5*a_pars[1]*a_min_A*a_min_A;
    return;
  case 2:
    *ypu  = a_ipars[1] - a_ipars[0]*new_wi - a_pars[1]*a_min_A;
    *yppu = a_ipars[0]*new_wi*(1 - new_wi) + a_pars[1];
    return;
  case 3:
    *yu  = a_ipars[1]*(*ai) - a_ipars[0]*log(new_sumexpa) - 0.5*a_pars[1]*a_min_A*a_min_A;
    *ypu = a_ipars[1] - a_ipars[0]*new_wi - a_pars[1]*a_min_A;
    return;
  default:
    throw returnR("Error in mcmc_Gspline2.cpp: full_a_logdens3(). Incorrect what argument.", 1);
  }
}

}   /*** end of the namespace mcmc_Gspline2 ***/

