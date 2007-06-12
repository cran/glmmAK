/*** mcmc_Gspline2.h ***/

#ifndef _MCMC_G_SPLINE_TWEE_H_
#define _MCMC_G_SPLINE_TWEE_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_Optim.h"

#include "Slice_sampler2.h"
#include "ARS2.h"

namespace mcmc_Gspline2 {
  enum a_sampling_scheme {Slice, ARS_quantile, ARS_mode, Block};              // sampling schemes for a

  const double _a_ceil = 10.0;                   // upper limit for 'a' coefficients 
                                                 // * in the UNIVARIATE case applied only if _aIdent = _Reference_
                                                 // * in the BIVARIATE case applied always
  const double _a_floor = -5.0;                  // lower limit for maximal 'a' coefficient
                                                 // * only applied in the BIVARIATE case

  /** constants for adaptive rejection sampling **/
  const int _nabscis = 3;                        // number of starting abscissae for each ARS
  const int _ns = 10;                            // maximal number of abscissae for each ARS
  const int _n_ARS_step = 5;                     // number of sampled points within each ARS (the last one is taken as the sampled value)
  const double _prob[3]   = {0.15, 0.50, 0.85};  // probs for quantiles of the upper hull to compute starting abscissae for the next iteration

  const int _liwv = 7 + mcmc_Gspline2::_ns;         // length of the working array _iwv
  const int _lrwv = 9 + 6*(1+mcmc_Gspline2::_ns);   // length of the working array _rwv


  /** constants for Newton-Raphson within ARS **/
  const int _maxiter_nr = 10;                  // maximum number of NR iterations
  const int _max_stephalf = 10;                // maximum number of step-halving steps
  const double _toler_nr = 1e-3;               // tolerance to detect convergence of NR

  /** constants for Newton-Raphson solver within the slice sampler **/
  const int _maxiter_solver_nr = 10;           // maximum number of NR iterations
  const double _toler_solver_nr = 1e-3;        // tolerance to detect convergence of NR

  /** some other constants **/
  const double _emax      = 64.0;              // exp(-_emax) = 0.0
  const double _epsilon   = exp(-_emax);
  const double _exp_emax  = exp(_emax); 
  const double _null_mass = 1e-6;              // area with this probability will be assigned a zero probability in some applications
                                               //   (point from such area would appear only every 1e6th iteration and I guess we will run 
                                               //   only rarely MCMC with more than 1 000 000 iterations)
  const double _log_inf   = _emax;
             // used in: 'full_a_logdens', 'Gspline2::update_a1' functions to detect problems


  /*** Update of 'a' coefficients ***/
  /*** ========================== ***/
  void
  sample_a_by_slice(double *newa,   double *Abscis,  double *_hx,       double *_hpx,
                    const int *ia,  const int *lia,  const double *aa,  const double *a_pars,  const int *a_ipars,  const int *overrelax);

  void
  sample_a_by_ARS(double *newa,   double *Abscis,  double *_hx,       double *_hpx,          double *_rwv,       int *_iwv,
                  const int *ia,  const int *lia,  const double *aa,  const double *a_pars,  const int *a_ipars,  const int *_type_update_a);


  void
  find_eval_abscis(double *Abscis,  double *_hx,     double *_hpx,  
                   const int *ia,   const int *lia,  const double *aa,  const double* a_pars,  const int* a_ipars);

  void 
  check_abscis(double *Abscis,        double *_hx,        double *_hpx,  
               const double* a_pars,  const int* a_ipars);


  void
  full_a_logdens0(const double *ai,  double *yu,  const double *a_pars,  const int *a_ipars);

  void
  full_a_logdens(const double* ai,  double* yu,  double* ypu,  const double* a_pars,  const int* a_ipars);

  void
  full_a_logdens2(const double* ai,  double* yu,  double* ypu,  double* yppu,  const double* a_pars,  const int* a_ipars);

  void
  full_a_logdens3(const double* ai,  double* yu,  double* ypu,  double* yppu,  const double* a_pars,  const int* a_ipars,  const int& what);


}  /*** end of the namespace mcmc_Gspline2  ***/

#endif


