/*** GMRF_Gspline.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  08/11/2006
//
//    PURPOSE: Functions that implement (intrinsic) Gaussian random fields
//             and also utilities for them useful in the G-spline context
//
//                 update:  18/11/2006,  22/11/2006
//                 ML_est:  09/11/2006,  22/11/2006
//                    ll0:  08/11/2006
//                    ll1:  08/11/2006,  21/11/2006
//                    ll2:  08/11/2006,  21/11/2006
//      mcmc_GMRF_Gspline:  22/11/2006
//                          
/* ********************************************************* */

#include "GMRF_Gspline.h"

namespace GMRF_Gspline {

/*** ================================================================================== ***/
/*** Iteration of the MCMC which updates both GMRF and its precision                    ***/
/***                                                                                    ***/
/***                                                                                    ***/
/*** ================================================================================== ***/
//
// INPUT:
// ======
//  a, lambda, expa, w, minw, sumexpa, Da, Qa, min_half_aQa:  
//     Current values of these quantities
//
// OUTPUT:
// =======
// accept[1]:        Acceptance indicator
// a[na]:            Values of the GMRF (transformed mixture weights)
// lambda[1]:        Value of the GMRF precision
// expa[na]:         Values of exp(a[j])
// w[na]:            Weights = exp(a[j])/sum(exp(a[j]))
// minw[1]:          Value of the lowest weight
// sumexpa[1]:       sum(exp(a[j]))
// Da[na]:           Value of D*a in the first na-order places, 
//                   for working purposes, array of length na is needed
// Qa[na]:           Value of Q*a
// min_half_aQa[1]:  Value of -0.5*a'Q*a (penalty)
//
// WORKING SPACE:
// ==============
// workML[]:       See workD argument of 'ML_est_GMRF_Gspline' function
// worka[]:        See below section "Set up pointers in worka"
// workGMRF[]:     See below section "Set up pointers in workGMRF"
//                 * length the same as in 'rGMRFR' function in GMRF.cpp
//
// OBSERVATIONS RELATED PARAMETERS:
// ================================
// allocN[na]:     Numbers of observations allocated in each mixture component 
//
// PARAMETERS OF THE LAMBDA-PRIOR AND ITS PROPOSAL:
// ================================================
// prior_for_lambda[1]:   Type of the prior for lambda
//                        See enum priorForLambdaTypes in GMRF_Gspline_Util.h
// par_lambda[2]:         if prior_for_lambda = _Gamma_:
//                           par_lambda[0] = prior shape
//                           par_lambda[1] = prior rate
//                        if prior_for_lambda = _SDUnif_:                        
//                           par_lambda[0] = arbitrary
//                           par_lambda[1] = 1/S^2, where S is the upper limit of the uniform prior on 1/sqrt(lambda)
// par_rscale[6]:         Parameters for the lambda-proposal.
//                        See functions 'dscale' and 'rscale' in GMRF.cpp
//
// PARAMETERS OF THE GMRF:
// =======================
//  Q[LT(na)]:          Lower triangle of the matrix Q = t(D) * D
//  order[1]:           Order of the differences in the GMRF prior
//  diffOper[order+1]:  Vector defining the difference operator, e.g.,
//                        order = 0:   1
//                        order = 1:   -1, 1
//                        order = 2:   1, -2, 1
//                        order = 3:   -1, 3, -3, 1
//                         etc.
//  epsw[1]:            Value which is added to the diagonal of the minus Hessian if (nobs*minw < epsw)
//
// OTHER PARAMETERS:
// =================
//  constraint[1]:      Type of the identifiability constraint (see GMRF_Gspline_Util.h for possible values)
//  iref[1]:            Index of a coefficient which is expressed as a function of d (remaining a coefficients)
//
// na[1]:               Length of the GMRG
// nobs[1]:             Number of observations used ti fit the GMRF
//                      Identity sum(allocN) = nobs should hold
//
// lambda_a_block[1]:   If <> 0 then (lambda, a) is sampled using the Metropolis-Hastings step and accepted/rejected jointly
//                      If = 0  then lambda is updated using the Gibbs move
//                              and a is sampled using the Metropolis-Hastings step
//
void
update(int *accept,               double *a,                    double *lambda,  
       double *expa,              double *sumexpa,              double *w,                  double *minw,
       double *Da,                double *Qa,                   double *min_half_aQa,   
       double *workML,            double *worka,                double *workGMRF,
       const int *allocN,         const int *prior_for_lambda,  const double *par_lambda,   const double* par_rscale,
       const double *Q,           const int *order,             const int *diffOper,        const double *epsw,
       const int *constraint,     const int *iref,              const int *na,              const int *nobs,
       const int *lambda_a_block)
{ 
  static const int _ZERO_ = 0;

  static const int nc[1] = {0};              /** Number of identifiability constraints in the vector d                                      **/
  static const int LTnc = 0;                 /** = (*nc)*(*nc+1)/2;                                                                         **/
  static const double A[1] = {0};            /** Left-hand side of the identifiability constraint (not used here)                           **/
  static const double e[1] = {0};            /** Right-hand side of the identifiability constraint (not used here)                          **/
  static const int e_nonZERO[1] = {0};       /** Indicator whether the RHS of the identifiability constraint is zero or not (not used here) **/
  //static const int mu_nonZERO0[1] = {0};     /** Indicator whether the prior mean of the GMRG is zero or not                                **/
  static const int mu_nonZERO1[1] = {1};     /** Indicator whether the prior mean of the GMRG is zero or not                                **/


  static int err[1], niter[1], Attempt[1]; 
  static double flambda[1], prop_lambda[1], ll[1], log_proposal_prop_a[1], log_proposal_a[1];
  static double lshape, lscale, log_A;

  const int na_order = (*na) - (*order);
  //const int LTna     = (*na)*(*na+1)/2;
  const int na_1     = *na - 1;
  const int LTna_1   = na_1*(*na)/2;

  /** Set-up pointers in workML **/
  static double *expaML, *sumexpaML, *wML, *minwML, *DaML, *QaML, *min_half_aQaML, *NR_stepML, *ddll_tempML, *workll2ML;
  expaML         = workML;
  sumexpaML      = expaML         + (*na);
  wML            = sumexpaML      + 1;
  minwML         = wML            + (*na);
  DaML           = minwML         + 1;
  QaML           = DaML           + (*na);
  min_half_aQaML = QaML           + (*na);
  NR_stepML      = min_half_aQaML + 1;
  ddll_tempML    = NR_stepML      + na_1;
  workll2ML      = ddll_tempML    + LTna_1;
  //next         = workll2ML + na + LTna;

  /** Set up pointers in worka **/
  static double *prop_mu, *prop_iVar, *dll, *prop_a, *prop_expa, *prop_sumexpa, *prop_w, *prop_minw, *prop_Da, *prop_Qa, *prop_min_half_aQa;
  prop_mu           = worka;
  prop_iVar         = prop_mu      + *na;
  dll               = prop_iVar    + LTna_1;
  prop_a            = dll          + na_1;
  prop_expa         = prop_a       + (*na);
  prop_sumexpa      = prop_expa    + (*na);
  prop_w            = prop_sumexpa + 1;
  prop_minw         = prop_w       + (*na);
  prop_Da           = prop_minw    + 1;
  prop_Qa           = prop_Da      + (*na);
  prop_min_half_aQa = prop_Qa      + (*na);
  //next         = prop_min_half_aQaML + 1;
  
  /** Set up pointers in workGMRF **/
  static double *log_dens_Ax_x, *log_dets, *U, *LW, *workGMRF2;
  log_dens_Ax_x = workGMRF;
  log_dets      = log_dens_Ax_x + 1;
  U             = log_dets      + 4;
  LW            = U             + (*nc)*(*na);
  workGMRF2     = LW            + LTnc;          /** working array for log_density_Ax_x, rGMRF_inputArgs, rGMRF, dGMRF_inputArgs, dGMRF **/

  /** Propose new value of lambda (lambda*) and compute appropriate part of the log(acceptance ratio) **/
  //Rprintf("\nInitial lambda=%g,", *lambda);
  switch (*prior_for_lambda){
  case GMRF_Gspline_Util::_Fixed_:
    *prop_lambda = *lambda;

    log_A = 0;
    break;

  case GMRF_Gspline_Util::_Gamma_:
    if (*lambda_a_block){
      GMRF::rscale(flambda, par_rscale);
      *prop_lambda = *flambda * (*lambda);

      log_A = (par_lambda[0] + 0.5*(na_order+1)) * (log_AK(prop_lambda[0]) - log_AK(*lambda));
      log_A += -par_lambda[1] * (*prop_lambda - (*lambda));
    }
    else{
      lshape = par_lambda[0] + 0.5*(na_order+1);
      lscale  = 1/(par_lambda[1] - (*min_half_aQa));
      *prop_lambda = rgamma(lshape, lscale);
      *lambda = *prop_lambda;

      log_A = 0;
    }
    break;

  case GMRF_Gspline_Util::_SDUnif_:
    if (*lambda_a_block){
      throw returnR("Error in GMRF_Gspline.cpp: update. Block update of (lambda, a) not implemented for _SDUnif_ prior on lambda.", 1);
    }
    else{
      lshape = 0.5*na_order;
      lscale  = 1/(- (*min_half_aQa));
      GMRF_Gspline_Util::rltruncGamma(prop_lambda, &lshape, &lscale, par_lambda+1);
      *lambda = *prop_lambda;

      log_A = 0;
    }
    break;

  default:
    throw returnR("Error in GMRF_Gspline.cpp: update. Unknown value of prior_for_lambda.", 1);
  }
  //Rprintf("   Proposed lambda=%g\n", *prop_lambda);

  /** Construct normal approximation to p(a* | lambda*, r), mean will be stored in 'prop_mu' (in first na-1 places) and precision in 'prop_iVar' **/
  GMRF_Gspline_Util::copy_within_update(expaML, sumexpaML, wML, minwML, DaML, QaML, min_half_aQaML, prop_mu,
                                        expa,   sumexpa,   w,   minw,   Da,   Qa,   min_half_aQa,   a,
                                        constraint,  iref,  na,  &na_1, GMRF_Gspline_Util::_a2a_);

  GMRF_Gspline::ll2(ll, dll, prop_iVar, workll2ML, prop_mu, allocN, prop_lambda, sumexpaML, min_half_aQaML, QaML, wML, Q, order, constraint, iref, na, nobs);
  //Rprintf("\nAllocN: ");
  //AK_BLAS_LAPACK::printIArray(allocN, *na);  
  //Rprintf("\nInitial a (log-lik=%g,  lambda=%g):\n", *ll, *lambda);
  //AK_BLAS_LAPACK::printArray(a, *na);  
  //Rprintf("\nInitial score (log-lik=%g):\n", *ll);
  //AK_BLAS_LAPACK::printArray(dll, na_1);
  //Rprintf("\nInitial Hessian:\n");
  //AK_BLAS_LAPACK::printLT4R(prop_iVar, na_1);
  GMRF_Gspline::ML_est(ll, dll, prop_iVar, prop_mu, workML, niter, err, allocN, prop_lambda, Q, order, diffOper, constraint, iref, na, nobs, &GMRF_Gspline::_maxiter, epsw);
  if (*err >= 3){
    REprintf("\na: ");
    AK_BLAS_LAPACK::printArray(a, *na);
    REprintf("w: ");
    AK_BLAS_LAPACK::printArray(w, *na);
    throw returnR("Trap in GMRF_Gspline.cpp: update. Not possible to find a normal approximation.", 1);
  }
  //Rprintf("\nMode of the normal approximation (niter=%d):\n", *niter);
  //AK_BLAS_LAPACK::printArray(prop_mu, *na);
  GMRF_Gspline_Util::a2d(prop_mu, iref, na);

  /** Propose new value of a by sampling from the normal approximation **/
  AK_BLAS_LAPACK::chol_dpptrfPD(prop_iVar, ddll_tempML, &na_1, Attempt, &GMRF_Gspline::_max_nAttempt, epsw, err);
  if (*err){
    REprintf("\na: ");
    AK_BLAS_LAPACK::printArray(a, *na);
    REprintf("w: ");
    AK_BLAS_LAPACK::printArray(w, *na);
    throw returnR("Trap in GMRF_Gspline.cpp: update. Singular precision matrix of the normal approximation.", 1);
  }
  *log_dens_Ax_x = 0;
  GMRF::rGMRF_inputArgs(log_dets, U, prop_mu, prop_iVar, A, e, &na_1, nc, mu_nonZERO1, e_nonZERO, workGMRF2);
  //Rprintf("\nNormal approximation.  log(|Q|^{1/2}) = %g,  -(nx - nc)*log(sqrt(2pi))= %g\n", log_dets[0], log_dets[1]);
  GMRF::rGMRF(prop_a, log_proposal_prop_a, prop_mu, prop_iVar, log_dets, A, e, U, log_dens_Ax_x, &na_1, nc, mu_nonZERO1, e_nonZERO, workGMRF2);
  GMRF_Gspline_Util::d2a(prop_a, constraint, iref, na);

  /** Construct reversal normal approximation, that is normal approximation to p(a | lambda, r) **/   
  GMRF_Gspline_Util::update4_ll12(prop_expa, prop_sumexpa, prop_Da, prop_min_half_aQa, prop_Qa, prop_w, prop_minw, prop_a, 
                                  order, diffOper, na);
  //Rprintf("\nProposed a:\n");
  //AK_BLAS_LAPACK::printArray(prop_a, *na);
  GMRF_Gspline_Util::copy_within_update(expaML,    sumexpaML,    wML,    minwML,    DaML,    QaML,    min_half_aQaML,    prop_mu,
                                        prop_expa, prop_sumexpa, prop_w, prop_minw, prop_Da, prop_Qa, prop_min_half_aQa, prop_a,
                                        constraint,  iref,  na,  &na_1,  GMRF_Gspline_Util::_a2a_);
  GMRF_Gspline::ll2(ll, dll, prop_iVar, workll2ML, prop_mu, allocN, lambda, sumexpaML, min_half_aQaML, QaML, wML, Q, order, constraint, iref, na, nobs);
  GMRF_Gspline::ML_est(ll, dll, prop_iVar, prop_mu, workML, niter, err, allocN, lambda, Q, order, diffOper, constraint, iref, na, nobs, &GMRF_Gspline::_maxiter, epsw);
  if (*err >= 3){
    REprintf("\na: ");
    AK_BLAS_LAPACK::printArray(a, *na);
    REprintf("prop_a: ");
    AK_BLAS_LAPACK::printArray(prop_a, *na);
    REprintf("w: ");
    AK_BLAS_LAPACK::printArray(w, *na);
    REprintf("prop_w: ");
    AK_BLAS_LAPACK::printArray(prop_w, *na);
    throw returnR("Trap in GMRF_Gspline.cpp: update. Not possible to find a normal approximation of the reversible proposal.", 1);
  }
  //Rprintf("\nMode of the reversible normal approximation (niter=%d):\n", *niter);
  //AK_BLAS_LAPACK::printArray(prop_mu, *na);
  GMRF_Gspline_Util::a2d(prop_mu, iref, na);

  /** Evaluate log-density of the reversible proposal in current 'a' **/
  AK_BLAS_LAPACK::chol_dpptrfPD(prop_iVar, ddll_tempML, &na_1, Attempt, &GMRF_Gspline::_max_nAttempt, epsw, err);
  if (*err){
    REprintf("\na: ");
    AK_BLAS_LAPACK::printArray(a, *na);
    REprintf("prop_a: ");
    AK_BLAS_LAPACK::printArray(prop_a, *na);
    REprintf("w: ");
    AK_BLAS_LAPACK::printArray(w, *na);
    REprintf("prop_w: ");
    AK_BLAS_LAPACK::printArray(prop_w, *na);
    throw returnR("Trap in GMRF_Gspline.cpp: update. Singular precision matrix of the normal approximation of the reversible proposal.", 1);
  }
  //*log_dens_Ax_x = 0;
  GMRF::dGMRF_inputArgs(log_dets, LW, U, prop_mu, prop_iVar, A, e, &na_1, nc, mu_nonZERO1, e_nonZERO, workGMRF2);
  //Rprintf("\nReversible normal approximation.  log(|Q|^{1/2}) = %g,  -(nx - nc)*log(sqrt(2pi))= %g\n", log_dets[0], log_dets[1]);
  GMRF_Gspline_Util::a2d2(dll, a, iref, na);                                                           /** Store d part of a(d) in dll **/  
  GMRF::dGMRF(log_proposal_a, dll, &_ZERO_, prop_mu, prop_iVar, log_dets, log_dens_Ax_x, &na_1, nc, mu_nonZERO1, workGMRF2);

  /** Log-acceptance ratio, see document 'block-update-GMRF.pdf' **/
  /** (its lambda part has already been computed above)          **/
  log_A += -(*nobs) * (log_AK(*prop_sumexpa) - log_AK(*sumexpa));

  AK_BLAS_LAPACK::c_aMinusb(prop_mu, prop_a, a, *na);                 /* a(star) - a      */
  AK_BLAS_LAPACK::ddot(ll, prop_mu, allocN, *na);                     /* (a(star) - a)'N  */
  log_A += *ll;

  if (*lambda_a_block){
    log_A += ((*prop_lambda)*(*prop_min_half_aQa) - (*lambda)*(*min_half_aQa));

    log_A += (*log_proposal_a - (*log_proposal_prop_a));
    if (!R_finite(log_A)) log_A = _AK_EMIN;
    //Rprintf("\nlog_proposal_a = %g,  log_proposal_prop_a = %g", *log_proposal_a, *log_proposal_prop_a);
    //Rprintf("\nLog(A)=%g,  A=%g\n", log_A, exp_AK(log_A));
  
    /** Accept/reject the proposal **/
    if (log_A < _AK_EMIN){
      *accept = 0;
      return;
    }
    if (log_A >= 0){
      *accept = 1;
    }
    else{  /* decide by sampling from exponential distribution */
      *ll = exp_rand();    
      *accept = (*ll > -log_A ? 1 : 0);
    }

    if (*accept){
      *lambda = *prop_lambda;
      GMRF_Gspline_Util::copy_within_update(expa,      sumexpa,      w,      minw,      Da,      Qa,      min_half_aQa,      a,
                                            prop_expa, prop_sumexpa, prop_w, prop_minw, prop_Da, prop_Qa, prop_min_half_aQa, prop_a,
                                            constraint,  iref,  na, &na_1,  GMRF_Gspline_Util::_a2a_);
    }
  }
  else{
    log_A += *lambda*((*prop_min_half_aQa) - (*min_half_aQa));

    log_A += (*log_proposal_a - (*log_proposal_prop_a));
    if (!R_finite(log_A)) log_A = _AK_EMIN;
    //Rprintf("\nlog_proposal_a = %g,  log_proposal_prop_a = %g", *log_proposal_a, *log_proposal_prop_a);
    //Rprintf("\nLog(A)=%g,  A=%g\n", log_A, exp_AK(log_A));
  
    /** Accept/reject the proposal **/
    if (log_A < _AK_EMIN){
      *accept = 0;
      return;
    }
    if (log_A >= 0){
      *accept = 1;
    }
    else{  /* decide by sampling from exponential distribution */
      *ll = exp_rand();    
      *accept = (*ll > -log_A ? 1 : 0);
    }

    if (*accept){
      GMRF_Gspline_Util::copy_within_update(expa,      sumexpa,      w,      minw,      Da,      Qa,      min_half_aQa,      a,
                                            prop_expa, prop_sumexpa, prop_w, prop_minw, prop_Da, prop_Qa, prop_min_half_aQa, prop_a,
                                            constraint,  iref,  na, &na_1,  GMRF_Gspline_Util::_a2a_);
    }
  }

  return;
}


/*** ================================================================================== ***/
/*** Maximum-likelihood, GMRF + G-spline                                                ***/
/***  Find the mode of the log-likelihood of the GMRF given allocations of the G-spline ***/
/***  and the precision of the GMRF                                                     ***/
/***                                                                                    ***/
/*** ================================================================================== ***/
//
// INPUT:
// ======
// ll[1]:            value of log-likelihood in the initial a(d)
// dll[na-1]:        first derivative of the log-likelihood in the initial d
// ddll[LT(na-1)]:   minus Hessian of the log-likelihood in the initial d
// a[na]:            initial value of a(d)
//
// OUTPUT:
// =======
// ll[1]:            value of log-likelihood in the mode
// dll[na-1]:        first derivative of the log-likelihood in the mode
// ddll[LT(na-1)]:   minus Hessian of the log-likelihood in the mode
// a[na]:            mode
//               
// workD[na+1+na+1+na+na+1+(na-1)+LT(na-1)+(na+LT(na))]:  working array (see section "Set-up pointers" below)
//               on INPUT, components of work except NR_step, ddll_temp, workll2  are assumed to be filled correctly with respect to a
// err[1]:      error flag:   0 = everything O.K.
//                            1 = maximum number of iterations was used without convergence
//                            2 = maximum number of step-halving steps was performed without increase of the objective function
//                            3 = minus second derivative matrix not positive definite and not possible to make it PD
//                                by inflation of the diagonal
//                            4 = NaN appeared in initial ll(x)
//
void
ML_est(double *ll,             double *dll,           double *ddll,
       double *a,              double *workD,         int *niter,          int *err,                 
       const int *allocN,      const double *lambda,  const double *Q,     const int *order,    const int *diffOper,
       const int *constraint,  const int *iref,
       const int *na,          const int *nobs,       const int *maxiter,  const double *epsw)
{
  *err = 0;

  //static const double HALF=0.5;

  static int na_1, info[1], Attempt[1], halfstep;
  static double old_ll, relat_diff;
  //static double *aP;  
  na_1 = *na - 1;


  /** Set-up pointers **/
  static double *expa, *sumexpa, *w, *minw, *Da, *Qa, *min_half_aQa, *NR_step, *ddll_temp, *workll2;
  expa         = workD;
  sumexpa      = expa         + (*na);
  w            = sumexpa      + 1;
  minw         = w            + (*na);
  Da           = minw         + 1;
  Qa           = Da           + (*na);
  min_half_aQa = Qa           + (*na);
  NR_step      = min_half_aQa + 1;
  ddll_temp    = NR_step      + na_1;
  workll2      = ddll_temp    + na_1*(*na)/2;

  /** Check initials for NaN **/
  if (!R_finite(*ll)){
    REprintf("Trap in GMRF_Gspline.cpp: ML_est. Initial values lead to the log-likelihood of -Inf.\n");
    *err = 4;
    return;
  }

  //Rprintf("\nallocN (nobs=%d):  ", *nobs);
  //AK_BLAS_LAPACK::printIArray(allocN, *na);
  //REprintf("\na(initial): ");
  //AK_BLAS_LAPACK::printArray(a, *na);
  //REprintf("\nw(initial): ");
  //AK_BLAS_LAPACK::printArray(w, *na);

  /*** Iterate (Newton-Raphson) ***/
  for (*niter = 0; *niter < *maxiter; (*niter)++){

    //REprintf("\nScore(NR iteration %d):  ", *niter);
    //AK_BLAS_LAPACK::printArray(dll, na_1);
    //REprintf("\nMinus Hessian (ll=%g,  NR iteration %d):\n", *ll, *niter);
    //AK_BLAS_LAPACK::printLT4R(ddll, na_1);

    /*** Cholesky decomposition of the Hessian, try to make it positive definite if it is not ***/
    AK_BLAS_LAPACK::chol_dpptrfPD(ddll, ddll_temp, &na_1, Attempt, &GMRF_Gspline::_max_nAttempt, epsw, info);
    if (*info){
      //REprintf("\nMinus Hessian (NR iteration %d):\n", *niter);
      //AK_BLAS_LAPACK::printLT(ddll, na_1);
      REprintf("Trap in GMRF_Gspline.cpp: ML_est. Singular Hessian.\n");
      *err = 3;
      return;
    }

    /*** Newton-Raphson step ***/
    AK_BLAS_LAPACK::copyArray(NR_step, dll, na_1);
    AK_BLAS_LAPACK::chol_solve_system(NR_step, ddll, &na_1, &_AK_ONE_INT);

    /*** Compute a new value of a ***/
    GMRF_Gspline_Util::NRstep(a, NR_step, constraint, iref, na);

    /*** Update derivatives ***/
    old_ll = *ll;
    GMRF_Gspline_Util::update4_ll12(expa, sumexpa, Da, min_half_aQa, Qa, w, minw, a, order, diffOper, na);
    GMRF_Gspline::ll2(ll, dll, ddll, workll2, a, allocN, lambda, sumexpa, min_half_aQa, Qa, w, Q, order, constraint, iref, na, nobs);    

    /*** Check convergence ***/
    relat_diff = R_finite(*ll) ? fabs(1 - old_ll/(*ll)) : R_PosInf;
    //Rprintf("NR-iteration %d,  old-ll=%g,  new-ll=%g,  relat-diff=%g\n", *niter, old_ll, *ll, relat_diff);
    //Rprintf("NR-step:  ");
    //AK_BLAS_LAPACK::printArray(NR_step, na_1);
    //Rprintf("New a:  ");
    //AK_BLAS_LAPACK::printArray(a, *na);
    //Rprintf("\n");
    if (relat_diff <= GMRF_Gspline::_toler){
      break;
    }

    /*** If not yet convergence, check whether the objective function increases ***/
    /***   if no increase perform step-halving                                  ***/
    if (!R_finite(*ll) || *ll < old_ll){      
      for (halfstep = 0; halfstep < GMRF_Gspline::_max_stephalf; halfstep++){
	GMRF_Gspline_Util::NRstephalve(a, NR_step, constraint, iref, na);
	GMRF_Gspline_Util::update4_ll0(expa, sumexpa, Da, min_half_aQa, a, order, na);
        GMRF_Gspline::ll0(ll, a, allocN, lambda, sumexpa, min_half_aQa, na, nobs);
        //Rprintf("Step-half %d, ll=%g\n", halfstep, *ll);
        if (*ll >= old_ll){
	  GMRF_Gspline_Util::update4_ll12(expa, sumexpa, Da, min_half_aQa, Qa, w, minw, a, order, diffOper, na);
          GMRF_Gspline::ll2(ll, dll, ddll, workll2, a, allocN, lambda, sumexpa, min_half_aQa, Qa, w, Q, order, constraint, iref, na, nobs);
          break;
        }
      }
      if (halfstep == GMRF_Gspline::_max_stephalf){
	GMRF_Gspline_Util::update4_ll12(expa, sumexpa, Da, min_half_aQa, Qa, w, minw, a, order, diffOper, na);
        GMRF_Gspline::ll2(ll, dll, ddll, workll2, a, allocN, lambda, sumexpa, min_half_aQa, Qa, w, Q, order, constraint, iref, na, nobs);
        *err = 2;
        break;
      }
    }
  }

  /*** Check niter ***/
  if (*maxiter && *niter == *maxiter) *err = 1;
  else                                (*niter)++;     /* to get the number of iterations really performed */

  return;
}


/*** ============================================================================= ***/
/*** Log-likelihood (except the additive constant)                                 ***/
/*** of the GMRF given allocations of the G-spline and the precision               ***/
/*** + its derivatives                                                             ***/
/***                                                                               ***/
/***   FUNCTIONS: ll0, ll1, ll2                                                    ***/
/***                                                                               ***/
/*** ============================================================================= ***/
//
// REMARK: All derivatives of the log-likelihood are with respect to the unconstrained parameter vector 'd'
//         where a = C*d, C = matrix na x (na-1)
//
//         It is assumed that d is equal to a with one omitted component
//         => after possible permutation of components in a, C' is equal to (eye(na-1), c) for some vector c
//
//         In this implementation, two types of identifiability constraints are considered:
//           _Mean_:        a[iref] = -sum(a[j], j<>iref)
//                                c = (-1, -1, ..., -1)'
//           _Reference_:   a[iref] = 0
//                                c = (0, 0, ..., 0)'
//
//  See document 'block-update-GMRF.pdf' for details.
//
//  ll[1]:              value of the log-likelihood
//  dll[na-1]:          gradient (w.r.t. d) of the log-likelihood
//  ddll[LT(na-1)]:     MINUS Hessian (w.r.t. d) of the log-likelihood (lower triangle stored in column major order)
//
//  a[na]:              values of the GMRF (ALL transformed mixture weights)
//  allocN[na]:         numbers of observations allocated in each mixture component 
//  lambda[1]:          value of the GMRF precision
//  sumexpa[1]:         sum(exp(a[j]))
//  min_half_aQa[1]:    value of -0.5*a'Q*a (penalty)
//  Qa[na]:             value of Q*a
//  w[na]:              weights = exp(a[j])/sum(exp(a[j]))
//  minw[1]:            value of the lowest weight
//  epsw[1]:            value which is added to the diagonal of the minus Hessian if (nobs*minw < epsw)
//  Q[LT(na)]:          lower triangle of the matrix Q = t(D)*D
//  order[1]:           order of the differences in the GMRF prior,
//                        that is, only diagonal and order subdiagonals of Q are non-zero
//
//  constraint[1]:      type of the identifiability constraint (see GMRF_Gspline_Util.h for possible values)
//  iref[1]:            index of a coefficient which is expressed as a function of d (remaining a coefficients)
//
//  na[1]:              dimension of the GMRF
//  nobs[1]:            number of observations upon which the allocN is based
//
//  workll2[]:          Working array for ll2_GMRF_Gspline function,
//                      array of length na + LT(na)
//

/*** Value of the log-likelihood ***/
void
ll0(double *ll,             
    const double *a,        const int *allocN,           const double *lambda,
    const double *sumexpa,  const double *min_half_aQa,  
    const int *na,          const int *nobs)
{
  static int j;
  static double aN;
  static const double *aP;
  static const int *NP;

  /*** Compute t(a)*N ***/
  aP = a;
  NP = allocN;
  aN = (*a)*(*NP);
  for (j = 1; j < *na; j++){
    aP++;
    NP++;
    aN += (*aP)*(*NP);
  }
 
  /*** Value of the log-likelihood ***/
  *ll = -(*nobs)*log_AK(*sumexpa) + (*lambda)*(*min_half_aQa) + aN;

  return;
}


/*** Value of the log-likelihood and the gradient ***/
void
ll1(double *ll,             double *dll,
    const double *a,        const int *allocN,           const double *lambda,
    const double *sumexpa,  const double *min_half_aQa,  const double *Qa,      const double *w,
    const int *constraint,  const int *iref,        
    const int *na,          const int *nobs)
{
  static int j;
  static double *dllP;
  static const double *wP, *QaP, *wRef, *QaRef;
  static const int *NP, *NRef;

  /*** Value of the log-likelihood ***/
  GMRF_Gspline::ll0(ll, a, allocN, lambda, sumexpa, min_half_aQa, na, nobs);

  /*** Gradient of the log-likelihood ***/
  dllP = dll;
  wP   = w;
  QaP  = Qa;
  NP   = allocN;
  j    = 0;
  switch (*constraint){  
  case GMRF_Gspline_Util::_Mean_:
    wRef  = w + (*iref);
    QaRef = Qa + (*iref);
    NRef  = allocN + (*iref);
    while (j < *iref){
      *dllP = -(*nobs)*(*wP - (*wRef)) - (*lambda)*(*QaP - (*QaRef)) + (*NP - (*NRef));
      dllP++;
      wP++;
      QaP++;
      NP++;
      j++;
    }
    wP++;
    QaP++;
    NP++;
    j++;
    while (j < *na){
      *dllP = -(*nobs)*(*wP - (*wRef)) - (*lambda)*(*QaP - (*QaRef)) + (*NP - (*NRef));
      dllP++;
      wP++;
      QaP++;
      NP++;
      j++;
    }
    return;

  case GMRF_Gspline_Util::_Reference_:
    while (j < *iref){
      *dllP = -(*nobs)*(*wP) - (*lambda)*(*QaP) + (*NP);
      dllP++;
      wP++;
      QaP++;
      NP++;
      j++;
    }
    wP++;
    QaP++;
    NP++;
    j++;
    while (j < *na){
      *dllP = -(*nobs)*(*wP) - (*lambda)*(*QaP) + (*NP);
      dllP++;
      wP++;
      QaP++;
      NP++;
      j++;
    }
    return;

  default:
    REprintf("constraint = %d\n", *constraint);
    throw returnR("Error in GMRF_Gspline.cpp: ll1(). Unknown value of the identifiability constraint.", 1);
  }
}


/*** Value of the log-likelihood, the gradient and the minus Hessian ***/
void
ll2(double *ll,             double *dll,                 double *ddll,          double *workll2,
    const double *a,        const int *allocN,           const double *lambda,
    const double *sumexpa,  const double *min_half_aQa,  const double *Qa,      const double *w,
    const double *Q,        const int *order,
    const int *constraint,  const int *iref,
    const int *na,          const int *nobs)
{
  static int j, i, iband;
  static double *ddllP;
  static const double *QP, *wP, *wP2;
  static double *m0, *M;

  /*** Set-up pointers for working arrays ***/
  m0 = workll2;
  M  = m0 + (*na);

  /*** Value of the log-likelihood and the gradient ***/
  GMRF_Gspline::ll1(ll, dll, a, allocN, lambda, sumexpa, min_half_aQa, Qa, w, constraint, iref, na, nobs);

  /*** Minus Hessian:  Compute lambda*Q, store it in M. ***/
  ddllP = M;
  QP    = Q;
  for (j = 0; j < *na; j++){
    iband = j + (*order) + 1;        /* index of the first zero row under the diagonal */
    if (iband > *na) iband = *na;
  
    /** non-zero elements **/
    for (i = j; i < iband; i++){
      *ddllP = (*lambda)*(*QP);
      ddllP++;
      QP++;
      }
  
    /** zero elements **/
    for (i = iband; i < *na; i++){
     *ddllP = 0.0;
      ddllP++;
      QP++;
    }
  }

  /*** Minus Hessian: Compute lambda*t(C)*Q*C, store it in ddll. ***/
  switch (*constraint){
  case GMRF_Gspline_Util::_Mean_:
    AK_BLAS_LAPACK::ALT_BLT_min1b_minb1_plusb(ddll, m0, M, *na, *iref);
    break;

  case GMRF_Gspline_Util::_Reference_:
    AK_BLAS_LAPACK::ALT_BLTremoveRowCol(ddll, m0, M, *na, *iref);
    break;

  default:
    REprintf("constraint = %d\n", *constraint);
    throw returnR("Error in GMRF_Gspline.cpp: ll2(). Unknown value of the identifiability constraint.", 1);
  }


  /*** Minus Hessian: Compute n*(diag(w) - w*t(w)), store it in M. ***/
  wP    = w;
  ddllP = M;
  for (j = 0; j < *na; j++){
    wP2 = wP;

    /** diagonal **/
    *ddllP = (*nobs)*((*wP) - (*wP)*(*wP2));
    ddllP++;
    wP2++;

    /** sub-diagonal **/
    for (i = j+1; i < *na; i++){
      *ddllP = -(*nobs)*(*wP)*(*wP2);
      ddllP++;
      wP2++;
    }    
    wP++;
  }

  /*** Minus Hessian: Compute n*t(C)*(diag(w) - w*t(w))*C, add it to ddll. ***/
  switch (*constraint){
  case GMRF_Gspline_Util::_Mean_:
    AK_BLAS_LAPACK::ALT_pp_BLT_min1b_minb1_plusb(ddll, m0, M, *na, *iref);
    return;

  case GMRF_Gspline_Util::_Reference_:
    AK_BLAS_LAPACK::ALT_pp_BLTremoveRowCol(ddll, m0, M, *na, *iref);
    return;

  default:
    REprintf("constraint = %d\n", *constraint);
    throw returnR("Error in GMRF_Gspline.cpp: ll2(). Unknown value of the identifiability constraint.", 1);
  }
}


/*** ================================================================================== ***/
/*** MCMC to test 'GMRF_Gspline::update' function                                       ***/
/***                                                                                    ***/
/***  * to be called from R                                                             ***/
/***  * sampled values are returned directly to R                                       ***/
/***    (they are not written to files)                                                 ***/
/***                                                                                    ***/
/*** ================================================================================== ***/
//
//  accept[niter]
//  aSample[na x niter]
//  wSample[na x niter]
//  lambdaSample[niter]
//  iter[1]:               Number of the null-th iteration
//  allocN[na]
//  par_lambda[2]
//  F[1]:                   Factor F>1 for the lambda-proposal
//  order[1]
//  A[1 x na]
//  na[1]
//  nobs[1]
//  lambda_a_block[1]
//  nsimul[3]
//  
//
extern "C"{
  void
  mcmc_GMRF_Gspline(int *acceptSample,   double *aSample,              double *wSample,   double *lambdaSample,  
                    int *iter,
                    const int *allocN,   const int *prior_for_lambda,  const double *par_lambda,  const double *F,  
                    const int *order,    const int *constraint,        const int *iref,
                    const int *na,       const int *nobs,              const int *lambda_a_block,
                    const int *nsimul)
  {
    try{
      int i;

      GetRNGstate(); 

      /*** Numbers of iterations etc. ***/
      const int niter = nsimul[0];
      const int nthin = nsimul[1];
      const int nwrite = nsimul[2];

      const double epsw[1] = {GMRF_Gspline::_null_mass};
      const int nc[1] = {0};
      const int LTnc = (*nc)*(*nc+1)/2;
      const int naTnc = (*na)*(*nc);
      const int na_1 = *na - 1;
      //const int na_order = *na - (*order);
      const int LTna = (*na)*(*na+1)/2;
      const int LTna_1 = na_1*(*na)/2;

      /*** Working arrays ***/
      const int nworkML   = (*na) + 1 + (*na) + 1 + 2*(*na) + 1 + na_1 + LTna_1 + (*na) + LTna;
      const int nworka    = (*na) + LTna_1 + na_1 + 2*(*na) + 1 + (*na) + 1 + 2*(*na) + 1;

      const int nworksGMRF2[5] = {
        (*nc)*(*nc),                   /** log_density_Ax_x **/
        LTnc + naTnc + (*nc),          /** rGMRF_inputArgs  **/
        (*na > *nc ? *na : *nc),       /** rGMRF            **/
        *nc,                           /** dGMRF_inputArgs  **/ 
        *na                            /** dGMRF            **/
      };
      const int nworkGMRF2 = maxArray(nworksGMRF2, 5);
      const int nworkGMRF = 1 + 4 + naTnc + LTnc + nworkGMRF2;

      double *workML   = (double*) calloc(nworkML, sizeof(double));
      double *worka    = (double*) calloc(nworka, sizeof(double));
      double *workGMRF = (double*) calloc(nworkGMRF, sizeof(double));
      if (!workML || !worka || !workGMRF) throw returnR("Out of memory in GMRF_Gspline.cpp: mcmc_GMRF_Gspline().", 99);

      /*** Difference operator etc. ***/
      int *diffOper = (int*) calloc(*order + 1, sizeof(int));
      double *Q     = (double*) calloc(LTna, sizeof(double));      
      if (!diffOper || !Q) throw returnR("Out of memory in GMRF_Gspline.cpp: mcmc_GMRF_Gspline().", 99);
      GMRF::diff_operator(diffOper, order);
      GMRF::Q_matrix(Q, order, na);

      /*** Parameters of the lambda proposal ***/
      double par_rscale[6];
      GMRF::dscale_norm_const(F, par_rscale);

      /*** Current values of parameters, related quantities and their initialization ***/
      double lambda[1] = {*lambdaSample};

      double *a     = (double*) calloc(*na, sizeof(double));      
      double *expa  = (double*) calloc(*na, sizeof(double));
      double sumexpa[1];
      double *w     = (double*) calloc(*na, sizeof(double));
      double minw[1];
      double *Da    = (double*) calloc(*na, sizeof(double));            
      double *Qa    = (double*) calloc(*na, sizeof(double));            
      double min_half_aQa[1];
      if (!a || !expa || !w || !Da || !Qa) throw returnR("Out of memory in GMRF_Gspline.cpp: mcmc_GMRF_Gspline().", 99);
      AK_BLAS_LAPACK::copyArray(a, aSample, *na);
      GMRF_Gspline_Util::update4_ll12(expa, sumexpa, Da, min_half_aQa, Qa, w, minw, a, order, diffOper, na);


      /*** MCMC ***/
      const double *aP, *wP;

      int accept[1];

      int *acceptSampleP    = acceptSample;
      double *aSampleP      = aSample;
      double *wSampleP      = wSample;
      double *lambdaSampleP = lambdaSample;

      int nullthIter   = *iter;
      int lastIter     = nullthIter + niter;
      int iterTotal    = nullthIter*nthin;
      int iterTotalNow = 0;
      int backs        = 0;
      int writeAll     = 0;
      int witer;
  
      Rprintf("Iteration ");    
      for (*iter=nullthIter + 1; *iter <= lastIter; (*iter)++){
        for (witer = 1; witer <= nthin; witer++){               /**  thinning cycle  **/
          iterTotal++;                                          /* = (iter-1)*nthin + witer                    */
          iterTotalNow++;                                       /* = (iter-1)*nthin + witer - nullthIter*nthin */

          GMRF_Gspline::update(accept, a, lambda, expa, sumexpa, w, minw, Da, Qa, min_half_aQa, workML, worka, workGMRF,
                               allocN, prior_for_lambda, par_lambda, par_rscale, Q, order, diffOper, 
                               epsw, constraint, iref, na, nobs, lambda_a_block);
          //goto CLEANING;
        }

        *acceptSampleP = *accept;
        acceptSampleP++;

        *lambdaSampleP = *lambda;
        lambdaSampleP++;

        aP = a;
        wP = w;
        for (i = 0; i < *na; i++){
          *aSampleP = *aP;
          *wSampleP = *wP;
          aP++;
          wP++;
          aSampleP++;
          wSampleP++;
        }                        

        if (!(*iter % nwrite) || *iter == lastIter){
          writeAll = 1;
          for (i = 0; i < backs; i++) Rprintf("\b");
            Rprintf("%d", *iter);
            backs = int(log10(double(*iter))) + 1;
        }
      }


      /*** Cleaning ***/
      //CLEANING:
      Rprintf("\n");
      free(Qa);
      free(Da);
      free(w);
      free(expa);
      free(a);
      free(Q);
      free(diffOper);
      free(workGMRF);
      free(worka);
      free(workML);

      PutRNGstate();
      return;
    }
    catch (returnR rr){
      PutRNGstate();
      return;
    }
  }
}

}  /*** end of namespace GMRF_Gspline ***/
