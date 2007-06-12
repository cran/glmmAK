/*** mcmc_common.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  18/09/2006
//        
//                    update_reg_gamermanCL:  18/09/2006
//                 update_reg_gamermanPoiss:  13/02/2007
//                        update_norm_mean1:  19/09/2006
//               update_norm_invvar_wishart:  21/09/2006
//             update_norm_invvar_sduniform:  27/10/2006
//
// PURPOSE: MCMC related functions useful for several purposes
//
/* ********************************************************************************* */

#include "mcmc_common.h"

namespace mcmc_common {

/* ----------------------------------------------------------------------------------------------------------------------------- */
/*                                                                                                                               */
/* update_reg_gamermanCL: Update of regression coefficients or random effects                                                    */
/*                      in the cumulative logit GALMM model using the Metropolis-Hastings algorithm                              */
/*                      with the proposal obtained by performing 1 Fisher scoring or                                             */
/*                      Newton-Raphson step (see Gamerman, 1997. Statistics and Computing, 7, pp. 57-68.                         */
/*                                                                                                                               */
/* ----------------------------------------------------------------------------------------------------------------------------- */
/*                                                                                                                               */
/*        accept:  indicator of acceptance (1 if accepted, 0 if not accepted)                                                    */
/*            ll:  value of the log-likelihood (it is updated)                                                                   */
/*             U:  working vector of length nTheta                                                                               */
/*             I:  working symmetric matrix of dimension nTheta x nTheta                                                         */
/*          etaX:  vector of length nObs, linear predictor proportional w.r.t. odds                                              */
/*      PropetaX:  working vector of length nObs                                                                                 */
/*          etaV:  matrix C x nObs, linear predictor not proportional w.r.t. odds                                                */
/*      PropetaV:  working vector of length C*nObs                                                                               */
/*           eta:  matrix C x nObs, linear predictor eta = etaX + etaV                                                           */
/*       Propeta:  working vector of length C*nObs                                                                               */
/*          prob:  matrix (C+1) x nObs, category probabilities for each observation                                              */
/*      Propprob:  working vector of length (C+1)*nObs                                                                           */
/*          work:  working vector for ll_* functions, length = 5*C                                                               */
/*        offset:  matrix C x nObs, offset to the linear predictor (may depend on c)                                             */
/*         Theta:  vector of length p + C*q, regression coefficients                                                             */
/*     PropTheta:  working vector of length p + C*q                                                                              */
/*             y:  vector of length nObs, response                                                                               */
/*             x:  matrix p x n, covariates proportional w.r.t. odds                                                             */
/*             v:  matrix q x n, covariates not proportional w.r.t. odds                                                         */
/*            xx:  array of length nObs, each element is a matrix holding x*t(x)                                                 */
/*            vv:  array of length nObs, each element is a matrix holding v*t(v)                                                 */
/*            xv:  array of length nObs, each element is a matrix holding x*t(v)                                                 */
/*          nObs:  number of observations                                                                                        */
/*             C:  number of response categories minus 1                                                                         */
/*             p:  number of covariates whose effect is proportional w.r.t. odds                                                 */
/*             q:  number of covariates whose effect is not proportional w.r.t. odds                                             */
/*     PriorMean:  vector a, where a = prior mean of Theta                                                                       */
/*    PriorWMean:  vector R^{-1}*a, where R = prior variance of Theta and a = prior mean of Theta                                */
/*   PriorInvVar:  if diagInvVar = true: vector of length nTheta with diagonal of R^{-1}                                         */
/*                 if diagInvVar = false: lower triangle with the matrix R^{-1}                                                  */
/*  PriorInvVarL:  if diagInvVar = true: NULL                                                                                    */
/*                 if diagInvVar = false: Cholesky decomposition of PriorInvVar                                                  */
/*    diagInvVar:  true/false                                                                                                    */
/*      PropMean:  working vector of length nTheta                                                                               */
/*        ll_FS2:  function that computes the log-likelihood, score and EXPECTED information matrix,                             */
/*                 e.g., ll_cumlogitFS2 from 'll_cumlogit.cpp'                                                                   */
/*        ll_NR2:  function that computes the log-likelihood, score and OBSERVED information matrix,                             */
/*                 e.g., ll_cumlogitNR2 from 'll_cumlogit.cpp'                                                                   */
/*        caller:  name of the main function that calls 'update_reg1' to provide informative error messages                      */
/*                                                                                                                               */
/* ----------------------------------------------------------------------------------------------------------------------------- */
void
update_reg_gamermanCL(int *accept,                       
                      double *ll,
                      MatrixRect<double> *U,
                      MatrixLT<double> *I,
                      double *etaX,                      double *PropetaX,
                      double *etaV,                      double *PropetaV,
                      double *eta,                       double *Propeta,
                      double *prob,                      double *Propprob,
                      double *work,
                      const double *offset,              
                      double *Theta,                     double *PropTheta,
                      const int *y,
                      const double *x,
                      const double *v,
                      const MatrixLT<double> *xx,
                      const MatrixLT<double> *vv,
                      const MatrixRect<double> *xv,
                      const int &nObs,
                      const int &C,
                      const int &p,
                      const int &q,
                      const double *PriorMean,
                      const double *PriorWMean,
                      const double *PriorInvVar,
                      const MatrixLT<double> *PriorInvVarL,
                      const bool &diagInvVar,
                      MatrixRect<double> *PropMean,
                      void (*ll_FS2)(double*, MatrixRect<double>*, MatrixLT<double>*, 
                                     double*, double*, double*, double*, double*, int*,
                                     const double*, const double*, const int*, const double*, const double*,
                                     const MatrixLT<double>*, const MatrixLT<double>*, const MatrixRect<double>*,
                                     const int&, const int&, const int&, const int&, const int&, const int&),
                      void (*ll_NR2)(double*, MatrixRect<double>*, MatrixLT<double>*, 
                                     double*, double*, double*, double*, double*, int*,
                                     const double*, const double*, const int*, const double*, const double*,
                                     const MatrixLT<double>*, const MatrixLT<double>*, const MatrixRect<double>*,
                                     const int&, const int&, const int&, const int&, const int&, const int&),
                      const char *caller)
{
  static int i, j;
  static int anyZero, Rank, Attempt;
  static double prop_ll, prop_ldprop, prop_lprior, ldprop, lprior, prop_lratio, erand;
  static double *probOld, *probNew, *etaXOld, *etaXNew, *etaVOld, *etaVNew, *etaOld, *etaNew, *valOld, *valNew;


  /*** Fisher-scoring step ***/
  /* If it is not possible to perform Fisher-scoring step due to zero probabilities of some categories, try Newton-Raphson step.   */
  /*== log-likelihood(theta) ==*/
  ll_FS2(ll, U, I, etaX, etaV, eta, prob, work, &anyZero, offset, Theta, y, x, v, xx, vv, xv, nObs, C, p, q, 0, 2);
  if (!R_finite(*ll)){
    throw returnR(caller, "Trap in update_reg_gamermanCL. Infinite log-likelihood encountered.", 1);
  }
  if (anyZero){
    ll_NR2(ll, U, I, etaX, etaV, eta, prob, work, &anyZero, offset, Theta, y, x, v, xx, vv, xv, nObs, C, p, q, 1, 2);
  }

  /*** Compute (partial) moments of the normal proposal            ***/
  /* * (partial) mean of the proposal will be stored in _PropMean    */
  /* * inverse variance of the proposal will be stored in _I         */
  Ab2(PropMean->a(), I, Theta);              /* I*Theta                */
  PropMean->add_A(U);                        /* I*Theta + U            */
  PropMean->add_A2(PriorWMean);              /* R^{-1}*a + I*Theta + U */

  if (diagInvVar) I->add_b2diag(PriorInvVar);           /* R^{-1} + I             */
  else            I->add_A2(PriorInvVar);

  /*** Propose a new value ***/  
  U->rNorm(0, 1);
  Rank = I->chol_solvePD(Attempt, PropMean->a(), PropTheta, mcmc_common::_AK_MAX_PD_ATTEMPTS_2, mcmc_common::_AK_EPS_PD_ATTEMPT_2, 1, 1, U->a());
         /* NOW: PropTheta = the proposed value                                      */
         /*              U = (PropTheta - mean of proposal)                          */ 
         /*              I = Cholesky factor of the inverse variance of the proposal */

  if (Rank < I->nrow()){
    throw returnR(caller, "Trap in update_reg_gamermanCL. Not possible to make the proposal variance matrix positive definite.", 1);
  }
  
  /*** Compute the logarithm of the proposal ratio ***/
    /*== log-q(theta, theta[proposed]) ==*/
  prop_ldprop = Mvtdist2::dmvtnorm2(*U, *I, 1, 0);  

    /*== log-prior(theta) and log-prior(theta[proposed]) ==*/
  if (diagInvVar){
    lprior      = Mvtdist2::dmvtnorm2b(Theta, PriorMean, PriorInvVar, I->nrow(), 1);
    prop_lprior = Mvtdist2::dmvtnorm2b(PropTheta, PriorMean, PriorInvVar, I->nrow(), 1);
  }
  else{
    lprior      = Mvtdist2::dmvtnorm2a(*U, Theta, PriorMean, PriorInvVarL, 1);
    prop_lprior = Mvtdist2::dmvtnorm2a(*U, PropTheta, PriorMean, PriorInvVarL, 1);
  }

    /*== log-likelihood(theta[proposed]) ==*/
  ll_FS2(&prop_ll, U, I, PropetaX, PropetaV, Propeta, Propprob, work, &anyZero, offset, PropTheta, y, x, v, xx, vv, xv, nObs, C, p, q, 1, 2);
  if (!R_finite(prop_ll)){
    *accept = 0;
    return;
  }
  if (anyZero){
    ll_NR2(&prop_ll, U, I, PropetaX, PropetaV, Propeta, Propprob, work, &anyZero, offset, PropTheta, y, x, v, xx, vv, xv, nObs, C, p, q, 1, 2);
  }

    /*== log-q(theta[proposed], theta) ==*/    
  Ab2(PropMean->a(), I, PropTheta);
  PropMean->add_A(U);
  PropMean->add_A2(PriorWMean);
  
  if (diagInvVar) I->add_b2diag(PriorInvVar);           /* R^{-1} + I             */
  else            I->add_A2(PriorInvVar);

  Rank = I->chol_solvePD(Attempt, PropMean->a(), U->a(), mcmc_common::_AK_MAX_PD_ATTEMPTS_2, mcmc_common::_AK_EPS_PD_ATTEMPT_2, 0, 0, NULL);
         /* NOW: U = the mean of the reversed proposal                                */
         /*      I = Cholesky factor of the inverse variance of the reversed proposal */

  if (Rank < I->nrow()){
    *accept = 0;
    return;
  }

  U->subtract_A2(Theta);                  /* NOW: U = (mean of reversed proposal - Theta) */
  ldprop = Mvtdist2::dmvtnorm2(*U, *I, 1, 0);  

    /*=== logarithm of the proposal ratio ===*/
  prop_lratio = prop_ll + prop_lprior + ldprop - (*ll) - lprior - prop_ldprop;
  if (prop_lratio < _AK_EMIN){
    *accept = 0;
    return;
  }
  if (prop_lratio >= 0){
    *accept = 1;
  }
  else{  /* decide by sampling from exponential distribution */
    erand = exp_rand();    
    *accept = (erand > -prop_lratio ? 1 : 0);
  }

  /*** Make the proposed value the new value (if accepted) ***/
  if (*accept){
    probOld = prob;
    probNew = Propprob;
    etaXOld = etaX;
    etaXNew = PropetaX;
    etaVOld = etaV;
    etaVNew = PropetaV;
    etaOld = eta;
    etaNew = Propeta;
    for (i = 0; i < nObs; i++){
      *etaXOld = *etaXNew;
      etaXOld++;
      etaXNew++;
      *probOld = *probNew;
      probOld++;
      probNew++;
      for (j = 0; j < C; j++){
        *probOld = *probNew;
        probOld++;
        probNew++;
        *etaVOld = *etaVNew;
        etaVOld++;
        etaVNew++;
        *etaOld = *etaNew;
        etaOld++;
        etaNew++;
      }
    }
    valOld = Theta;
    valNew = PropTheta;
    for (i = 0; i < p + C*q; i++){
      *valOld = *valNew;
      valOld++;
      valNew++;
    }
    *ll = prop_ll;
  }

  return;  
}


/* ----------------------------------------------------------------------------------------------------------------------------- */
/*                                                                                                                               */
/* update_reg_gamermanPoiss: Update of regression coefficients or random effects                                                 */
/*                      in the Poisson log-linear GALMM model using the Metropolis-Hastings algorithm                            */
/*                      with the proposal obtained by performing 1 Fisher scoring or                                             */
/*                      Newton-Raphson step (see Gamerman, 1997. Statistics and Computing, 7, pp. 57-68.                         */
/*                                                                                                                               */
/* ----------------------------------------------------------------------------------------------------------------------------- */
/*                                                                                                                               */
/*        accept:  indicator of acceptance (1 if accepted, 0 if not accepted)                                                    */
/*            ll:  value of the log-likelihood (it is updated)                                                                   */
/*             U:  working vector of length nTheta                                                                               */
/*             I:  working symmetric matrix of dimension nTheta x nTheta                                                         */
/*           eta:  vector of length nObs, linear predictor eta                                                                   */
/*       Propeta:  working vector of length nObs                                                                                 */
/*            mu:  vector of length nObs, expected counts                                                                        */
/*        Propmu:  working vector of length nObs                                                                                 */
/*          work:  working vector for ll_* functions, length = LT(p)                                                             */
/*        offset:  vector of length nObs, offset to the linear predictor                                                         */
/*         Theta:  vector of length p, regression coefficients                                                                   */
/*     PropTheta:  working vector of length p                                                                                    */
/*             y:  vector of length nObs, response                                                                               */
/*    log_y_factor:                                                                                                              */
/*             x:  matrix p x n, covariates proportional w.r.t. odds                                                             */
/*            xx:  array of length nObs, each element is a matrix holding x*t(x)                                                 */
/*          nObs:  number of observations                                                                                        */
/*             p:  number of covariates (length of Theta)                                                                        */
/*     PriorMean:  vector a, where a = prior mean of Theta                                                                       */
/*    PriorWMean:  vector R^{-1}*a, where R = prior variance of Theta and a = prior mean of Theta                                */
/*   PriorInvVar:  if diagInvVar = true: vector of length nTheta with diagonal of R^{-1}                                         */
/*                 if diagInvVar = false: lower triangle with the matrix R^{-1}                                                  */
/*  PriorInvVarL:  if diagInvVar = true: NULL                                                                                    */
/*                 if diagInvVar = false: Cholesky decomposition of PriorInvVar                                                  */
/*                                        (lower triangle)                                                                       */
/*    diagInvVar:  true/false                                                                                                    */
/*      PropMean:  working vector of length nTheta                                                                               */
/*        caller:  name of the main function that calls 'update_reg1' to provide informative error messages                      */
/*                                                                                                                               */
/* ----------------------------------------------------------------------------------------------------------------------------- */
void
update_reg_gamermanPoiss(int *accept,                       
                         double *ll,
                         double *U,
                         double *I,
                         double *eta,                       double *Propeta,
                         double *mu,                        double *Propmu,
                         double *work,
                         const double *offset,              
                         double *Theta,                     double *PropTheta,
                         const int *y,
                         const double *log_y_factor,
                         const double *x,
                         const MatrixLT<double> *xx,
                         const int &nObs,
                         const int &p,
                         const double *PriorMean,
                         const double *PriorWMean,
                         const double *PriorInvVar,
                         const double *PriorInvVarL,
                         const bool &diagInvVar,
                         double *PropMean,
                         const char *caller)
{
  static int i;
  static double prop_ll, prop_ldprop, prop_lprior, ldprop, lprior, prop_lratio, erand;
  static double *muOld, *muNew, *etaOld, *etaNew, *valOld, *valNew;
  static int LTp;
  static int Attempt, info;

  LTp = (p*(p+1))/2;

  //Rprintf("\n================================\n");
  //Rprintf("PriorMean=");
  //AK_BLAS_LAPACK::printArray(PriorMean, p);
  //Rprintf("PriorInvVar=");
  //if (diagInvVar) AK_BLAS_LAPACK::printArray(PriorInvVar, p);
  //else            AK_BLAS_LAPACK::printLT4R(PriorInvVar, p);

  /*** Fisher-scoring step ***/
  /*== log-likelihood(theta) ==*/
  Ll_poisson::ll_poisson(ll, U, I, eta, mu, offset, Theta, y, log_y_factor, x, xx, nObs, p, 2);
  if (!R_finite(*ll)) throw returnR(caller, "Trap in update_reg_gamermanPoiss. Infinite log-likelihood encountered.", 1);
  //Rprintf("Theta=");
  //AK_BLAS_LAPACK::printArray(Theta, p);
  //Rprintf("ll=%g,  score=", *ll);
  //AK_BLAS_LAPACK::printArray(U, p);
  //Rprintf("I=");
  //AK_BLAS_LAPACK::printLT4R(I, p);

  /*** Compute (partial) moments of the normal proposal            ***/
  /* * (partial) mean of the proposal will be stored in _PropMean    */
  /* * inverse variance of the proposal will be stored in _I         */
  AK_BLAS_LAPACK::c_ALTb(PropMean, I, Theta, &p);            /* PropMean = I*Theta                 */
  AK_BLAS_LAPACK::a_aPlusb(PropMean, U, p);                  /* PropMean = I*Theta + U             */
  AK_BLAS_LAPACK::a_aPlusb(PropMean, PriorWMean, p);         /* PropMean = R^{-1}*a + I*Theta + U  */
  //Rprintf("PropMean (canonical)=");
  //AK_BLAS_LAPACK::printArray(PropMean, p);

  if (diagInvVar) AK_BLAS_LAPACK::ALT_addb2diag(I, PriorInvVar, p);        /* I = R^{-1} + I             */
  else            AK_BLAS_LAPACK::a_aPlusb(I, PriorInvVar, LTp);
  //Rprintf("PropPrecision=");
  //AK_BLAS_LAPACK::printLT4R(I, p);

  /*** Propose a new value ***/  
  AK_BLAS_LAPACK::chol_dpptrfPD(I, work, &p, &Attempt, &mcmc_common::_AK_MAX_PD_ATTEMPTS_2, &mcmc_common::_AK_EPS_PD_ATTEMPT_2, &info);
  //Rprintf("Cholesky decomposition of PropPrecision=");
  //AK_BLAS_LAPACK::printLT4R(I, p);
  if (info){
    Rprintf("\nLog-likelihood = %g  \nTheta = ", *ll);
    AK_BLAS_LAPACK::printArray(Theta, p);    
    throw returnR(caller, "Trap in update_reg_gamermanPoiss. Not possible to make the proposal variance matrix positive definite.", 1);
  }
  Mvtdist3::rmvnormC2006b(PropTheta, U, PropMean, I, &p);   /* NOW:  PropTheta = the proposed value                                                  */
                                                            /*        PropMean = mean of the proposal = (R^{-1} + I){-1} * (R^{-1}*a + I*Theta + U)  */
                                                            /*               U = t(L)*(PropTheta - mean of proposal)                                 */
                                                            /*               I = L = Cholesky factor of the inverse variance of the proposal         */ 
  //Rprintf("PropMean=");
  //AK_BLAS_LAPACK::printArray(PropMean, p);
  //Rprintf("PropTheta=");
  //AK_BLAS_LAPACK::printArray(PropTheta, p);

  /*** Compute the logarithm of the proposal ratio ***/
    /*== log-q(theta, theta[proposed]) ==*/
  Mvtdist3::ldmvnorm2006b(&prop_ldprop, U, I, &p);

    /*== log-prior(theta) and log-prior(theta[proposed]) ==*/
  if (diagInvVar){
    Mvtdist3::ldmvnorm2007b(&lprior, Theta, PriorMean, PriorInvVar, &p);
    Mvtdist3::ldmvnorm2007b(&prop_lprior, PropTheta, PriorMean, PriorInvVar, &p);
  }
  else{
    Mvtdist3::ldmvnorm2007a(&lprior, U, Theta, PriorMean, PriorInvVarL, &p);
    Mvtdist3::ldmvnorm2007a(&prop_lprior, U, PropTheta, PriorMean, PriorInvVarL, &p);
  }

    /*== log-likelihood(theta[proposed]) ==*/
  Ll_poisson::ll_poisson(&prop_ll, U, I, Propeta, Propmu, offset, PropTheta, y, log_y_factor, x, xx, nObs, p, 2);
  if (!R_finite(prop_ll)){
    *accept = 0;
    return;
  }
  //Rprintf("-------------------\n");
  //Rprintf("PropTheta=");
  //AK_BLAS_LAPACK::printArray(PropTheta, p);
  //Rprintf("ll(iprop)=%g,  score=", prop_ll);
  //AK_BLAS_LAPACK::printArray(U, p);
  //Rprintf("I(iprop)=");
  //AK_BLAS_LAPACK::printLT4R(I, p);

    /*== log-q(theta[proposed], theta) ==*/    
  AK_BLAS_LAPACK::c_ALTb(PropMean, I, PropTheta, &p);        /* PropMean = I*PropTheta                 */
  AK_BLAS_LAPACK::a_aPlusb(PropMean, U, p);                  /* PropMean = I*PropTheta + U             */
  AK_BLAS_LAPACK::a_aPlusb(PropMean, PriorWMean, p);         /* PropMean = R^{-1}*a + I*PropTheta + U  */
  //Rprintf("PropMean(iprop) (canonical)=");
  //AK_BLAS_LAPACK::printArray(PropMean, p);

  if (diagInvVar) AK_BLAS_LAPACK::ALT_addb2diag(I, PriorInvVar, p);        /* I = R^{-1} + I             */
  else            AK_BLAS_LAPACK::a_aPlusb(I, PriorInvVar, LTp);
  //Rprintf("PropPrecision(iprop)=");
  //AK_BLAS_LAPACK::printLT4R(I, p);

  AK_BLAS_LAPACK::chol_dpptrfPD(I, work, &p, &Attempt, &mcmc_common::_AK_MAX_PD_ATTEMPTS_2, &mcmc_common::_AK_EPS_PD_ATTEMPT_2, &info);
  //Rprintf("Cholesky decomposition of PropPrecision(iprop)=");
  //AK_BLAS_LAPACK::printLT4R(I, p);
  if (info){
    *accept = 0;
    return;
  }
  Mvtdist3::ldmvnormC2006(&ldprop, PropMean, Theta, I, &p);

    /*=== logarithm of the proposal ratio ===*/
  prop_lratio = prop_ll + prop_lprior + ldprop - (*ll) - lprior - prop_ldprop;
  //Rprintf("prop_lratio=%g;\n   prop_ll=%g, prop_lprior=%g, ldprop=%g;\n   ll=%g, lprior=%g, prop_ldprop=%g\n", prop_lratio, prop_ll, prop_lprior, ldprop, *ll, lprior, prop_ldprop);
  if (prop_lratio < _AK_EMIN){
    *accept = 0;
    return;
  }
  if (prop_lratio >= 0){
    *accept = 1;
  }
  else{  /* decide by sampling from exponential distribution */
    erand = exp_rand();    
    *accept = (erand > -prop_lratio ? 1 : 0);
  }

  /*** Make the proposed value the new value (if accepted) ***/
  if (*accept){
    muOld   = mu;
    muNew   = Propmu;
    etaOld = eta;
    etaNew = Propeta;
    for (i = 0; i < nObs; i++){
      *muOld = *muNew;
      muOld++;
      muNew++;
      *etaOld = *etaNew;
      etaOld++;
      etaNew++;
    }

    valOld = Theta;
    valNew = PropTheta;
    for (i = 0; i < p; i++){
      *valOld = *valNew;
      valOld++;
      valNew++;
    }
    *ll = prop_ll;
  }

  return;  
}


/* ----------------------------------------------------------------------------------------------------------------------------- */
/*                                                                                                                               */
/* update_norm_mean1: Update of means of normal effects where the mean itself has a conjugate normal prior                       */
/*                                                                                                                               */
/*    p(b[i] | alpha, Sigma) = N(alpha, Sigma), i = 1,...,N                                                                      */
/*                  p(alpha) = N(a, R)                                                                                           */
/*                                                                                                                               */
/*                                                                                                                               */
/* Let I = N*Sigma^{-1} + R^{-1}                                                                                                 */
/* Then p(alpha | ...) = Normal with                                                                                             */
/*                       var(alpha | ...) = I^{-1}                                                                               */
/*                         E(alpha | ...) = I^{-1} * (N*Sigma^{-1}*b(bar) + R^{-1}*a)                                            */
/*                                          where b(bar) = (1/N)*sum(b[i])                                                       */
/*                       That is, E(alpha | ...) solves I*x = (N*Sigma^{-1}*b(bar) + R^{-1}*a)                                   */
/*                                                                                                                               */
/* ----------------------------------------------------------------------------------------------------------------------------- */
/*                                                                                                                               */
/*         alpha:  means of normal effects to be updated                                                                         */
/*             I:  working matrix dim x dim                                                                                      */
/*             U:  working vector of length dim                                                                                  */
/*      PropMean:  working vector of length dim                                                                                  */
/*        ObsBar:  array of length dim holding b(bar)                                                                            */
/*        InvVar:  LT(dim) holding Sigma^{-1}                                                                                    */
/*    PriorWMean:  array of length dim with R^{-1}*a                                                                             */
/*   PriorInvVar:  array with R^{-1}                                                                                             */
/*           if diagPriorVar = true, this is array of length dim                                                                 */
/*           if diagPriorVar = false, this is array of length LT(dim) holding the lower triangle of R^{-1}                       */
/*  diagPriorVar:  if true then R^{-1} is diagonal, if false then R^{-1} is general symmetric PD matrix                          */
/*                                                                                                                               */
/*          nObs:  = N, number of observations                                                                                   */
/*           dim:  dimension ofthe normal distribution                                                                           */
/*                                                                                                                               */
/* ----------------------------------------------------------------------------------------------------------------------------- */
void
update_norm_mean1(double *alpha, 
                  MatrixLT<double> *I,
                  MatrixRect<double> *U,
                  MatrixRect<double> *PropMean,
                  const double *ObsBar,
                  const double *InvVar,
                  const double *PriorWMean,
                  const double *PriorInvVar,
                  const bool &diagPriorVar,
                  const int &nObs,
                  const int &dim,
                  const char *caller)
{
  static int i, Rank;
  static double *IP;
  static const double *dP;
  
  /*** N*Sigma^{-1}:  ***/
  IP = I->a();
  dP = InvVar;
  for (i = 0; i < I->length(); i++){
    *IP = nObs * (*dP);
    IP++;
    dP++;
  }

  /*** N*Sigma^{-1}*b(bar):  ***/
  Ab2(PropMean->a(), I, ObsBar);

  /*** Add R^{-1} to N*Sigma^{-1}: ***/
  if (diagPriorVar) I->add_b2diag(PriorInvVar);
  else              I->add_A2(PriorInvVar);
    
  /*** Add R^{-1}*a to N*Sigma^{-1}*b(bar):  ***/
  PropMean->add_A2(PriorWMean);

  /*** Sample a new value ***/  
  U->rNorm(0, 1);
  Rank = I->chol_solve(PropMean->a(), alpha, 1, 0, U->a());
  if (Rank < I->nrow()){
    throw returnR(caller, "Trap in update_norm_mean1. Variance matrix of the full conditional is not PD.", 1);
  }

  return;
}


/* ----------------------------------------------------------------------------------------------------------------------------- */
/*                                                                                                                               */
/* update_norm_invvar_wishart: Update of inverse variances of normal effects where the inverse variance itself has               */
/*                      a conjugate Wishart prior with 'nu' degrees of freedom and scale matrix 'S'                              */
/*                      (parametrization of Gelman, 2005)                                                                        */
/*                                                                                                                               */
/*    p(b[i] | alpha, Sigma) = N(alpha, Sigma), i = 1,...,N                                                                      */
/*             p(Sigma^{-1}) = Wishart(nu, S),                  that is E(Sigma^{-1}) = nu*S                                     */
/*                                                                                                                               */
/*       p(Sigma^{-1}) \propto det(Sigma^{-1})^((nu-k-1)/2) * exp(-0.5*tr(S^{-1}*Sigma^{-1})),                                   */
/*                             where k = nrow(Sigma) = ncol(Sigma)                                                               */
/*                                                                                                                               */
/*       Sigma^{-1} | ... is Wishart with (N + nu) degrees of freedom                                                            */
/*                        and scale matrix equal to (sum((b[i]-alpha)*(b[i]-alpha)') + S^{-1})^{-1}                              */
/*                                                                                                                               */
/* ----------------------------------------------------------------------------------------------------------------------------- */
/*                                                                                                                               */
/*        InvVar:  array of length LT(dim) with the Sigma^{-1} matrix to be updated                                              */
/*     InvS_Full:  working matrix LT(dim)                                                                                        */
/*     Obs_alpha:  working array of length dim*N                                                                                 */
/* work_rwishart:  working array of length 2*LT(dim) + dim*dim                                                                   */
/*           Obs:  matrix dim x N stored in an array with b[i], i=1,...,N in columns                                             */
/*         alpha:  array of length dim with the mean of b[i]                                                                     */
/*    df_Wishart:  degrees of freedom of the prior Wishart distribution                                                          */
/*  InvS_Wishart:  array of length LT(dim) with the S^{-1} matrix (inverse scale matrix of the Wishart distribution)             */
/*                                                                                                                               */
/*          nObs:  = N, number of observations                                                                                   */
/*           dim:  dimension of the normal distribution                                                                          */
/*                                                                                                                               */
/* ----------------------------------------------------------------------------------------------------------------------------- */
void
update_norm_invvar_wishart(double *InvVar,
                           MatrixLT<double> *InvS_Full,
                           double *Obs_alpha,
                           double *work_rwishart,
                           const double *Obs,
                           const double *alpha,
                           const double *df_Wishart,
                           const double *InvS_Wishart,
                           const int &nObs,
                           const int &dim,
                           const char *caller)
{
  static const double df_Full = *df_Wishart + nObs;

  /*** Inverse scale matrix of the full conditional Wishart distribution ***/
  NSampleVar(InvS_Full, Obs_alpha, Obs, alpha, nObs);
  InvS_Full->add_A2(InvS_Wishart);

  /*** Sample from the Wishart distribution  ***/
  Mvtdist3::rwishart3(InvVar, work_rwishart, &df_Full, InvS_Full->a(), &dim, 1);
  //rwishart2(InvVar, work_rwishart, &df_Full, InvS_Full);      /** original: removed on 12/01/2007, work_wishart had to be of length 2*LT(dim) + dim*dim **/

  return;
}


/* ----------------------------------------------------------------------------------------------------------------------------- */
/*                                                                                                                               */
/* update_norm_invvar_sduniform: Update of inverse covariance matrix  of normal effects where the inverse covariance matrix      */
/*                               is diagonal and square roots of diagonal ellemets are assumed to have uniform prior             */
/*                               Unif(0, S)                                                                                      */
/*                                                                                                                               */
/*    p(b[i] | alpha, Sigma) = N(alpha, Sigma), i = 1,...,N                                                                      */
/*                     Sigma = diag(sigma[1]^2,...,sigma[q]^2)                                                                   */
/*               p(sigma[i]) = Unif(0, S[i]),                                                                                    */
/*          p(sigma[i]^{-2}) = Gamma(-0.5, 0)*I[sigma[i]^{-2} > S[i]^{-2}]                                                       */
/*                                                                                                                               */
/*    p(sigma[i]^{-2} | ...) = Gamma(-0.5+0.5*N, 0.5*sum(b[i] - alpha[i])^2)*I[sigma[i]^{-2} > S[i]^{-2}]                        */
/*                                                                                                                               */
/* ----------------------------------------------------------------------------------------------------------------------------- */
/*                                                                                                                               */
/*        InvVar[LT(dim)]:  array with the Sigma^{-1} matrix to be updated                                                       */
/*                 it is assumed that this matrix is diagonal, off-diagonal elements are not modified                            */
/*              rate[dim]:  working array of length dim                                                                          */
/*                                                                                                                               */
/*             Obs[dim*N]:  matrix dim x N stored in an array with b[i], i=1,...,N in columns                                    */
/*             alpha[dim]:  array with the mean of b[i]                                                                          */
/*             invS2[dim]:  array of length dim with S[i]^{-2}                                                                   */
/*                                                                                                                               */
/*          nObs:  = N, number of observations                                                                                   */
/*           dim:  dimension of the normal distribution                                                                          */
/*                                                                                                                               */
/* ----------------------------------------------------------------------------------------------------------------------------- */
void
update_norm_invvar_sduniform(double *InvVar,
                             double *rate,
                             const double *Obs,
                             const double *alpha,
                             const double *invS2,
                             const int &nObs,
                             const int &dim,
                             const char *caller)
{
  static const double shape = 0.5*(nObs - 1);

  int j;
  double tempd;

  /*** Reset rate ***/
  double *rateP = rate;
  *rateP = 0.0;
  for (j = 1; j < dim; j++){
    rateP++;
    *rateP = 0.0;
  }

  /*** Compute rate ***/
  const double *alphaP;
  const double *obsP = Obs - 1;
  for (int i = 0; i < nObs; i++){
    obsP++;
    alphaP = alpha;
    rateP = rate;
  
    tempd = *obsP - (*alphaP);
    *rateP += tempd * tempd;
    for (j = 1; j < dim; j++){
      obsP++;
      alphaP++;
      rateP++;
      tempd = *obsP - (*alphaP);
      *rateP += tempd * tempd;
    }
  }

  /*** Sample new inverse variances ***/
  double *invsigma2P = InvVar;
  const double *invS2P = invS2;
  rateP = rate;
  for (j = dim; j > 0; j--){
    *rateP *= 0.5;
    Random::rltruncGamma(invsigma2P, &shape, rateP, invS2P, &_AK_ONE_INT);
    invsigma2P += j;
    rateP++;
    invS2P++;
  }

  return;
}

}  /*** end of the namespace mcmc_common ***/
