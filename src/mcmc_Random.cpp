/*** mcmc_Random.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  30/03/2007
//                            
// PURPOSE: Functions needed for the MCMC in the models with random effects
//          * used mainly by objects of classes RandomCL and RandomPoiss
//
//                   updateMeanRE1:  26/09/2006, 01/04/2007
//                 updateInvVarRE1:  26/09/2006, 01/04/2007
//
//                   updateMeanRE2:  29/10/2006, 01/04/2007
//                 updateMeanRE2Bi:  30/03/2007, 01/04/2007
//                 updateInvVarRE2:  30/10/2006, 01/04/2007
//               updateInvVarRE2Bi:  01/04/2007
//
//        full_Gspline_InvVar_pars:  30/10/2006, 30/03/2007
//    full_Gspline_InvVar_logdens0:  30/10/2006, 30/03/2007
//    full_Gspline_InvVar_logdens3:  30/10/2006, 30/03/2007
//
/* ********************************************************************************* */

#include "mcmc_Random.h"

namespace mcmc_Random{


/***** =============================================== *****/
/*****                                                 *****/
/***** MODELS WITH NORMALLY DISTRIBUTED RANDOM EFFECTS *****/
/*****                                                 *****/
/***** =============================================== *****/

/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* updateMeanRE1: Update of means of normally distributed random effects                                                    */
/*                                                                                                                          */
/*       Means are assumed to be fixed or have a normal prior with diagonal covariance matrix                               */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
//
//  REMean[nTheta]:             Vector of means of random effects to be updated
//
//  ThetaBar[nTheta]:           Working space to store the sample average of Theta components
//  PropMean[nTheta]:           Working space to store the mean of the proposal
//  U[nTheta]:                  Working space to store the score
//  I[LT(nTheta)]:              Working space to store the Hessian
//
//  Theta[nTheta, N]:           Matrix of values of random effects
//
//  REInvVar[LT(nTheta)]:       Covariance matrix of the random effects
//  REMeanPriorWMean[nTheta]:   Vector of prior inverse variances of means of random effects weighted by means
//  REMeanPriorInvVar[nTheta]:  Vector of prior inverse variances of means of random effects
//  prior_for_REMean:           Value from mcmc_Random::priorForMeanTypes
//
//  nTheta:                     Dimension of the random effects
//  N:                          Number of clusters
//
void
updateMeanRE1(MatrixRect<double> *REMean,                  MatrixRect<double> *ThetaBar,      
              MatrixRect<double> *PropMean,                MatrixRect<double> *U,                         MatrixLT<double> *I,
              const MatrixRect<double> *Theta,             const MatrixLT<double> *REInvVar,
              const MatrixRect<double> *REMeanPriorWMean,  const MatrixRect<double> *REMeanPriorInvVar,   const int &prior_for_REMean,
              const int &nTheta,                           const int &N)
{
  /*** Update _ThetaBar ***/
  Theta->mean(ThetaBar->a(), 1);

  
  switch (prior_for_REMean){
  case mcmc_Random::_Fixed_:
    return;  

  case mcmc_Random::_Normal_:
    /*** Update REMean = a ***/
    mcmc_common::update_norm_mean1(REMean->a(), I, U, PropMean, ThetaBar->aconst(), REInvVar->aconst(), 
                                   REMeanPriorWMean->aconst(), REMeanPriorInvVar->aconst(), true, 
                                   N, nTheta, "mcmc_Random::updateMeanRE1");
    return;

  default:
    REprintf("prior_for_REMean = %d\n", prior_for_REMean);
    throw returnR("Error in mcmc_Random.cpp: updateMeanRE1(). Update is not implemented for required prior_for_REMean", 1);
  }
}


/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* updateInvVarRE1: Update of inverse variance of normally distributed random effects                                       */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
//
// REInvVar[LT(nTheta)]:       Inverse covariance matrix of the random effects to be updated
// REInvVarL[LT(nTheta)]:      Cholesky decomposition of the inverse covariance matrix (updated as well)
//
// Theta_REMean[nTheta, N]:    Working space to store random effects minus their mean
// U[nTheta]:                  Working space
// I[LT(nTheta)]:              Working space
// work_rwishart[??]:          Working space, see RandomCL.h or RandomPoiss.h for its minimal length
//
// Theta[nTheta, N]:           Matrix of values of random effects
// REMean[nTheta]:             Vector of means of random effects
//
// REInvVarPriorDF[??]:        Parameters of the prior for the inverse variance of random effects, see RandomCL.h or RandomPoiss.h for explanation
// REInvVarPriorInvScale[??]:  Parameters of the prior for the inverse variance of random effects, see RandomCL.h or RandomPoiss.h for explanation
// prior_for_REMean:           Value from mcmc_Random::priorForInvVarTypes
// 
//  nTheta:                    Dimension of the random effects
//  N:                         Number of clusters
//
void
updateInvVarRE1(MatrixLT<double> *REInvVar,                 MatrixLT<double> *REInvVarL,
                MatrixRect<double> *Theta_REMean,           MatrixRect<double> *U,        MatrixLT<double> *I,  MatrixRect<double> *work_rwishart,
                const MatrixRect<double> *Theta,            const MatrixRect<double> *REMean,
                const MatrixRect<double> *REInvVarPriorDF,  const MatrixLT<double> *REInvVarPriorInvScale,      const int &prior_for_REInvVar,
                const int &nTheta,                          const int &N)
{
  static int j, Rank;
  static double *invsigma2P;

  switch (prior_for_REInvVar){
  case mcmc_Random::_Fixed:
    return;

  case mcmc_Random::_Wishart:

    /*** Update _REInvVar = D^{-1}                                 ***/
    /*** as by product, also _Theta_REMean = b - E(b) is updated   ***/
    mcmc_common::update_norm_invvar_wishart(REInvVar->a(), I, Theta_REMean->a(), work_rwishart->a(), 
                                            Theta->aconst(), REMean->aconst(), REInvVarPriorDF->aconst(), REInvVarPriorInvScale->aconst(),
                                            N, nTheta, "mcmc_Random::updateInvVarRE1");

    /*** Compute Cholesky decomposition of D^{-1} ***/
    REInvVarL->array2mat(REInvVar->aconst());
    Rank = REInvVarL->cholesky(0);
    if (Rank < nTheta){
        Rprintf("WARNING: mcmc_Random::updateInvVarRE1(...), sampled _REInvVar is not of full rank\n");
        Rprintf("_REInvVar:\n");
        REInvVar->print(0);
    }
    return;

  case mcmc_Random::_SDUnif:

    /*** Update REInvVar = D^{-1}, diagonal elements  ***/
    mcmc_common::update_norm_invvar_sduniform(REInvVar->a(), U->a(), 
                                              Theta->aconst(), REMean->aconst(), REInvVarPriorInvScale->aconst(),
                                              N, nTheta, "mcmc_Random::updateInvVarRE1");

    /*** Compute Cholesky decomposition of D^{-1}, that is compute square roots of diagonal elements ***/
    REInvVarL->array2mat(REInvVar->aconst());
    invsigma2P = REInvVarL->a();
    for (j = nTheta; j > 0; j--){
      if (*invsigma2P <= 0){
        Rprintf("WARNING: mcmc_Random::updateInvVarRE1(...), sampled _REInvVar is not of full rank\n");
        Rprintf("_REInvVar:\n");
        REInvVar->print(0);
      }
      *invsigma2P = sqrt(*invsigma2P);
      invsigma2P += j;
    }
    return;

  case mcmc_Random::_GammaIndep:
    throw returnR("Error in mcmc_Random.cpp: mcmc_Random::updateInvVarRE1(). Independent gamma priors for normal random effects inverse variances not implemented.", 1);   

  default:
    REprintf("prior_for_REInvVar=%d\n", prior_for_REInvVar);
    throw returnR("Error in mcmc_Random.cpp: mcmc_Random::updateInvVarRE1(). Update is not implemented for required prior_for_REInvVar", 1);
  }
}


/***** =============================================== *****/
/*****                                                 *****/
/***** MODELS WITH G-SPLINE DISTRIBUTED RANDOM EFFECTS *****/
/*****                                                 *****/
/***** =============================================== *****/

/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* updateMeanRE2: Update of intercepts of G-spline distributed random effects                                               */
/*                                                                                                                          */
/*      * Intercepts are assumed to have a normal prior with diagonal covariance matrix                                     */
/*      * Margins of the G-spline are assumed to be independent, given the other G-spline parameters                        */
/*        (e.g., copula)                                                                                                    */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
//
//  REMean[nTheta]:             Vector of means of random effects to be updated
//
//  ThetaBar[nTheta]:           Working space to store the sample average of Theta components
//
//  Theta[nTheta, N]:           Matrix of values of random effects
//
//  REMeanPriorWMean[nTheta]:   Vector of prior inverse variances of means of random effects weighted by means
//  REMeanPriorInvVar[nTheta]:  Vector of prior inverse variances of means of random effects
//  prior_for_REMean:           Value from mcmc_Random::priorForMeanTypes
//
//  Gspl:                       G-spline distribution of random effects
//  alloc[nTheta,N]:            Allocations for each margin on the scale -K,...,K                          
//
//  nTheta:                     Dimension of the random effects
//  N:                          Number of clusters
//
void
updateMeanRE2(MatrixRect<double> *REMean,                  MatrixRect<double> *ThetaBar,
              const MatrixRect<double> *Theta,
              const MatrixRect<double> *REMeanPriorWMean,  const MatrixRect<double> *REMeanPriorInvVar,   const int &prior_for_REMean,
              const Gspline2 *Gspl,                        const MatrixRect<int> *alloc,
              const int &nTheta,                           const int &N)
{
  static int i, j;
  static double *ThetaBarP, *REMeanP;
  static const double *ThetaP, *inv_sigma2_d2P, *REMeanPriorWMeanP, *REMeanPriorInvVarP;
  static const int *allocP, *KP;
  static const MatrixRect<double> *d_knotsM;

  static MatrixLT<double>   Ione        = MatrixLT<double>(1);
  static MatrixRect<double> Uone        = MatrixRect<double>(1, 1);
  static MatrixRect<double> PropMeanone = MatrixRect<double>(1, 1);

  switch (prior_for_REMean){
  case mcmc_Random::_Fixed_:
    return;  

  case mcmc_Random::_Normal_:

    /*** Compute sum(b[j,i] - d[j]*knot[alloc[j,i]]) = average of "observations" ***/
    ThetaBarP = ThetaBar->a();
    for (j = 0; j < nTheta; j++){
      *ThetaBarP = 0.0;
      ThetaBarP++;
    }

    ThetaP = Theta->aconst();
    allocP = alloc->aconst();
    for (i = 0; i < N; i++){
      KP = Gspl->KAconst();
      d_knotsM = Gspl->d_knotsconst();
      ThetaBarP = ThetaBar->a();
      for (j = 0; j < nTheta; j++){
        *ThetaBarP += (*ThetaP) - d_knotsM->aconst()[*allocP + (*KP)];
 
        ThetaP++;
        allocP++;

        KP++;
        d_knotsM++;
        ThetaBarP++;
      }
    }

    /*** Compute (1/N)*ThetaBar and update REMean, which must be done one-by-one! ***/
    ThetaBarP = ThetaBar->a();
    inv_sigma2_d2P     = Gspl->inv_sigma2_d2Aconst();
    REMeanPriorWMeanP  = REMeanPriorWMean->aconst();
    REMeanPriorInvVarP = REMeanPriorInvVar->aconst();
    REMeanP            = REMean->a();
    for (j = 0; j < nTheta; j++){  
      *ThetaBarP /= N;
      mcmc_common::update_norm_mean1(REMeanP, &Ione, &Uone, &PropMeanone, ThetaBarP, inv_sigma2_d2P, 
                                     REMeanPriorWMeanP, REMeanPriorInvVarP, true, 
                                     N, 1, "mcmc_Random::updateMeanRE2");

      ThetaBarP++;
      inv_sigma2_d2P++;
      REMeanPriorWMeanP++;
      REMeanPriorInvVarP++;
      REMeanP++;
    }

    return;  

  default:
    REprintf("prior_for_REMean = %d\n", prior_for_REMean);
    throw returnR("Error in mcmc_Random.cpp: updateMeanRE2(). Update is not implemented for required prior_for_REMean", 1);
  }
}


/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* updateMeanRE2Bi: Update of intercepts of BIVARIATE G-spline distributed random effects                                   */
/*                                                                                                                          */
/*      * Intercepts are assumed to have a normal prior with diagonal covariance matrix                                     */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
//
//  REMean[nTheta]:             Vector of means of random effects to be updated
//
//  ThetaBar[nTheta]:           Working space to store the sample average of Theta components
//  PropMean[nTheta]:           Working space to store the mean of the proposal
//  U[nTheta]:                  Working space to store the score
//  I[LT(nTheta)]:              Working space to store the Hessian
//
//  Theta[nTheta, N]:           Matrix of values of random effects
//
//  REMeanPriorWMean[nTheta]:   Vector of prior inverse variances of means of random effects weighted by means
//  REMeanPriorInvVar[nTheta]:  Vector of prior inverse variances of means of random effects
//  prior_for_REMean:           Value from mcmc_Random::priorForMeanTypes
//
//  Gspl:                       G-spline distribution of random effects
//  alloc[N]:                   Vector of single index (0,...,Gspl->total_length-1) allocations
//
//  nTheta:                     Dimension of the random effects
//  N:                          Number of clusters
//
void
updateMeanRE2Bi(MatrixRect<double> *REMean,                  MatrixRect<double> *ThetaBar,
                MatrixRect<double> *PropMean,                MatrixRect<double> *U,                         MatrixLT<double> *I,
                const MatrixRect<double> *Theta,
                const MatrixRect<double> *REMeanPriorWMean,  const MatrixRect<double> *REMeanPriorInvVar,   const int &prior_for_REMean,
                const BiGspline2 *Gspl,                      const MatrixRect<int> *alloc,
                const int &nTheta,                           const int &N)
{
  static int i, j, length0, i0, i1;
  static double *ThetaBarP;
  static double InvVar[3];
  static const double *ThetaP, *d_knotsP0, *d_knotsP1, *inv_sigma2_d2P;
  static const int *allocP;

  if (nTheta != BiGspline2A::_dim){
    REprintf("nTheta=%d,  BiGspline2A::_dim=%d\n", nTheta, BiGspline2A::_dim);
    throw returnR("Error in mcmc_Random.cpp: updateMeanRE2Bi(). Not implemented for this dimension.", 1);
  }

  switch (prior_for_REMean){
  case mcmc_Random::_Fixed_:
    return;  

  case mcmc_Random::_Normal_:

    /*** Compute sum(b[j,i] - d[j]*knot[alloc[j,i]]) = average of "observations" ***/
    ThetaBarP    = ThetaBar->a();
    ThetaBarP[0] = 0.0;
    ThetaBarP[1] = 0.0;

    ThetaP    = Theta->aconst();
    allocP    = alloc->aconst();
    d_knotsP0 = Gspl->d_knotsconst()[0].aconst();
    d_knotsP1 = Gspl->d_knotsconst()[1].aconst();
    length0   = Gspl->lengthconst()[0];
    for (i = 0; i < N; i++){
      i0 = *allocP % length0;
      i1 = *allocP / length0;
      ThetaBarP[0] += *ThetaP - d_knotsP0[i0];
      ThetaP++;
      ThetaBarP[1] += *ThetaP - d_knotsP1[i1];
      ThetaP++;

      allocP++;
    }

    /*** Compute (1/N)*ThetaBar ***/
    ThetaBarP[0] /= N;
    ThetaBarP[1] /= N;

    /*** Create inverse covariance matrix of "observations" ***/
    inv_sigma2_d2P = Gspl->inv_sigma2_d2const();
    InvVar[0] = inv_sigma2_d2P[0];
    InvVar[1] = 0;
    InvVar[2] = inv_sigma2_d2P[1];  
  
    /*** Update (jointly) REMean ***/
    mcmc_common::update_norm_mean1(REMean->a(), I, U, PropMean, ThetaBarP, InvVar,
                                   REMeanPriorWMean->aconst(), REMeanPriorInvVar->aconst(), true, 
                                   N, nTheta, "mcmc_Random::updateMeanRE2Bi");
 
    return;

  default:
    REprintf("prior_for_REMean = %d\n", prior_for_REMean);
    throw returnR("Error in mcmc_Random.cpp: updateMeanRE2Bi(). Update is not implemented for required prior_for_REMean", 1);
  }
}


/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* updateInvVarRE2: Update of inverse variance of G-spline distributed random effects                                       */
/*                                                                                                                          */
/*      * Margins of the G-spline are assumed to be independent, given the other G-spline parameters                        */
/*        (e.g., copula)                                                                                                    */
/*                                                                                                                          */
/*      * Only diagonal elements of REInvVar are updated, the rest is left unchanged and Cholesky decomposition             */
/*        of REInvVar assumes that the off-diagonal elements are equal to zero                                              */
/*                                                                                                                          */
/*      * Update is performed by the slice sampling from the full conditional distribution                                  */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
//
// REInvVar[LT(nTheta)]:       Inverse covariance matrix of the random effects to be updated
// REInvVarL[LT(nTheta)]:      Cholesky decomposition of the inverse covariance matrix (updated as well)
// REVar[LT(nTheta)]:          Covariance matrix of the random effects (updated as well)
//
// Gspl:                       G-spline distribution of random effects
//                             * its elements _d_knots, _sigma_d, _inv_sigma2_d2 are modified
//
// work_invVarSlice[??]:       Working space, see RandomCL.h or RandomPoiss.h for its minimal length
//
// Theta[nTheta, N]:           Matrix of values of random effects
// REMean[nTheta]:             Vector of means of random effects
//
// REInvVarPriorDF[??]:        Parameters of the prior for the inverse variance of random effects, see RandomCL.h or RandomPoiss.h for explanation
// REInvVarPriorInvScale[??]:  Parameters of the prior for the inverse variance of random effects, see RandomCL.h or RandomPoiss.h for explanation
// prior_for_REMean:           Value from mcmc_Random::priorForInvVarTypes
// 
// alloc[nTheta,N]:            Allocations for each margin on the scale -K,...,K
// 
// nTheta:                     Dimension of the random effects
// N:                          Number of clusters
//
// overrelax:                  Should overrelaxation be performed?
//                             ORIGINALLY:  overrelax = 1*((*iter/_k_overrelax_scale[j]) != 0);
//
void
updateInvVarRE2(MatrixLT<double> *REInvVar,                 MatrixLT<double> *REInvVarL,                    MatrixLT<double> *REVar,
                Gspline2 *Gspl,                             MatrixRect<double> *work_invVarSlice,           
                const MatrixRect<double> *Theta,            const MatrixRect<double> *REMean,
                const MatrixRect<double> *REInvVarPriorDF,  const MatrixLT<double> *REInvVarPriorInvScale,  const int &prior_for_REInvVar,
                const MatrixRect<int> *alloc,
                const int &nTheta,                          const int &N,                                   const int &overrelax)
{
  static int j, j_, itmp, err_nr, iter_nr;
  static int ipars[1];
  static double slice[3];
  static double gx[2];
  static double gx0, horiz, tmp, dgx, newInvVar;
  static double *invVarP, *invVarLP, *VarP, *pars;
  static const double *zeta_1, *sqrt_eta, *xi_half;

  static const MatrixRect<double> *knotsM;
  static const double *sigmaP, *invsigma2P, *knotsP;
  static const int *KP;
  static MatrixRect<double> *d_knotsM;
  static double *sigma_dP, *inv_sigma2_d2P, *d_knotsP;

  pars = work_invVarSlice->a();

  ipars[0] = 0;
  switch (prior_for_REInvVar){
  case mcmc_Random::_Fixed:
    return;

  case mcmc_Random::_Wishart:
    throw returnR("Error in mcmc_Random.cpp: updateInvVarRE2. Not implemented for _Wishart prior.", 1);

  case mcmc_Random::_SDUnif:
    ipars[0] = 1;
  case mcmc_Random::_GammaIndep:

    /** Compute parameters of the full conditional distributions for all inverse variances  **/
    mcmc_Random::full_Gspline_InvVar_pars(pars, 0, alloc, 
                                          Gspl->knotsconst(), Gspl->KAconst(), Gspl->invsigma2Aconst(),
                                          Theta, nTheta, N, 
                                          REMean->aconst(),
                                          prior_for_REInvVar, REInvVarPriorDF->aconst(), REInvVarPriorInvScale->aconst());

    /** Loop over dimensions **/
    zeta_1    = pars;
    sqrt_eta  = zeta_1 + 1;
    xi_half   = sqrt_eta + 1;
    invVarP   = REInvVar->a();
    invVarLP  = REInvVarL->a();
    VarP      = REVar->a();

    knotsM         = Gspl->knotsconst();
    sigmaP         = Gspl->sigmaAconst();
    invsigma2P     = Gspl->invsigma2Aconst();
    KP             = Gspl->KAconst();
    d_knotsM       = Gspl->d_knots();
    sigma_dP       = Gspl->sigma_dA();
    inv_sigma2_d2P = Gspl->inv_sigma2_d2A();
   
    for (j = nTheta; j > 0; j--){

      /** Evaluate the full conditional distribution in the current point and sample the level defining the slice **/
      mcmc_Random::full_Gspline_InvVar_logdens0(invVarP, &gx0, zeta_1, ipars);
      horiz = gx0 - rexp(1);       

      /** Compute mode of the full conditional distribution, store it in slice[2] **/    
      if ((*zeta_1) <= 0) throw returnR("Trap in mcmc_Random.cpp: updateInvVarRE2(). Zeta parameter for scale update <= 1. Is your sample size > 2?", 1);
      tmp = -(*xi_half) + sqrt((*xi_half)*(*xi_half) + 4*(*sqrt_eta)*(*sqrt_eta)*(*zeta_1));
      if (tmp < _AK_expEMIN) throw returnR("Trap in mcmc_Random.cpp: updateInvVarRE2(). Cannot allocate mode of the full conditional.", 1);
      slice[2] = (4*(*zeta_1)*(*zeta_1))/(tmp*tmp);

      /** Give initial guesses for the interval defining the slice,                   **/
      /** starting left limit x0 of the interval must be such that g(x0) < horiz,     **/
      /** otherwise the Newton_raphson solver fails due to going below zero           **/
      /** Set initial right limit either to *invVarP or to mode + 2*std. dev.         **/
      /** but always to the right from the mode                                       **/
      if (*invVarP < slice[2]){
        dgx = (*zeta_1)/(slice[2]*slice[2]) + (*xi_half)/(2*slice[2]*sqrt(slice[2]));        /* l"(x_mode) */
        slice[1] = slice[2] + 2/sqrt(dgx);
      }
      else{
        slice[1] = *invVarP;
      }
      slice[0] = 0.5*slice[2];
      mcmc_Random::full_Gspline_InvVar_logdens0(slice, gx, zeta_1, ipars);
      while (gx[0] >= horiz && slice[0] > mcmc_Random::_zero_invvariance){
        slice[0] *= 0.5;
        mcmc_Random::full_Gspline_InvVar_logdens0(slice, gx, zeta_1, ipars);
      }

      /* Compute the interval defining the slice */
      itmp = (slice[0] <= mcmc_Random::_zero_invvariance ? 1 : 0);
      for (j_ = 1; j_ >= itmp; j_--){
        mcmc_Random::full_Gspline_InvVar_logdens3(slice+j_, gx+j_, &dgx, &tmp, zeta_1, ipars, 3);
	AK_Optim::solver_newton_raphson(slice+j_, gx+j_, &dgx, &horiz, zeta_1, ipars, mcmc_Random::full_Gspline_InvVar_logdens3, 
                                        &iter_nr, &mcmc_Random::_maxiter_solver_nr, &mcmc_Random::_toler_solver_nr, &_AK_expEMIN, &err_nr);
        if (err_nr >= 3){
          REprintf("\nerr_nr = %d\n", err_nr);
          REprintf("inverse variance[%d] = %f\n", nTheta-j, *invVarP);
          REprintf("mode = %f, horizontal = %f\n", slice[2], horiz);
          REprintf("xi[0]-1 = %f,  sqrt(xi[1]) = %f, xi[2]/2 = %f\n", *zeta_1, *sqrt_eta, *xi_half);
          throw returnR("Trap in mcmc_Random.cpp: updateInvVarRE2(). Unable to find an interval defining the slice.", 1);
        }
      }
      if (ipars[0]){
        if (slice[0] <= zeta_1[3]) slice[0] = zeta_1[3];     /* left limit of the interval below the truncation limit       */
        if (slice[1] <= zeta_1[3]){                          /* also right limit of the interval below the truncation limit */
          *invVarP = zeta_1[3];
          goto labelPOST_SAMPLE;
        }
      }

      /* Sample the new point */
      if (overrelax){
	Slice_sampler::ss_exact_overrelax(&newInvVar, slice, invVarP, &horiz, mcmc_Random::full_Gspline_InvVar_logdens0, zeta_1, ipars);
      }
      else{
        Slice_sampler::ss_exact_sample(&newInvVar, slice, gx, invVarP, &horiz, mcmc_Random::full_Gspline_InvVar_logdens0, zeta_1, ipars);
      }

      /* Update necessary quantities */
      *invVarP = newInvVar;

      labelPOST_SAMPLE:

      /*** Update Cholesky decomposition and covariance matrix itself ***/
      *invVarLP = sqrt(*invVarP);              
      *VarP     = 1/(*invVarP);

      /*** Update G-spline related quantities ***/
      *sigma_dP       = (*sigmaP)/(*invVarLP);
      *inv_sigma2_d2P = (*invsigma2P)*(*invVarP);
      knotsP          = knotsM->aconst();
      d_knotsP        = d_knotsM->a();
      for (j_ = -(*KP); j_ <= (*KP); j_++){
        *d_knotsP = *knotsP/(*invVarLP);
        knotsP++;
        d_knotsP++;
      }

      /*** Increase pointers ***/
      zeta_1    = xi_half + 2;
      sqrt_eta  = zeta_1 + 1;
      xi_half   = sqrt_eta + 1;
      invVarP   += j;
      invVarLP  += j;
      VarP      += j;

      knotsM++;
      sigmaP++;
      invsigma2P++;
      KP++;
      d_knotsM++;
      sigma_dP++;
      inv_sigma2_d2P++;
    }

    return;
  }
}


/* ------------------------------------------------------------------------------------------------------------------------ */
/*                                                                                                                          */
/* updateInvVarRE2Bi: Update of inverse variance of BIVARIATE G-spline distributed random effects                           */
/*                                                                                                                          */
/*      * Only diagonal elements of REInvVar are updated, the rest is left unchanged and Cholesky decomposition             */
/*        of REInvVar assumes that the off-diagonal elements are equal to zero                                              */
/*                                                                                                                          */
/*      * Update is performed by the slice sampling from the full conditional distribution                                  */
/*                                                                                                                          */
/* ------------------------------------------------------------------------------------------------------------------------ */
//
// REInvVar[LT(nTheta)]:       Inverse covariance matrix of the random effects to be updated
// REInvVarL[LT(nTheta)]:      Cholesky decomposition of the inverse covariance matrix (updated as well)
// REVar[LT(nTheta)]:          Covariance matrix of the random effects (updated as well)
//
// Gspl:                       G-spline distribution of random effects
//                             * its elements _d_knots, _sigma_d, _inv_sigma2_d2 are modified
//
// work_invVarSlice[??]:       Working space, see RandomCL.h or RandomPoiss.h for its minimal length
//
// Theta[nTheta, N]:           Matrix of values of random effects
// REMean[nTheta]:             Vector of means of random effects
//
// REInvVarPriorDF[??]:        Parameters of the prior for the inverse variance of random effects, see RandomCL.h or RandomPoiss.h for explanation
// REInvVarPriorInvScale[??]:  Parameters of the prior for the inverse variance of random effects, see RandomCL.h or RandomPoiss.h for explanation
// prior_for_REMean:           Value from mcmc_Random::priorForInvVarTypes
// 
// alloc[N]:                   Vector of single index (0,...,Gspl->total_length-1) allocations
// 
// nTheta:                     Dimension of the random effects
// N:                          Number of clusters
//
// overrelax:                  Should overrelaxation be performed?
//                             ORIGINALLY:  overrelax = 1*((*iter/_k_overrelax_scale[j]) != 0);
//
void
updateInvVarRE2Bi(MatrixLT<double> *REInvVar,                 MatrixLT<double> *REInvVarL,                    MatrixLT<double> *REVar,
                  BiGspline2 *Gspl,                           MatrixRect<double> *work_invVarSlice,           
                  const MatrixRect<double> *Theta,            const MatrixRect<double> *REMean,
                  const MatrixRect<double> *REInvVarPriorDF,  const MatrixLT<double> *REInvVarPriorInvScale,  const int &prior_for_REInvVar,
                  const MatrixRect<int> *alloc,
                  const int &nTheta,                          const int &N,                                   const int &overrelax)
{
  static int j, j_, itmp, err_nr, iter_nr;
  static int ipars[1];
  static double slice[3];
  static double gx[2];
  static double gx0, horiz, tmp, dgx, newInvVar;
  static double *invVarP, *invVarLP, *VarP, *pars;
  static const double *zeta_1, *sqrt_eta, *xi_half;

  static const MatrixRect<double> *knotsM;
  static const double *sigmaP, *invsigma2P, *knotsP;
  static const int *KP;
  static MatrixRect<double> *d_knotsM;
  static double *sigma_dP, *inv_sigma2_d2P, *d_knotsP;

  if (nTheta != BiGspline2A::_dim){
    REprintf("nTheta=%d,  BiGspline2A::_dim=%d\n", nTheta, BiGspline2A::_dim);
    throw returnR("Error in mcmc_Random.cpp: updateInvVarRE2Bi(). Not implemented for this dimension.", 1);
  }

  pars = work_invVarSlice->a();

  ipars[0] = 0;
  switch (prior_for_REInvVar){
  case mcmc_Random::_Fixed:
    return;

  case mcmc_Random::_Wishart:
    throw returnR("Error in mcmc_Random.cpp: updateInvVarRE2Bi. Not implemented for _Wishart prior.", 1);

  case mcmc_Random::_SDUnif:
    ipars[0] = 1;
  case mcmc_Random::_GammaIndep:

    /** Compute parameters of the full conditional distributions for all inverse variances  **/
    mcmc_Random::full_Gspline_InvVar_pars(pars, 1, alloc, 
                                          Gspl->knotsconst(), Gspl->Kconst(), Gspl->invsigma2const(),
                                          Theta, nTheta, N, 
                                          REMean->aconst(),
                                          prior_for_REInvVar, REInvVarPriorDF->aconst(), REInvVarPriorInvScale->aconst());

    /** Loop over dimensions **/
    zeta_1    = pars;
    sqrt_eta  = zeta_1 + 1;
    xi_half   = sqrt_eta + 1;
    invVarP   = REInvVar->a();
    invVarLP  = REInvVarL->a();
    VarP      = REVar->a();

    knotsM         = Gspl->knotsconst();
    sigmaP         = Gspl->sigmaconst();
    invsigma2P     = Gspl->invsigma2const();
    KP             = Gspl->Kconst();
    d_knotsM       = Gspl->d_knots();
    sigma_dP       = Gspl->sigma_d();
    inv_sigma2_d2P = Gspl->inv_sigma2_d2();
   
    for (j = nTheta; j > 0; j--){

      /** Evaluate the full conditional distribution in the current point and sample the level defining the slice **/
      mcmc_Random::full_Gspline_InvVar_logdens0(invVarP, &gx0, zeta_1, ipars);
      horiz = gx0 - rexp(1);       

      /** Compute mode of the full conditional distribution, store it in slice[2] **/    
      if ((*zeta_1) <= 0) throw returnR("Trap in mcmc_Random.cpp: updateInvVarRE2Bi(). Zeta parameter for scale update <= 1. Is your sample size > 2?", 1);
      tmp = -(*xi_half) + sqrt((*xi_half)*(*xi_half) + 4*(*sqrt_eta)*(*sqrt_eta)*(*zeta_1));
      if (tmp < _AK_expEMIN) throw returnR("Trap in mcmc_Random.cpp: updateInvVarRE2Bi(). Cannot allocate mode of the full conditional.", 1);
      slice[2] = (4*(*zeta_1)*(*zeta_1))/(tmp*tmp);

      /** Give initial guesses for the interval defining the slice,                   **/
      /** starting left limit x0 of the interval must be such that g(x0) < horiz,     **/
      /** otherwise the Newton_raphson solver fails due to going below zero           **/
      /** Set initial right limit either to *invVarP or to mode + 2*std. dev.         **/
      /** but always to the right from the mode                                       **/
      if (*invVarP < slice[2]){
        dgx = (*zeta_1)/(slice[2]*slice[2]) + (*xi_half)/(2*slice[2]*sqrt(slice[2]));        /* l"(x_mode) */
        slice[1] = slice[2] + 2/sqrt(dgx);
      }
      else{
        slice[1] = *invVarP;
      }
      slice[0] = 0.5*slice[2];
      mcmc_Random::full_Gspline_InvVar_logdens0(slice, gx, zeta_1, ipars);
      while (gx[0] >= horiz && slice[0] > mcmc_Random::_zero_invvariance){
        slice[0] *= 0.5;
        mcmc_Random::full_Gspline_InvVar_logdens0(slice, gx, zeta_1, ipars);
      }

      /* Compute the interval defining the slice */
      itmp = (slice[0] <= mcmc_Random::_zero_invvariance ? 1 : 0);
      for (j_ = 1; j_ >= itmp; j_--){
        mcmc_Random::full_Gspline_InvVar_logdens3(slice+j_, gx+j_, &dgx, &tmp, zeta_1, ipars, 3);
	AK_Optim::solver_newton_raphson(slice+j_, gx+j_, &dgx, &horiz, zeta_1, ipars, mcmc_Random::full_Gspline_InvVar_logdens3, 
                                        &iter_nr, &mcmc_Random::_maxiter_solver_nr, &mcmc_Random::_toler_solver_nr, &_AK_expEMIN, &err_nr);
        if (err_nr >= 3){
          REprintf("\nerr_nr = %d\n", err_nr);
          REprintf("inverse variance[%d] = %f\n", nTheta-j, *invVarP);
          REprintf("mode = %f, horizontal = %f\n", slice[2], horiz);
          REprintf("xi[0]-1 = %f,  sqrt(xi[1]) = %f, xi[2]/2 = %f\n", *zeta_1, *sqrt_eta, *xi_half);
          throw returnR("Trap in mcmc_Random.cpp: updateInvVarRE2Bi(). Unable to find an interval defining the slice.", 1);
        }
      }
      if (ipars[0]){
        if (slice[0] <= zeta_1[3]) slice[0] = zeta_1[3];     /* left limit of the interval below the truncation limit       */
        if (slice[1] <= zeta_1[3]){                          /* also right limit of the interval below the truncation limit */
          *invVarP = zeta_1[3];
          goto labelPOST_SAMPLE;
        }
      }

      /* Sample the new point */
      if (overrelax){
	Slice_sampler::ss_exact_overrelax(&newInvVar, slice, invVarP, &horiz, mcmc_Random::full_Gspline_InvVar_logdens0, zeta_1, ipars);
      }
      else{
        Slice_sampler::ss_exact_sample(&newInvVar, slice, gx, invVarP, &horiz, mcmc_Random::full_Gspline_InvVar_logdens0, zeta_1, ipars);
      }

      /* Update necessary quantities */
      *invVarP = newInvVar;

      labelPOST_SAMPLE:

      /*** Update Cholesky decomposition and covariance matrix itself ***/
      *invVarLP = sqrt(*invVarP);              
      *VarP     = 1/(*invVarP);

      /*** Update G-spline related quantities ***/
      *sigma_dP       = (*sigmaP)/(*invVarLP);
      *inv_sigma2_d2P = (*invsigma2P)*(*invVarP);
      knotsP          = knotsM->aconst();
      d_knotsP        = d_knotsM->a();
      for (j_ = -(*KP); j_ <= (*KP); j_++){
        *d_knotsP = *knotsP/(*invVarLP);
        knotsP++;
        d_knotsP++;
      }

      /*** Increase pointers ***/
      zeta_1    = xi_half + 2;
      sqrt_eta  = zeta_1 + 1;
      xi_half   = sqrt_eta + 1;
      invVarP   += j;
      invVarLP  += j;
      VarP      += j;

      knotsM++;
      sigmaP++;
      invsigma2P++;
      KP++;
      d_knotsM++;
      sigma_dP++;
      inv_sigma2_d2P++;
    }

    return;
  }
}


/***** ----------------------------------------------------------------------------------------------------------------------- *****/
/***** full_Gspline_InvVar_pars: Compute parameters of the full conditional distribution of                                    *****/
/*****                           all (overall)-2} parameters inverse variances                                                 *****/
/*****                           when the random effects (_Theta) have G-spline distribution                                   *****/
/*****                                                                                                                         *****/
/***** Taken (partially) from Gspline_update_Scale.cpp (full_Scale_pars function) of bayesSurv package                         *****/
/*****                                                                                                                         *****/
/***** ----------------------------------------------------------------------------------------------------------------------- *****/
//
// pars[4*dim] .............. computed parameters of the full conditional distribution for each d^{-2}
//          * raw parameters of the full conditional distribution are:
//                          xi[0] = shape + N/2		     
//		  	    xi[1]  = rate + 0.5*_invsigma2*sum(y-alpha)^2
//			    xi[2]   = invsigma2*sum[mu*(y-alpha)]        
//          * we return:
//                        pars[0+j*4] ..... xi[0] - 1
//                        pars[1+j*4] ..... sqrt(xi[1])
//                        pars[2+j*4] ..... xi[2]/2
//                        pars[3+j*4] ..... possibly truncation point 1/S^2
//
// BiGspline           if <> 0 then it is assumed that the G-spline is bivariate
//                     with a bivariate grid of knots
// alloc               if (BiGspline) alloc[1 x N] with single index component labels taking values 0, ..., total_length-1
//                     else           alloc[nTheta x N] with component labels taking values -K[i],...,K[i]
//                     
// knots[nTheta]       pointer to matrices 1 x (2*K[i] + 1) with knots for each margin               
// K[nTheta]           for each margin, number of knots on each side of the reference knot
// invsigma2[nTheta]   inverse of basis variances of the G-spline
//
//
// Theta[nTheta x N]   values of random effects      
// nTheta              dimension of the random effects
// N                   number of observations (clusters)
//
// REMean[nTheta]      mean of the random effects (_REMean component of the objects RandomCL or RandomPoiss)
//
// prior_for_REInvVar  specification of the prior for 'tau' parameters
//                     (value from mcmc_Random::priorForInvVarTypes)
// REInvVarPriorDF        array corresponding to _REInvVarPriorDF components of objects RandomCL or RandomPoiss
// REInvVarPriorInvScale  array corresponding to _REInvVarPriorInvScale components of objects RandomCL or RandomPoiss
//                
void
full_Gspline_InvVar_pars(double *pars, 
                         const int &BiGspline,             const MatrixRect<int> *alloc,  
                         const MatrixRect<double> *knots,  const int *K,                   const double *invsigma2,
                         const MatrixRect<double> *Theta,  const int &nTheta,              const int &N,
                         const double *REMean,
                         const int &prior_for_REInvVar,    const double *REInvVarPriorDF,  const double *REInvVarPriorInvScale)
{
  //static const int npars = 4;
  static int j, obs, length0, ind0, ind1;
  static double y_alpha;

  static const double *y_obs, *intcptP0, *intcptP1, *invsigma2P, *knotsP0, *knotsP1, *GshapeP, *GinvscaleP;
  static double *xi1_0, *xi2_0, *xi1_1, *xi2_1;
  static const int *rp;
  static double *dP;
  
  /*** Reset xi[1] and xi[2] ***/
  dP = pars + 1;
  for (j = 0; j < nTheta; j++){
    *dP = 0.0;      //  pars[j*npars+1] = 0.0
    dP++;
    *dP = 0.0;      //  pars[j*npars+2] = 0.0
    dP += 3;
  }

  /*** Compute sum(y[i]-alpha)^2 = part of xi[1] and sum(knot[r[i]]*(y-alpha)) = part of xi[2] ***/
  y_obs = Theta->aconst();
  rp    = alloc->aconst();

  if (nTheta == 1){
    intcptP0 = REMean;    
    knotsP0  = knots->aconst();
    xi1_0    = pars + 1;
    xi2_0    = xi1_0 + 1;
    for (obs = 0; obs < N; obs++){
      y_alpha = (*y_obs) - (*intcptP0);      
      *xi1_0 += y_alpha * y_alpha;
      *xi2_0 += knotsP0[*rp + (*K)] * y_alpha;

      y_obs++;
      rp++;
    }
  }
  else{   /** nTheta > 1 **/
    if (BiGspline){
      length0 = 2*K[0] + 1;

      intcptP0 = REMean;
      intcptP1 = REMean + 1;
      knotsP0  = knots[0].aconst();
      knotsP1  = knots[1].aconst(); 
      xi1_0    = pars + 1;
      xi2_0    = xi1_0 + 1;
      xi1_1    = xi2_0 + 3;
      xi2_1    = xi1_1 + 1;     

      for (obs = 0; obs < N; obs++){
        ind0 = *rp % length0;
        ind1 = *rp / length0;

        y_alpha = (*y_obs) - (*intcptP0);
        *xi1_0 += y_alpha * y_alpha;
        *xi2_0 += knotsP0[ind0] * y_alpha;
        y_obs++;

        y_alpha = (*y_obs) - (*intcptP1);
        *xi1_1 += y_alpha * y_alpha;
        *xi2_1 += knotsP1[ind1] * y_alpha;
        y_obs++;
      
        rp++;
      }  
    }
    else{
      REprintf("nTheta=%d\n", nTheta);
      throw returnR("Error in mcmc_Random.cpp: full_Gspline_InvVar_pars. Only implemented for uni- and bivariate G-splines.", 1);      
    }
  }

  /*** Compute final values ***/
  switch (prior_for_REInvVar){
  case mcmc_Random::_Fixed:
    return;

  case mcmc_Random::_Wishart:
    throw returnR("Error in mcmc_Random.cpp: full_Gspline_InvVar_pars. Not implemented for _Wishart prior.", 1);
    
  case mcmc_Random::_SDUnif:
    dP = pars;
    GinvscaleP = REInvVarPriorInvScale;          /* 1/S[i]^2  (inverted squares of the truncation points) */
    invsigma2P = invsigma2;
    for (j = 0; j < nTheta; j++){
      *dP = -0.5 + N/2 - 1;                      /* xi[0] - 1 */
      dP++;

      *dP *= 0.5 * (*invsigma2P);
      *dP = sqrt(*dP);                           /* sqrt(xi[1]) */
      dP++;

      *dP *= 0.5 * (*invsigma2P);                /* xi[2]/2 */
      dP++;

      *dP = *GinvscaleP;                        /* 1/square of the truncation point */
      dP++;

      GinvscaleP++;
      invsigma2P++;
    }
    return;

  case mcmc_Random::_GammaIndep:
    dP = pars;
    GshapeP    = REInvVarPriorDF;
    GinvscaleP = REInvVarPriorInvScale;
    invsigma2P = invsigma2;
    for (j = 0; j < nTheta; j++){
      *dP = *GshapeP + N/2 - 1;                  /* xi[0] - 1 */
      dP++;

      *dP *= 0.5 * (*invsigma2P);
      *dP += *GinvscaleP;
      *dP = sqrt(*dP);                           /* sqrt(xi[1]) */
      dP++;

      *dP *= 0.5 * (*invsigma2P);                /* xi[2]/2 */
      dP+= 2;

      GshapeP++;
      GinvscaleP++;
      invsigma2P++;
    }
    return;

  default:
    throw returnR("Error in mcmc_Random.cpp: full_Gspline_InvVar_pars. Unimplemented prior.", 1);
  }
}


/***** ----------------------------------------------------------------------------------------- *****/
/***** full_Gspline_InvVar_logdens0: Compute log-density of a full conditional distribution      *****/
/*****                               of d^{-2} (overall variances)                               *****/
/***** ----------------------------------------------------------------------------------------- *****/
//
//  log(p(d^{-2}|...)) = (xi[0]-1)*log(d^{-2}) - (xi[1]*d^{-2} - xi[2]*sqrt(d^{-2})) + C1
//                     = (xi[0]-1)*log(d^{-2}) - (sqrt(xi[1]*d^{-2}) - 0.5*xi[2]/sqrt(xi[1]))^2 + C2
//
//  This function returns (xi[0]-1)*log(d^{-2}) - (sqrt(xi[1]*d^{-2}) - 0.5*xi[2]/sqrt(xi[1]))^2
//       
//   *** TRUNCATION: taken into account by      full_sigma_logdens0
//
// pars[4] ..... pars[0] = xi[0] - 1
//               pars[1] = sqrt(xi[1])
//               pars[2] = xi[2]/2
//               pars[3] = truncation point 1/S^2
// ipars[1] .... ipars[0] = 0/1 is truncated?
//
void
full_Gspline_InvVar_logdens0(const double* x,  double* yu,  const double* pars,  const int* ipars)
{
  static double tmp, sqrt_x;

  if (ipars[0] && (*x) <= pars[3]){
    *yu = R_NegInf;
    return;
  }

  if (*x <= _AK_expEMIN){
    *yu = R_NegInf;
    return;
  }

  sqrt_x = sqrt(*x);
  tmp = pars[1]*sqrt_x - pars[2]/pars[1];
  *yu = pars[0]*log(*x) - tmp*tmp;
  return;
}


// ***** full_Gspline_InvVar_logdens3: Compute log-density (and possibly also some derivatives) of a full conditional distribution  *****//
//                                     of d^{-2} (overall inverse variances)                                                        *****//
//       
//   *** Curently: full_Gspline_InvVar_logdens3 is nowhere used
//
//   *** TRUNCATION: taken into account
//
// pars[3] ..... pars[0] = xi[0] - 1
//               pars[1] = sqrt(xi[1])
//               pars[2] = xi[2]/2
// ipars[1] .... ipars[0] = 0/1 is truncated?
// what ........ 0 = compute l(x), l'(x), -l''(x)
//               1 = compute only l(x)
//               2 = compute only l'(x), -l''(x)
//               3 = compute only l(x), l'(x)
//
void
full_Gspline_InvVar_logdens3(const double* x,     double* yu,        double* ypu,     double* yppu,  
                             const double* pars,  const int* ipars,  const int& what)
{  
  static double tmp, sqrt_x, invx, inv_sqrt_x, zeta_1_invx;

  sqrt_x = sqrt(*x);
  switch (what){
  case 0:
    if ((ipars[0] && (*x) <= pars[3]) || (*x <= _AK_expEMIN)){
      *yu = R_NegInf;
      *ypu = R_NegInf;
      *yppu = R_NegInf;
    }
    else{
      tmp = pars[1]*sqrt_x - pars[2]/pars[1];
      *yu = pars[0]*log(*x) - tmp*tmp;

      invx = 1/(*x);
      zeta_1_invx = pars[0]*invx;
      inv_sqrt_x = 1/sqrt_x;
      *ypu = zeta_1_invx - pars[1]*pars[1] + pars[2]*inv_sqrt_x;

      *yppu = zeta_1_invx*invx + (pars[2]/2)*invx*inv_sqrt_x;
    }
    return;

  case 1:
    if ((ipars[0] && (*x) <= pars[3]) || (*x <= _AK_expEMIN)){
      *yu = R_NegInf;
    }
    else{
      tmp = pars[1]*sqrt_x - pars[2]/pars[1];
      *yu = pars[0]*log(*x) - tmp*tmp;
    }
    return;

  case 2:
    if ((ipars[0] && (*x) <= pars[3]) || (*x <= _AK_expEMIN)){
      *ypu = R_NegInf;
      *yppu = R_NegInf;
    }
    else{
      invx = 1/(*x);
      zeta_1_invx = pars[0]*invx;
      inv_sqrt_x = 1/sqrt_x;
      *ypu = zeta_1_invx - pars[1]*pars[1] + pars[2]*inv_sqrt_x;

      *yppu = zeta_1_invx*invx + (pars[2]/2)*invx*inv_sqrt_x;
    }
    return;

  case 3:
    if ((ipars[0] && (*x) <= pars[3]) || (*x <= _AK_expEMIN)){
      *yu = R_NegInf;
      *ypu = R_NegInf;
    }
    else{
      tmp = pars[1]*sqrt_x - pars[2]/pars[1];
      *yu = pars[0]*log(*x) - tmp*tmp;

      invx = 1/(*x);
      zeta_1_invx = pars[0]*invx;
      inv_sqrt_x = 1/sqrt_x;
      *ypu = zeta_1_invx - pars[1]*pars[1] + pars[2]*inv_sqrt_x;
    }
    return;

  default:
    throw returnR("Error in mcmc_Random.cpp: full_Gspline_InvVar_logdens3(). Incorrect what argument", 1);
  }
}

}  /*** end of the namespace mcmc_Random ***/


