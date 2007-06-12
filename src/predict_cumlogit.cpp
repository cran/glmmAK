/*** predict_cumlogit.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:                           22/01/2007
//    WORKING VERSION:                   31/01/2007
//    VERSION WITH BIVARIATE G-SPLINES:  12/04/2007
//
//
// PURPOSE: Cumulative logit model with random effects
//          * predictive probabilities
//
/* ********************************************************************************* */

#include "predict_cumlogit.h"

namespace predict_cumlogitA{

extern "C"{

// prob_value:               If (retValue || nqprob) prob_value[n*(C+1), niter] = values of all evaluated predictive probabilities
//                           Else                    prob_value[n*(C+1)]        = values of the last evaluated predictive probabilities
//                           
// prob_summary[(C+1)*n, 1+nqprob]:  Pointwise averages of the predictive probabilities
//                                   + quantiles (if nqprob>0)
// retValue:                         Do we want to return the values of the predictive probabilities at each iteration?
//
// qprob[nprob]:                Probabilities of required quantiles
// nqprob:                      Number of quantiles required
//
// X[p,n]:                          Covariates for fixed effects whose effect is proportional w.r.t. odds
// V[q,n]:                          Covariates for fixed effects whose effect is not proportional w.r.t. odds
// p[1]:                            Number of fixed effects whose effect is proportional w.r.t. odds
// q[1]:                            Number of fixed effects whose effect is not proportional w.r.t. odds
//
// G_dim[]:              Number of random effects and numbers of knots on each side of the zero knot for each margin of the G-spline
//                       (see argument dimPar of the parametric Gspline2 or BiGspline2 constructor)
//
// G_dpar[]:             Basis standard deviations and knots
//                       (see argument penMix_dPar of the parametric Gspline2 or BiGspline2 constructor)
//
// betaF[p+C*q, niter]:       sampled regression parameters for fixed effects
// betaR[pRE+C*qRE, niter]:   sampled means of random effects (if hierarCenter)
//                            * this will also contain G-spline intercepts
// ivarR_varR[???,niter]:     sampled variance components or inverse variance components for random effects
//                            * are decomposed within this procedure
//                            * For NORMAL random effects, it contains lower triangles of corresponding matrices and ??? = LT(nRandom)
//                            * For G-spline random effects, it contains only overall G-spline (inverse) variances tau^2/tau^(-2) and ??? = nRandom
// is_varR[1]:                If <> 0 then it is assumed that ivarR_varR contains variance components of the random effects
//                            and not inverse variance components
//
// skip:                     How many rows are to be skipped in the file with (log-)weights
//                           If header is included in the file, then skip should be >= 1
// by:                       Possible thinning of the file with (log-)weights
//
// logw:                     0 if weights are given directly, <>0 if log-weights are supplied
// is_indfile:               APPLICABLE only when random-effects are bivariate G-splines
//                           If <> 0 then file indpath is used
//                           Else         file indpath is ignored and it is assumed that wpath contains on each row
//                                        the weights of all (zero as well) components    
//
void
predict_cumlogit(double *prob_value,               double *prob_summary,                
                 const double *qprob,              const int *nqprob,              const int *retValue,
                 const int *C,                     const int *n,                   const int *N,                   const int *ni,
                 const double *X,                  const double *V,                const int *p,                   const int *q,
                 const double *XRE,                const double *VRE,              const int *pRE,                 const int *qRE,
                 const int *REdist,                const int *hierarCenter,
                 const int *G_dim,                 const double *G_dPar,
                 const double *betaF,              const double *betaR,            double *ivarR_varR,             const int *is_varR,
                 const int *niter,                 
                 char **wpath,                     char **indpath,                 const int *skip,                const int *by,                  
                 const int *logw,                  const int *is_indfile,          const int *nwrite,              int *err)
{
  try{
    *err = 0;

    GetRNGstate();

    const double _ZERO_ = 0.0;
    const int _ONE_ = 1;

    int j, k, ix, nwrite2;
    double *sumP, *dP;

    const char *fname = "predict_cumlogit.cpp";
    *err = 0;

    /*** Some input checks ***/
    if (*niter <= 0){
      REprintf("niter=%d\n", *niter);
      throw returnR(fname, "niter must be positive", 1);
    }
    if (*skip < 0){
      REprintf("skip=%d\n", *skip);
      throw returnR(fname, "Incorrect value of skip", 1);
    }
    if (*by <= 0 || *by > *niter){
      REprintf("niter=%d,  by=%d\n", *niter, *by);
      throw returnR(fname, "Incorrect value of by", 1);
    }
    nwrite2 = *nwrite > *niter ? *niter : *nwrite;

    /*** Some dimensionality parameters ***/
    int nFixed  = (*C)*(*q) + (*p);
    int nRandom = (*C)*(*qRE) + (*pRE);
    int LT_varR = (nRandom * (1 + nRandom))/2;
    int ngrid   = (*C+1)*(*n);                   // number of computed predictive probabilities at each iteration

    /*** Offset term (now constantly equal to zero) ***/
    double *offset = (double*) calloc(*C*(*n), sizeof(double));
    if (!offset) throw returnR(fname, "Not enough memory available (offset)", 99);  
    AK_BLAS_LAPACK::fillArray(offset, &_ZERO_, *C*(*n));

    /*** Linear predictors based on the FIXED effects ***/    
    double *etaV = (double*) calloc(*C*(*n), sizeof(double));              /*** Parts of the linear predictor based on the fixed effects not proportional w.r.t. odds ***/
    double *etaX = (double*) calloc(*n, sizeof(double));                   /*** Parts of the linear predictor based on the fixed effects proportional w.r.t. odds     ***/
    if (!etaV || !etaX) throw returnR(fname, "Not enough memory available (etaV, etaX)", 99);
    AK_BLAS_LAPACK::fillArray(etaV, &_ZERO_, *C*(*n));
    AK_BLAS_LAPACK::fillArray(etaX, &_ZERO_, *n);

    /*** Linear predictors based on the RANDOM effects ***/    
    double *etaVRE = (double*) calloc(*C*(*n), sizeof(double));              /*** Parts of the linear predictor based on the random effects not proportional w.r.t. odds ***/
    double *etaXRE = (double*) calloc(*n, sizeof(double));                   /*** Parts of the linear predictor based on the random effects proportional w.r.t. odds     ***/
    if (!etaVRE || !etaXRE) throw returnR(fname, "Not enough memory available (etaVRE, etaXRE)", 99);
    AK_BLAS_LAPACK::fillArray(etaVRE, &_ZERO_, *C*(*n));
    AK_BLAS_LAPACK::fillArray(etaXRE, &_ZERO_, *n);

    /*** Overall linear predictor ***/
    double *eta  = (double*) calloc(*C*(*n), sizeof(double));    
    if (!eta) throw returnR(fname, "Not enough memory available (eta)", 99);
    AK_BLAS_LAPACK::fillArray(eta, &_ZERO_, *C*(*n));

    /*** Reset averages ***/
    AK_BLAS_LAPACK::fillArray(prob_summary, &_ZERO_, ngrid);

    /*** Different actions for different random effects distributions ***/
    double *varRP;                                /** storage of the covariance components of the normal random effects **/
    double *bb;                                   /** storage of the values of random effects                           **/
    double *weight, *sigma;                       /** storage of G-spline related quantities                            **/    
    int *nknots, total_nknots;                    /** storage of G-spline related quantities                            **/
    const double *knots, *knots0, *knots1;        /** pointers to G-spline related quantities                           **/
    double *sigma0, *sigma1;                      /** pointers to G-spline related quantities                           **/

    int k_effect[1];
    int total_length[1];
    int *ind_w_effect;                            /** needed by models with bivariate G-splines                         **/

    std::string swpath, sindpath;                           /** needed by G-spline random effects                                 **/
    std::ifstream wfile, indfile;                           /** needed by G-spline random effects                                 **/

    const double *betaFP, *betaRP;

    double *probP  = prob_value;

    int iter;
    int backs      = 0;
    int by_1       = *by - 1;
    int jump_value = (*retValue || *nqprob) ? ngrid : 0;
    int anyZero[1];

    Rprintf("Values evaluated: ");
    switch (*REdist){      /*** distribution switch ***/

    /******* No random effects *******/
    /******* ================= *******/
    case mcmc_Random::_None:
      if (nRandom) throw returnR(fname, "Positive nRandom and _None random effects?", 1);

      betaFP = betaF;
      for (iter = 0; iter < *niter; iter++){
	util_cumlogit::prob_etas_cumlogit2(anyZero, probP, etaX, etaV, eta, offset, betaFP, X, V, *n, *C, *p, *q);

	AK_BLAS_LAPACK::a_aPlusb(prob_summary, probP, ngrid);

        betaFP += nFixed;
        probP  += jump_value;

        predict_common::print_iter_info2(backs, iter, nwrite2, *niter);
      }      
      break;

    /******* Normal random effects *******/
    /******* ===================== *******/
    case mcmc_Random::_Normal:
      bb = (double*) calloc(nRandom*(*N), sizeof(double));
      if (!bb) throw returnR(fname, "Not enough memory available (bb)", 99);  

      varRP = ivarR_varR;
      betaFP = betaF;

      if (*hierarCenter){          /*** means of random effects are in betaR ***/
        betaRP = betaR;
        for (iter = 0; iter < *niter; iter++){
	  predict_common::predict_Nrandom_effect(bb, varRP, betaRP, &nRandom, N, is_varR);
	  util_cumlogit::prob_etas_etaREs_cumlogit2(anyZero, probP, etaX, etaV, etaXRE, etaVRE, eta, offset, betaFP, bb, X, V, XRE, VRE, *N, ni, *C, *p, *q, *pRE, *qRE);

  	  AK_BLAS_LAPACK::a_aPlusb(prob_summary, probP, ngrid);

          betaFP += nFixed;
          betaRP += nRandom;
          varRP  += LT_varR;
          probP  += jump_value;

          predict_common::print_iter_info2(backs, iter, nwrite2, *niter);     
        }      
      }
      else{                        /*** random effects themselves have zero means ***/
        for (iter = 0; iter < *niter; iter++){
	  predict_common::predict_Nrandom_effect_Zero(bb, varRP, &nRandom, N, is_varR);
	  util_cumlogit::prob_etas_etaREs_cumlogit2(anyZero, probP, etaX, etaV, etaXRE, etaVRE, eta, offset, betaFP, bb, X, V, XRE, VRE, *N, ni, *C, *p, *q, *pRE, *qRE);

  	  AK_BLAS_LAPACK::a_aPlusb(prob_summary, probP, ngrid);

          betaFP += nFixed;
          varRP  += LT_varR;
          probP  += jump_value;

          predict_common::print_iter_info2(backs, iter, nwrite2, *niter);     
        }      
      }
      free(bb);
      break;

    /******* G-spline random effects *******/
    /******* ======================= *******/    
    case mcmc_Random::_Gspline:
      if (nRandom != G_dim[0]){
        REprintf("qRE=%d,  pRE=%d,  C=%d,  nRandom=%d;   G_dim[0]=%d\n", *qRE, *pRE, *C, nRandom, G_dim[0]);
        throw returnR(fname, "nRandom does not coincide with G_dim[0]", 1);
      }

      bb = (double*) calloc(nRandom*(*N), sizeof(double));
      if (!bb) throw returnR(fname, "Not enough memory available (bb)", 99);  

      nknots = (int*) calloc(nRandom, sizeof(int));
      if (!nknots) throw returnR(fname, "Not enough memory available (nknots)", 99);  
      total_nknots = 0;
      for (j = 0; j < nRandom; j++){
        nknots[j] = 2*G_dim[j+1] + 1;
        total_nknots += nknots[j];
      }

      knots = G_dPar + nRandom;

      sigma = (double*) calloc(total_nknots, sizeof(double));
      if (!sigma) throw returnR(fname, "Not enough memory available (sigma)", 99);  
      dP    = sigma;
      for (j = 0; j < nRandom; j++){
        for (k = 0; k < nknots[j]; k++){
          *dP = G_dPar[j];
          dP++;
        }
      }      
  
      /*** Open file with simulated G-spline (log-)weights and skip rows at the beginning of the file that are to be skipped ***/
      swpath = *wpath;
      In_Output::open_File_toRead(wfile, swpath, *skip);    /** Possible skipping of the header should be included in the value of skip **/

      switch (nRandom){    /*** dimension switch ***/
      case 1:
        weight    = (double*) calloc(total_nknots, sizeof(double));        
        if (!weight) throw returnR(fname, "Not enough memory available (weight)", 99);

        /*** Read the first set of G-spline weights ***/
        summary_Gspline1A::readWeightsFromFiles1(weight, 0, *skip, nknots, wfile, swpath);
        varRP = ivarR_varR;
        betaFP = betaF;

        if (*hierarCenter){          /*** means of random effects are in betaR ***/
          betaRP = betaR;

          /*** Iteration 0 ***/
          predict_common::predict_G1random_effect(bb, weight, varRP, betaRP, N, is_varR, knots, sigma, nknots, logw);
  	  util_cumlogit::prob_etas_etaREs_cumlogit2(anyZero, probP, etaX, etaV, etaXRE, etaVRE, eta, offset, betaFP, bb, X, V, XRE, VRE, *N, ni, *C, *p, *q, *pRE, *qRE);

    	  AK_BLAS_LAPACK::a_aPlusb(prob_summary, probP, ngrid);

          /*** Remaining iterations ***/
          for (iter = 1; iter < *niter; iter++){
            summary_Gspline1A::readWeightsFromFiles1(weight, by_1, iter, nknots, wfile, swpath);
            betaFP += nFixed;
            betaRP += nRandom;
            varRP++;
            probP  += jump_value;

	    predict_common::predict_G1random_effect(bb, weight, varRP, betaRP, N, is_varR, knots, sigma, nknots, logw);
  	    util_cumlogit::prob_etas_etaREs_cumlogit2(anyZero, probP, etaX, etaV, etaXRE, etaVRE, eta, offset, betaFP, bb, X, V, XRE, VRE, *N, ni, *C, *p, *q, *pRE, *qRE);

    	    AK_BLAS_LAPACK::a_aPlusb(prob_summary, probP, ngrid);

            predict_common::print_iter_info2(backs, iter, nwrite2, *niter);     
          }      
        }
        else{                        /*** random effects themselves have zero means ***/
          /*** Iteration 0 ***/
          predict_common::predict_G1random_effect_Zero(bb, weight, varRP, N, is_varR, knots, sigma, nknots, logw);
	  util_cumlogit::prob_etas_etaREs_cumlogit2(anyZero, probP, etaX, etaV, etaXRE, etaVRE, eta, offset, betaFP, bb, X, V, XRE, VRE, *N, ni, *C, *p, *q, *pRE, *qRE);

  	  AK_BLAS_LAPACK::a_aPlusb(prob_summary, probP, ngrid);

          /*** Remaining iterations ***/
          for (iter = 1; iter < *niter; iter++){
            summary_Gspline1A::readWeightsFromFiles1(weight, by_1, iter, nknots, wfile, swpath);
            betaFP += nFixed;
            varRP++;
            probP  += jump_value;

	    predict_common::predict_G1random_effect_Zero(bb, weight, varRP, N, is_varR, knots, sigma, nknots, logw);
	    util_cumlogit::prob_etas_etaREs_cumlogit2(anyZero, probP, etaX, etaV, etaXRE, etaVRE, eta, offset, betaFP, bb, X, V, XRE, VRE, *N, ni, *C, *p, *q, *pRE, *qRE);

  	    AK_BLAS_LAPACK::a_aPlusb(prob_summary, probP, ngrid);

            predict_common::print_iter_info2(backs, iter, nwrite2, *niter);     
          }      
        }

        free(weight);
        break;

      case 2:
        *total_length = nknots[0] * nknots[1];
        sigma0 = sigma;
        sigma1 = sigma + nknots[0];
        knots0 = knots;
        knots1 = knots + nknots[0];

        sindpath = *indpath;           /** Path to the file with indeces of non-zero components **/

         /*** Allocate the space for ind_w_effect and fill it by 1 ***/
        ind_w_effect = (int*) calloc(*total_length, sizeof(int));
        if (!ind_w_effect) throw returnR(fname, "Not enough memory available (ind_w_effect)", 99);
        AK_BLAS_LAPACK::fillArrayI(ind_w_effect, &_ONE_, *total_length);

        /*** Allocate the space for weights ***/
        weight    = (double*) calloc(*total_length, sizeof(double));        
        if (!weight) throw returnR(fname, "Not enough memory available (weight)", 99);

        /*** Open file with indeces of non-zero G-spline components and skip rows at the beginning of the file that are to be skipped ***/
        if (*is_indfile) In_Output::open_File_toRead(indfile, sindpath, *skip);    /** Possible skipping of the header should be included in the value of skip **/

        /*** Read the first set of G-spline weights ***/
        if (*is_indfile) summary_BiGsplineA::readWeightsFromFilesBi(k_effect, ind_w_effect, weight, 0, *skip, total_length, indfile, sindpath, wfile, swpath);
        else             summary_Gspline1A::readWeightsFromFiles1(weight, 0, *skip, total_length, wfile, swpath);
        varRP = ivarR_varR;
        betaFP = betaF;

        if (*hierarCenter){          /*** means of random effects are in betaR ***/
          betaRP = betaR;

          /*** Iteration 0 ***/
          predict_common::predict_GBirandom_effect(bb, weight, ind_w_effect, varRP, betaRP, N, is_varR, knots0, knots1, sigma0, sigma1, nknots, total_length, logw, is_indfile);
	  util_cumlogit::prob_etas_etaREs_cumlogit2(anyZero, probP, etaX, etaV, etaXRE, etaVRE, eta, offset, betaFP, bb, X, V, XRE, VRE, *N, ni, *C, *p, *q, *pRE, *qRE);

  	  AK_BLAS_LAPACK::a_aPlusb(prob_summary, probP, ngrid);

          /*** Remaining iterations ***/
          for (iter = 1; iter < *niter; iter++){
            if (*is_indfile) summary_BiGsplineA::readWeightsFromFilesBi(k_effect, ind_w_effect, weight, by_1, iter, total_length, indfile, sindpath, wfile, swpath);
            else             summary_Gspline1A::readWeightsFromFiles1(weight, by_1, iter, total_length, wfile, swpath);
            betaFP += nFixed;
            betaRP += nRandom;
            varRP++;
            probP  += jump_value;

            predict_common::predict_GBirandom_effect(bb, weight, ind_w_effect, varRP, betaRP, N, is_varR, knots0, knots1, sigma0, sigma1, nknots, total_length, logw, is_indfile);
	    util_cumlogit::prob_etas_etaREs_cumlogit2(anyZero, probP, etaX, etaV, etaXRE, etaVRE, eta, offset, betaFP, bb, X, V, XRE, VRE, *N, ni, *C, *p, *q, *pRE, *qRE);

  	    AK_BLAS_LAPACK::a_aPlusb(prob_summary, probP, ngrid);

            predict_common::print_iter_info2(backs, iter, nwrite2, *niter);     
          }      
        }
        else{                        /*** random effects themselves have zero means ***/
          /*** Iteration 0 ***/
          predict_common::predict_GBirandom_effect_Zero(bb, weight, ind_w_effect, varRP, N, is_varR, knots0, knots1, sigma0, sigma1, nknots, total_length, logw, is_indfile);
	  util_cumlogit::prob_etas_etaREs_cumlogit2(anyZero, probP, etaX, etaV, etaXRE, etaVRE, eta, offset, betaFP, bb, X, V, XRE, VRE, *N, ni, *C, *p, *q, *pRE, *qRE);

  	  AK_BLAS_LAPACK::a_aPlusb(prob_summary, probP, ngrid);

          /*** Remaining iterations ***/
          for (iter = 1; iter < *niter; iter++){
            if (*is_indfile) summary_BiGsplineA::readWeightsFromFilesBi(k_effect, ind_w_effect, weight, by_1, iter, total_length, indfile, sindpath, wfile, swpath);
            else             summary_Gspline1A::readWeightsFromFiles1(weight, by_1, iter, total_length, wfile, swpath);
            betaFP += nFixed;
            varRP++;
            probP  += jump_value;

            predict_common::predict_GBirandom_effect_Zero(bb, weight, ind_w_effect, varRP, N, is_varR, knots0, knots1, sigma0, sigma1, nknots, total_length, logw, is_indfile);
	    util_cumlogit::prob_etas_etaREs_cumlogit2(anyZero, probP, etaX, etaV, etaXRE, etaVRE, eta, offset, betaFP, bb, X, V, XRE, VRE, *N, ni, *C, *p, *q, *pRE, *qRE);

  	    AK_BLAS_LAPACK::a_aPlusb(prob_summary, probP, ngrid);

            predict_common::print_iter_info2(backs, iter, nwrite2, *niter);     
          }      
        }

        free(weight);
        free(ind_w_effect);
        break;

      default:
        REprintf("nRandom=%d\n", nRandom);
        throw returnR(fname, "Not implemented for G-spline random effects of dimension > 2", 1);
      }                    /*** end of the dimension switch ***/

      free(sigma);
      free(nknots);
      free(bb);
      break;     

    default:
      throw returnR(fname, "Unknown REdist argument", 1);
    }                      /*** end of the distribution switch ***/
    Rprintf("\n");

    /*** Cleaning ***/
    free(eta);
    free(etaVRE);
    free(etaXRE);
    free(etaV);
    free(etaX);
    free(offset);

    /*** MCMC averages ***/
    sumP = prob_summary;
    for (ix = 0; ix < ngrid; ix++){
      *sumP /= *niter;
      sumP++;
    }

    /*** Quantiles (if required) ***/
    if (*nqprob){
      Rprintf("Computing quantiles\n");
      Quantile::Quantile(sumP, prob_value, &ngrid, niter, qprob, nqprob);
    }    

    PutRNGstate();
    return;
  }
  catch(returnR rr){
    *err = rr.errflag();
    PutRNGstate();
    return;
  }
}  /*** end of predict_cumlogit() ***/

}    /*** end of extern "C" ***/

}    /*** end of namespace predict_cumlogitA ***/
