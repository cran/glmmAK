/*** predict_poisson.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:                           16/02/2007
//    VERSION WITH BIVARIATE G-SPLINES:  12/04/2007
//
//
// PURPOSE: Poisson log-linear model with random effects
//          * predictive counts
//
/* ********************************************************************************* */

#include "predict_poisson.h"

namespace predict_poissonA{

extern "C"{

// ecount_value:             If (retValue || nqprob) ecount_value[n, niter] = values of all evaluated predictive counts
//                           Else                    ecount_value[n]        = values of the last evaluated predictive counts
//                           
// ecount_summary[(C+1)*n, 1+nqprob]:  Pointwise averages of the predictive counts
//                                     + quantiles (if nqprob>0)
// retValue:                         Do we want to return the values of the predictive probabilities at each iteration?
//
// qprob[nprob]:                Probabilities of required quantiles
// nqprob:                      Number of quantiles required
//
// offset[n]:                   Possible offset term
// X[p,n]:                      Covariates for fixed effects
// p[1]:                        Number of fixed effects
//
// G_dim[]:              Number of random effects and numbers of knots on each side of the zero knot for each margin of the G-spline
//                       (see argument dimPar of the parametrix Gspline2 constructor)
//
// G_dpar[]:             Basis standard deviations and knots
//                       (see argument penMix_dPar of the parametric Gspline2 constructor)
//
// betaF[p, niter]:           sampled regression parameters for fixed effects
// betaR[pRE, niter]:         sampled means of random effects (if hierarCenter)
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
predict_poisson(double *ecount_value,             double *ecount_summary,                
                const double *qprob,              const int *nqprob,              const int *retValue,
                const int *n,                     const int *N,                   const int *ni,
                const double *offset,             const double *X,                const int *p,
                const double *XRE,                const int *pRE,
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

    const char *fname = "predict_poisson.cpp";
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
    int nFixed  = *p;
    int nRandom = *pRE;
    int LT_varR = (nRandom * (1 + nRandom))/2;
    int ngrid   = *n;                                 // number of computed predictive counts at each iteration

    /*** Linear predictors based on the FIXED effects ***/    
    double *etaF = (double*) calloc(*n, sizeof(double));
    if (!etaF) throw returnR(fname, "Not enough memory available (etaF)", 99);
    AK_BLAS_LAPACK::fillArray(etaF, &_ZERO_, *n);

    /*** Linear predictors based on the RANDOM effects ***/    
    double *etaRE = (double*) calloc(*n, sizeof(double));
    if (!etaRE) throw returnR(fname, "Not enough memory available (etaRE)", 99);
    AK_BLAS_LAPACK::fillArray(etaRE, &_ZERO_, *n);

    /*** Overall linear predictor ***/
    double *eta  = (double*) calloc(*n, sizeof(double));    
    if (!eta) throw returnR(fname, "Not enough memory available (eta)", 99);
    AK_BLAS_LAPACK::fillArray(eta, &_ZERO_, *n);

    /*** Reset averages ***/
    AK_BLAS_LAPACK::fillArray(ecount_summary, &_ZERO_, ngrid);

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

    double *ecountP  = ecount_value;

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
	util_poisson::ecount_etas_poisson(ecountP, eta, offset, betaFP, X, *n, *p);

	AK_BLAS_LAPACK::a_aPlusb(ecount_summary, ecountP, ngrid);

        betaFP  += nFixed;
        ecountP += jump_value;

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
	  util_poisson::ecount_etas_etaREs_poisson(ecountP, etaF, etaRE, eta, offset, betaFP, bb, X, XRE, *N, ni, *p, *pRE);

  	  AK_BLAS_LAPACK::a_aPlusb(ecount_summary, ecountP, ngrid);

          betaFP += nFixed;
          betaRP += nRandom;
          varRP  += LT_varR;
          ecountP  += jump_value;

          predict_common::print_iter_info2(backs, iter, nwrite2, *niter);     
        }      
      }
      else{                        /*** random effects themselves have zero means ***/
        for (iter = 0; iter < *niter; iter++){
	  predict_common::predict_Nrandom_effect_Zero(bb, varRP, &nRandom, N, is_varR);
	  util_poisson::ecount_etas_etaREs_poisson(ecountP, etaF, etaRE, eta, offset, betaFP, bb, X, XRE, *N, ni, *p, *pRE);

  	  AK_BLAS_LAPACK::a_aPlusb(ecount_summary, ecountP, ngrid);

          betaFP += nFixed;
          varRP  += LT_varR;
          ecountP  += jump_value;

          predict_common::print_iter_info2(backs, iter, nwrite2, *niter);     
        }      
      }
      free(bb);
      break;

    /******* G-spline random effects *******/
    /******* ======================= *******/    
    case mcmc_Random::_Gspline:
      if (nRandom != G_dim[0]){
        REprintf("pRE=%d,  nRandom=%d;   G_dim[0]=%d\n", *pRE, nRandom, G_dim[0]);
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
	  util_poisson::ecount_etas_etaREs_poisson(ecountP, etaF, etaRE, eta, offset, betaFP, bb, X, XRE, *N, ni, *p, *pRE);

    	  AK_BLAS_LAPACK::a_aPlusb(ecount_summary, ecountP, ngrid);

          /*** Remaining iterations ***/
          for (iter = 1; iter < *niter; iter++){
            summary_Gspline1A::readWeightsFromFiles1(weight, by_1, iter, nknots, wfile, swpath);
            betaFP += nFixed;
            betaRP += nRandom;
            varRP++;
            ecountP += jump_value;

	    predict_common::predict_G1random_effect(bb, weight, varRP, betaRP, N, is_varR, knots, sigma, nknots, logw);
  	    util_poisson::ecount_etas_etaREs_poisson(ecountP, etaF, etaRE, eta, offset, betaFP, bb, X, XRE, *N, ni, *p, *pRE);

    	    AK_BLAS_LAPACK::a_aPlusb(ecount_summary, ecountP, ngrid);

            predict_common::print_iter_info2(backs, iter, nwrite2, *niter);     
          }      
        }
        else{                        /*** random effects themselves have zero means ***/
          /*** Iteration 0 ***/
          predict_common::predict_G1random_effect_Zero(bb, weight, varRP, N, is_varR, knots, sigma, nknots, logw);
          util_poisson::ecount_etas_etaREs_poisson(ecountP, etaF, etaRE, eta, offset, betaFP, bb, X, XRE, *N, ni, *p, *pRE);

  	  AK_BLAS_LAPACK::a_aPlusb(ecount_summary, ecountP, ngrid);

          /*** Remaining iterations ***/
          for (iter = 1; iter < *niter; iter++){
            summary_Gspline1A::readWeightsFromFiles1(weight, by_1, iter, nknots, wfile, swpath);
            betaFP += nFixed;
            varRP++;
            ecountP += jump_value;

	    predict_common::predict_G1random_effect_Zero(bb, weight, varRP, N, is_varR, knots, sigma, nknots, logw);
            util_poisson::ecount_etas_etaREs_poisson(ecountP, etaF, etaRE, eta, offset, betaFP, bb, X, XRE, *N, ni, *p, *pRE);

  	    AK_BLAS_LAPACK::a_aPlusb(ecount_summary, ecountP, ngrid);

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
	  util_poisson::ecount_etas_etaREs_poisson(ecountP, etaF, etaRE, eta, offset, betaFP, bb, X, XRE, *N, ni, *p, *pRE);

    	  AK_BLAS_LAPACK::a_aPlusb(ecount_summary, ecountP, ngrid);

          /*** Remaining iterations ***/
          for (iter = 1; iter < *niter; iter++){
            if (*is_indfile) summary_BiGsplineA::readWeightsFromFilesBi(k_effect, ind_w_effect, weight, by_1, iter, total_length, indfile, sindpath, wfile, swpath);
            else             summary_Gspline1A::readWeightsFromFiles1(weight, by_1, iter, total_length, wfile, swpath);
            betaFP += nFixed;
            betaRP += nRandom;
            varRP++;
            ecountP += jump_value;

            predict_common::predict_GBirandom_effect(bb, weight, ind_w_effect, varRP, betaRP, N, is_varR, knots0, knots1, sigma0, sigma1, nknots, total_length, logw, is_indfile);
  	    util_poisson::ecount_etas_etaREs_poisson(ecountP, etaF, etaRE, eta, offset, betaFP, bb, X, XRE, *N, ni, *p, *pRE);

    	    AK_BLAS_LAPACK::a_aPlusb(ecount_summary, ecountP, ngrid);

            predict_common::print_iter_info2(backs, iter, nwrite2, *niter);     
          }      
        }
        else{                        /*** random effects themselves have zero means ***/
          /*** Iteration 0 ***/
          predict_common::predict_GBirandom_effect_Zero(bb, weight, ind_w_effect, varRP, N, is_varR, knots0, knots1, sigma0, sigma1, nknots, total_length, logw, is_indfile);
          util_poisson::ecount_etas_etaREs_poisson(ecountP, etaF, etaRE, eta, offset, betaFP, bb, X, XRE, *N, ni, *p, *pRE);

  	  AK_BLAS_LAPACK::a_aPlusb(ecount_summary, ecountP, ngrid);

          /*** Remaining iterations ***/
          for (iter = 1; iter < *niter; iter++){
            if (*is_indfile) summary_BiGsplineA::readWeightsFromFilesBi(k_effect, ind_w_effect, weight, by_1, iter, total_length, indfile, sindpath, wfile, swpath);
            else             summary_Gspline1A::readWeightsFromFiles1(weight, by_1, iter, total_length, wfile, swpath);
            betaFP += nFixed;
            varRP++;
            ecountP += jump_value;

            predict_common::predict_GBirandom_effect_Zero(bb, weight, ind_w_effect, varRP, N, is_varR, knots0, knots1, sigma0, sigma1, nknots, total_length, logw, is_indfile);
            util_poisson::ecount_etas_etaREs_poisson(ecountP, etaF, etaRE, eta, offset, betaFP, bb, X, XRE, *N, ni, *p, *pRE);

  	    AK_BLAS_LAPACK::a_aPlusb(ecount_summary, ecountP, ngrid);

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
    free(etaF);
    free(etaRE);
    free(eta);

    /*** MCMC averages ***/
    sumP = ecount_summary;
    for (ix = 0; ix < ngrid; ix++){
      *sumP /= *niter;
      sumP++;
    }

    /*** Quantiles (if required) ***/
    if (*nqprob){
      Rprintf("Computing quantiles\n");
      Quantile::Quantile(sumP, ecount_value, &ngrid, niter, qprob, nqprob);
    }    

    PutRNGstate();
    return;
  }
  catch(returnR rr){
    *err = rr.errflag();
    PutRNGstate();
    return;
  }
}  /*** end of predict_poisson() ***/

}    /*** end of extern "C" ***/

}    /*** end of namespace predict_poissonA ***/
