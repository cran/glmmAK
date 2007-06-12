/*** Quantile2.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  19/01/2007
//
//           Quantile:   19/01/2007
//
// PURPOSE:  Various routines to compute quantiles from the sample
//           Partially motivated by routines in quantile.cpp, predictive_GS.cpp in bayesSurv package
//
/* ********************************************************************************* */
#include "Quantile2.h"

namespace Quantile{

extern "C"{

/*** Quantile:   Compute (pointwise) quantiles for sampled functional (e.g. predictive density)   ***/
/***                                                                                              ***/
/*** ============================================================================================ ***/
//
// quantile[ngrid, nprob]:   Computed quantiles
//                             each quantile in 1 column, values for a specific grid-point in rows
// sample[ngrid, nsample]:   Sampled values of the functional,
//                             each iteration in 1 column, values for a specific grid-point in rows
// prob[nprob]:              Probabilities for quantiles we require
//
void
Quantile(double *quantile,
         const double *sample,  const int *ngrid,  const int *nsample,
         const double *prob,    const int *nprob)
{
  try{
    if (*nprob <= 0) return;

    int i, ix, j;
    double tmpd, lower;
    int *indquant1, *indquant2, *ind1P, *ind2P;
    double *value, *valP, *quantStart, *quantP;
    const double *probP, *sampStart, *sampP;   


    /*** Indeces of quantile values in sampled chain (indexing starting from 0)                    ***/
    /***    indquant1, indquant2 ..... quantile = q*sample[indquant1] + (1-q)sample[indquant2]     ***/
    indquant1  = (int*) calloc(*nprob, sizeof(int));
    indquant2  = (int*) calloc(*nprob, sizeof(int));
    if (!indquant1 || !indquant2) throw returnR("Not enough memory available in Quantile2.cpp: Quantile().", 99);

    probP = prob;
    ind1P = indquant1;
    ind2P = indquant2;
    for (i = 0; i < *nprob; i++){
      if (*probP < 0 || *probP > 1){
        REprintf("prob[%d]=%g\n", i, *probP);
        throw returnR("Error in Quantile2.cpp: Quantile(). prob must lie between 0 and 1", 1);
      }
      if (*probP <= 0){
        *ind1P = *ind2P = 0;
      }
      else{
        if (*probP >= 1){
          *ind1P = *ind2P = *nsample - 1;
        }
        else{
          tmpd = *probP * double(*nsample);
          if (fabs(tmpd - floor(tmpd + 1e-8)) < 1e-8){
            *ind2P = int(floor(tmpd + 1e-8));
            *ind1P = *ind2P - 1;
          }
          else{
            *ind1P = *ind2P = int(floor(tmpd));
          }
        }        
      }
      //Rprintf("prob[%d]=%g,  ind1=%d,  ind2=%d\n", i, *probP, *ind1P, *ind2P);
      probP++;
      ind1P++;
      ind2P++;
    }

    
    /*** Compute quantiles ***/
    value = (double*) calloc(*nsample, sizeof(double));
    if (!value) throw returnR("Not enough memory available in Quantile2.cpp: Quantile().", 99);

    sampStart  = sample;
    quantStart = quantile;
    for (ix = 0; ix < *ngrid; ix++){

      /** Copy sampled values from the ix-th grid-point **/
      valP  = value;
      sampP = sampStart;
      for (i = 0; i < *nsample; i++){
        *valP = *sampP;
        valP++;
        sampP += *ngrid;
      }
      sampStart++;

      /** Partial sorting and extracting quantiles for the ix-th grid point **/      
      probP = prob;
      ind1P = indquant1;
      ind2P = indquant2;
      quantP = quantStart;
      for (i = 0; i < *nprob; i++){
        rPsort(value, *nsample, *ind1P);                                 /***  from R/include/R-ext/Utils.h  ***/
        valP = value;
        for (j = 0; j < *ind1P; j++){
          valP++;
        }
        lower = *valP;
        valP++;

        if (*ind2P != *ind1P){
          rPsort(valP, *nsample-(*ind1P)-1, 0);                               /***  from R/include/R-ext/Utils.h  ***/
          *quantP = *probP*lower + (1 - *probP)*(*valP);
        }
        else{
          *quantP = lower;
        }

        probP++;
        ind1P++;
        ind2P++;
        quantP += *ngrid;
      }
      quantStart++;
    }


    /*** Cleaning ***/
    free(value);
    free(indquant1);
    free(indquant2);

    return;
  }
  catch(returnR rr){
    return;
  }
}

}  /*** end of extern "C" ***/

}  /*** end of the namespace Quantile ***/

