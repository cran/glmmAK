/*** util_cumlogit.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  09/08/2006
//    MAJOR MODIFICATIONS:  10/10/2006
//                          (order of 'p' and 'q' parameters in theta vector changed)
//
//     Functions to be used with ML estimation, possible offset term may not depend on c
//     =================================================================================
//         prob_cumlogit:  09/08/2006
//       prob01_cumlogit:  17/08/2006
//                         10/10/2006, order of 'p' and 'q' parameters in theta vector changed
//
//     Functions to be used with MCMC, possible offset term may depend on c
//     ====================================================================
//        prob_cumlogit2:        19/09/2006
//        etas_cumlogit2:        19/09/2006
//                               10/10/2006, order of 'p' and 'q' parameters in theta vector changed
//   prob_etas_cumlogit2:        19/09/2006
//                               10/10/2006, order of 'p' and 'q' parameters in theta vector changed
//      prob01_cumlogit2:        19/09/2006
//  prob_etas_etaREs_cumlogit2:  30/01/2007
//
// PURPOSE: Cumulative logit model, several smaller functions related to this model
//
/* ********************************************************************************* */

#include "util_cumlogit.h"

namespace util_cumlogit{

/************ prob_cumlogit ****************************************************************************/
/*                                                                                                     */
/* Predictive category probabilities                                                                   */
/*                                                                                                     */
/* *************************************************************************************************** */
//
// prob[C+1,n]:    predictive probabilities
// offset[1,n]:    offset terms which does not depend on c
// etaX[1,n]:      linear predictor from the part of the model which is proportional w.r.t. odds
// etaV[C,n]:      linear predictor from the part of the model which is not proportional w.r.t. odds
//
void
prob_cumlogit(MatrixRect<double> *prob,
              const MatrixRect<double> *offset,
              const MatrixRect<double> *etaX,
              const MatrixRect<double> *etaV)
{
  int i, c;
  double *prob_now;
  const double *etaX_now, *etaV_now, *offset_now;
  double offset_etaX, eta, h0, h1, tprob;

  /*** Loop over observations ***/
  prob_now = prob->a();
  offset_now = offset->aconst();
  etaX_now = etaX->aconst();
  etaV_now = etaV->aconst();

  for (i = 0; i < prob->ncol(); i++){
   
    /*** prob[0,i] ***/
    offset_etaX = (*offset_now) + (*etaX_now);
    eta = offset_etaX + (*etaV_now);
    h1 = invlogit_AK(eta);
    tprob = 1 - h1;
    *prob_now = (tprob < _AK_ZERO ? 0 : tprob);
    prob_now++;
    etaV_now++;

    /*** prob[c,i], c=1,...,C-1 ***/
    for (c = 1; c < prob->nrow()-1; c++){
      h0 = h1;
      eta = offset_etaX + (*etaV_now);
      h1 = invlogit_AK(eta);
      tprob = h0 - h1;
      *prob_now = (tprob < _AK_ZERO ? 0 : tprob);
      prob_now++;
      etaV_now++;
    }
    
    /*** prob[C,i] ***/
    tprob = h1;
    *prob_now = (tprob < _AK_ZERO ? 0 : tprob);
    prob_now++;
   
    offset_now++;
    etaX_now++;
  }  /* end of the loop over observations */
  
  return;
}


/************ prob01_cumlogit *****************************************************************************/
/*                                                                                                        */
/* Predictive category probabilities for 1 observation                                                    */
/* Function also updates linear predictors and computes some quantities needed for the computation        */
/*  of the observed information matrix                                                                    */
/*                                                                                                        */
/* For the usage and understanding of the arguments see 'll_cumlogitFS' function from 'll_cumlogit.cpp'   */
/*                                                                                                        */
/* ****************************************************************************************************** */
//
// anyZero[1]:   on OUTPUT: equal to 0 if there are no probabilities equal to zero
//                          equal to 1 if at least one category probability is equal to zero (less than _AK_ZERO)
//
// prob[C+1]:     predictive probabilities
// etaX[1]:       linear predictor from the part of the model which is proportional w.r.t. odds
// etaV[C]:       linear predictor from the part of the model which is not proportional w.r.t. odds
// offset[1]:     offset that does not depend on c
// theta[C*q+p]:  regression parameters
//
void
prob01_cumlogit(int *anyZero,
                double *prob,
                double *etaX,
                double *etaV,
                double *h,
                double *dh,
                const double *offset,
                const double *theta,
                const double *x,
                const double *v,
                const int &p,
                const int &q,
                const int &C)
{
  static int j, s;

  double *etaV_now, *h_now, *dh_now, *prob_now;
  double offset_etaX, eta0, prob_sum;

  const double *thetaX = theta + C*q;  
  const double *thetaV = theta;

  const double *x_now = x;
  const double *v_now = v;

  *anyZero = 0;

  /* Update etaX */
  /***************/
  *etaX = 0.0;
  for (j = 0; j < p; j++){
    *etaX += (*x_now) * (*thetaX);
    x_now++;
    thetaX++;
  }

  /* Update etaV, prob, h and dh */
  /*******************************/
  etaV_now = etaV;
  h_now = h;
  dh_now = dh;
  prob_now = prob;

  /** prob[0], h[1], dh[1] and etaV[1] **/
  *etaV_now = 0.0;
  for (j = 0; j < q; j++){
    *etaV_now += (*v_now) * (*thetaV);
    v_now++;
    thetaV++;
  }
  v_now -= q;                                               /* return back to the beginning of the v vector */

  offset_etaX = (*offset) + (*etaX);
  eta0 = offset_etaX + (*etaV_now);
  *h_now = invlogit_AK(eta0); 
  *dh_now = *h_now * (1 - (*h_now));
  *prob_now = 1 - (*h_now);                                     /* prob[0] */
  if (*prob_now < _AK_ZERO){
    *anyZero = 1;
    *prob_now = 0;
  }
  prob_sum = (*prob_now);

  etaV_now++;
  h_now++;
  dh_now++;
  prob_now++;

  /** prob[s],h[s+1], dh[s+1] and etaV[s+1], s=1,...,C-1 **/
  for (s = 1; s < C; s++){
    *etaV_now = 0.0;
    for (j = 0; j < q; j++){
      *etaV_now += (*v_now) * (*thetaV);
      v_now++;
      thetaV++;
    }
    v_now -= q;                                          /* return back to the beginning of the v vector */

    eta0 = offset_etaX + (*etaV_now);
    *h_now = invlogit_AK(eta0); 
    *dh_now = *h_now * (1 - (*h_now));
    *prob_now = (*(h_now-1)) - (*h_now);                      /* prob[s] */
    if (*prob_now < _AK_ZERO){
      *anyZero = 1;
      *prob_now = 0;
    }
    prob_sum += (*prob_now);

    etaV_now++;
    h_now++;
    dh_now++;
    prob_now++;
  }

  /** prob[C] **/
  *prob_now = *(h_now - 1);
  if (*prob_now < _AK_ZERO){
    *anyZero = 1;
    *prob_now = 0;
  }
  prob_sum += (*prob_now);

  /*** Scale probabilities to sum up to 1 ***/
  prob_now = prob;
  for (s = 0; s <= C; s++){
    *prob_now /= prob_sum;
    prob_now++;
  }

  return;
}


/************ prob_cumlogit2 ***************************************************************************/
/*                                                                                                     */
/* Predictive category probabilities                                                                   */
/* * offset term may depend of c                                                                       */
/*                                                                                                     */
/* *************************************************************************************************** */
//
// anyZero[1]:     indicator whether any of the predictive probabilities is equal to zero
// prob[C+1,n]:    predictive probabilities
// offset[C,n]:    offset terms that may depend on c
// etaX[1,n]:      linear predictor from the part of the model which is proportional w.r.t. odds
// etaV[C,n]:      linear predictor from the part of the model which is not proportional w.r.t. odds
//
void
prob_cumlogit2(int *anyZero,
               double *prob,
               const double *offset,
               const double *etaX,
               const double *etaV,
               const int &nObs,
               const int &C)
{
  int i, c;
  double eta0, h0, h1, tprob;

  *anyZero = 0;

  /*** Loop over observations ***/
  const double *offset_now = offset;
  const double *etaX_now = etaX;
  const double *etaV_now = etaV;
  double *prob_now = prob;

  for (i = 0; i < nObs; i++){
   
    /*** prob[0,i] ***/
    eta0 = *offset_now + (*etaX_now) + (*etaV_now);
    h1 = invlogit_AK(eta0);
    tprob = 1 - h1;
    if (tprob < _AK_ZERO){
      *prob_now = 0;
      *anyZero = 1;
    }
    else{
      *prob_now = tprob;
    }
    prob_now++;
    etaV_now++;
    offset_now++;

    /*** prob[c,i], c=1,...,C-1 ***/
    for (c = 1; c < C; c++){
      h0 = h1;
      eta0 = *offset_now + (*etaX_now) + (*etaV_now);
      h1 = invlogit_AK(eta0);
      tprob = h0 - h1;
      if (tprob < _AK_ZERO){
        *prob_now = 0;
        *anyZero = 1;
      }
      else{
        *prob_now = tprob;
      }
      prob_now++;
      etaV_now++;
      offset_now++;
    }
    
    /*** prob[C,i] ***/
    tprob = h1;
    if (tprob < _AK_ZERO){
      *prob_now = 0;
      *anyZero = 1;
    }
    else{
      *prob_now = tprob;
    }
    prob_now++;
   
    etaX_now++;
  }  /* end of the loop over observations */
  
  return;
}


/************ etas_cumlogit2 ***************************************************************************/
/*                                                                                                     */
/* Update linear predictors                                                                            */
/*                                                                                                     */
/* *************************************************************************************************** */
//
// etaX[nObs]:       linear predictor from the part of the model which is proportional w.r.t. odds
// etaV[C,nObs]:     linear predictor from the part of the model which is not proportional w.r.t. odds
// eta[C,nObs]:      linear predictor without the offset term, that is eta[c] = etaX + etaV[c]
// theta[C*q+p]:     regression parameters
// x[p,nObs]:        covariates for effects proportional w.r.t. odds
// v[q,nObs]:        covariates for effects not proportional w.r.t. odds
//
void
etas_cumlogit2(double *etaX,
               double *etaV,
               double *eta,
               const double *theta,
               const double *x,
               const double *v,
               const int &nObs,
               const int &C,
               const int &p,
               const int &q)
{
  static int i, j, s;

  const double *thetaX = theta + C*q;

  const double *x_now = x;
  const double *v_now = v;
  const double *thetaP;
  double *eta_now = eta;
  double *etaX_now = etaX;
  double *etaV_now = etaV;  

  for (i = 0; i < nObs; i++){

    /* Update etaX */
    /***************/
    thetaP = thetaX;
    *etaX_now = 0.0;
    for (j = 0; j < p; j++){
      *etaX_now += (*x_now) * (*thetaP);
      x_now++;
      thetaP++;
    }

    /* Update etaV, eta  */
    /*********************/
    thetaP = theta;
    for (s = 0; s < C; s++){
      *etaV_now = 0.0;
      for (j = 0; j < q; j++){
        *etaV_now += (*v_now) * (*thetaP);
        v_now++;
        thetaP++;
      }
      v_now -= q;                                          /* return back to the beginning of the v vector */

      *eta_now = *etaX_now + (*etaV_now);
      eta_now++;
      etaV_now++;
    }

    etaX_now++;
    v_now += q;
  }

  return;
}


/************ prob_etas_cumlogit2 *************************************************************************/
/*                                                                                                        */
/* Predictive category probabilities for all observations                                                 */
/* Function also updates linear predictors and computes some quantities needed for the computation        */
/*  of the observed information matrix                                                                    */
/* Offset term is allowed to depend on c                                                                  */
/*                                                                                                        */
/* ****************************************************************************************************** */
//
// anyZero[1]:       indicator whether any of the predictive probabilities is equal to zero
// prob[C+1,nObs]:   predictive probabilities
// etaX[nObs]:       linear predictor from the part of the model which is proportional w.r.t. odds
// etaV[C,nObs]:     linear predictor from the part of the model which is not proportional w.r.t. odds
// eta[C,nObs]:      linear predictor without the offset term, that is eta[c] = etaX + etaV[c]
// offset[C,nObs]:   offset that may depend on c
// theta[C*q+p]:     regression parameters
// x[p,nObs]:        covariates for effects proportional w.r.t. odds
// v[q,nObs]:        covariates for effects not proportional w.r.t. odds
//
void
prob_etas_cumlogit2(int *anyZero,
                    double *prob,
                    double *etaX,
                    double *etaV,
                    double *eta,
                    const double *offset,
                    const double *theta,
                    const double *x,
                    const double *v,
                    const int &nObs,
                    const int &C,
                    const int &p,
                    const int &q)
{
  static int i, j, s;

  double prob_sum, eta0, h0, h1;

  *anyZero = 0;

  const double *thetaX = theta + C*q;

  const double *x_now = x;
  const double *v_now = v;
  const double *offset_now = offset;
  const double *thetaP;
  double *prob_now = prob;
  double *eta_now = eta;
  double *etaX_now = etaX;
  double *etaV_now = etaV;  

  for (i = 0; i < nObs; i++){

    /* Update etaX */
    /***************/
    thetaP = thetaX;
    *etaX_now = 0.0;
    for (j = 0; j < p; j++){
      *etaX_now += (*x_now) * (*thetaP);
      x_now++;
      thetaP++;
    }

    /* Update etaV, eta, prob, h and dh */
    /************************************/

    /** prob[0], etaV[1] **/
    thetaP = theta;
    *etaV_now = 0.0;
    for (j = 0; j < q; j++){
      *etaV_now += (*v_now) * (*thetaP);
      v_now++;
      thetaP++;
    }

    *eta_now = *etaX_now + (*etaV_now);
    eta0 = *offset_now + (*eta_now);
    h1 = invlogit_AK(eta0); 
    *prob_now = 1 - h1;                                       /* prob[0] */
    if (*prob_now < _AK_ZERO){
      *prob_now = 0;
      *anyZero = 1;
    }
    prob_sum = (*prob_now);

    eta_now++;
    etaV_now++;
    offset_now++;
    prob_now++;

    /** prob[s], etaV[s+1], s=1,...,C-1 **/
    for (s = 1; s < C; s++){
      v_now -= q;                                          /* return back to the beginning of the v vector */
      h0 = h1;
      *etaV_now = 0.0;
      for (j = 0; j < q; j++){
        *etaV_now += (*v_now) * (*thetaP);
        v_now++;
        thetaP++;
      }


      *eta_now = *etaX_now + (*etaV_now);
      eta0 = *offset_now + (*eta_now);
      h1 = invlogit_AK(eta0); 
      *prob_now = h0 - h1;                                /* prob[s] */
      if (*prob_now < _AK_ZERO){
        *prob_now = 0;
        *anyZero = 1;
      }
      prob_sum += (*prob_now);

      eta_now++;
      etaV_now++;
      offset_now++;
      prob_now++;
    }

    /** prob[C] **/
    *prob_now = h1;
    if (*prob_now < _AK_ZERO){
      *prob_now = 0;
      *anyZero = 1;
    }
    prob_sum += (*prob_now);

    /*** Scale probabilities to sum up to 1 ***/
    prob_now -= C;
    for (s = 0; s <= C; s++){
      *prob_now /= prob_sum;
      prob_now++;
    }

    etaX_now++;
  }

  return;
}


/************ prob01_cumlogit2 ****************************************************************************/
/*                                                                                                        */
/* For 1 observation, compute quantities needed for derivatives                                           */
/*                                                                                                        */
/* ****************************************************************************************************** */
//
// h[C]:
// dh[C]:
// prob[C+1]:    predictive probabilities
//
void
prob01_cumlogit2(double *h,
                 double *dh,
                 const double *prob,
                 const int &C)
{
  static int s;

  double *h_now = h;
  double *dh_now = dh;
  const double *prob_now = prob;

  /** prob[0], h[1], dh[1] **/
  *h_now = 1 - (*prob_now);
  *dh_now = *h_now * (1 - (*h_now));
  h_now++;
  dh_now++;
  prob_now++;

  /** prob[s], h[s+1], dh[s+1], s=1,...,C-1 **/
  for (s = 1; s < C; s++){
    *h_now = *(h_now-1) - (*prob_now);
    *dh_now = *h_now * (1 - (*h_now));
    h_now++;
    dh_now++;
    prob_now++;
  }

  return;
}


/************ prob_etas_etaREs_cumlogit2 *************************************************************************/
/*                                                                                                               */
/* Predictive category probabilities for all observations in the model with random-effects                       */
/* Function also updates linear predictors based on both fixed-effects and random-effects                        */
/*                                                                                                               */
/* Offset term is allowed to depend on c                                                                         */
/*                                                                                                               */
/* ************************************************************************************************************* */
//
// anyZero[1]:                    indicator whether any of the predictive probabilities is equal to zero
// prob[C+1,nObs]:                predictive probabilities
// etaX[nObs]:                    linear predictor from the part of the model which is proportional w.r.t. odds
// etaV[C,nObs]:                  linear predictor from the part of the model which is not proportional w.r.t. odds
// etaXRE[nObs]:                  linear predictor from the random-effect part of the model which is proportional w.r.t. odds
// etaVRE[C,nObs]:                linear predictor from the random-effect part of the model which is not proportional w.r.t. odds
// eta[C,nObs]:                   linear predictor without the offset term, that is eta[c] = etaX + etaV[c]
// offset[C,nObs]:                offset that may depend on c
//
// theta[C*q+p]:                  regression parameters (fixed-effects)
// thetaRE[C*qRE+pRE,nCluster]:   values of random-effects   
//
// x[p,nObs]:                     covariates for fixed-effects proportional w.r.t. odds
// v[q,nObs]:                     covariates for fixed-effects not proportional w.r.t. odds
// xRE[pRE,nObs]:                 covariates for random-effects proportional w.r.t. odds
// vRE[qRE,nObs]:                 covariates for random-effects not proportional w.r.t. odds
//
void
prob_etas_etaREs_cumlogit2(int *anyZero,
                           double *prob,
                           double *etaX,
                           double *etaV,
                           double *etaXRE,
                           double *etaVRE,
                           double *eta,
                           const double *offset,
                           const double *theta,
                           const double *thetaRE,
                           const double *x,
                           const double *v,
                           const double *xRE,
                           const double *vRE,
                           const int &nCluster,
                           const int *ni,
                           const int &C,
                           const int &p,
                           const int &q,
                           const int &pRE,
                           const int &qRE)
{
  static int i, k, j, s;
  static int C_qRE;
  C_qRE = C*qRE;

  double prob_sum, eta0, h0, h1;

  *anyZero = 0;

  const double *thetaX = theta + C*q;
  const double *thetaXRE;

  const double *x_now       = x;
  const double *v_now       = v;
  const double *xRE_now     = xRE;
  const double *vRE_now     = vRE;
  const double *offset_now  = offset;
  const double *thetaRE_now = thetaRE;
  const double *thetaP;
  const double *thetaREP;
  double *prob_now   = prob;
  double *eta_now    = eta;
  double *etaX_now   = etaX;
  double *etaV_now   = etaV;  
  double *etaXRE_now = etaXRE;
  double *etaVRE_now = etaVRE;  


  const int *niP = ni;

  for (i = 0; i < nCluster; i++){
    thetaXRE = thetaRE_now + C_qRE;
    thetaREP = thetaXRE;

    for (k = 0; k < *niP; k++){

      /* Update etaX */
      /***************/
      thetaP    = thetaX;
      *etaX_now = 0.0;
      for (j = 0; j < p; j++){
        *etaX_now += (*x_now) * (*thetaP);
        x_now++;
        thetaP++;
      }

      /* Update etaXRE */
      /*****************/
      thetaREP    = thetaXRE;
      *etaXRE_now = 0.0;
      for (j = 0; j < pRE; j++){
        *etaXRE_now += (*xRE_now) * (*thetaREP);
        xRE_now++;
        thetaREP++;
      }


      /* Update etaV, etaVRE, eta, prob  */
      /***********************************/

      /** prob[0], etaV[1], etaVRE[1] **/
      thetaP    = theta;
      *etaV_now = 0.0;
      for (j = 0; j < q; j++){
        *etaV_now += (*v_now) * (*thetaP);
        v_now++;
        thetaP++;
      }

      thetaREP    = thetaRE_now;
      *etaVRE_now = 0.0;
      for (j = 0; j < qRE; j++){
        *etaVRE_now += (*vRE_now) * (*thetaREP);
        vRE_now++;
        thetaREP++;
      }

      *eta_now  = *etaX_now + (*etaV_now) + (*etaXRE_now) + (*etaVRE_now);
      eta0      = *offset_now + (*eta_now);
      h1        = invlogit_AK(eta0); 
      *prob_now = 1 - h1;                                       /* prob[0] */
      if (*prob_now < _AK_ZERO){
        *prob_now = 0;
        *anyZero  = 1;
      }
      prob_sum = (*prob_now);

      eta_now++;
      etaV_now++;
      etaVRE_now++;
      offset_now++;
      prob_now++;

      /** prob[s], etaV[s+1], s=1,...,C-1 **/
      for (s = 1; s < C; s++){
        v_now   -= q;                                          /* return back to the beginning of the v vector */
        *etaV_now = 0.0;
        for (j = 0; j < q; j++){
          *etaV_now += (*v_now) * (*thetaP);
          v_now++;
          thetaP++;
        }

        vRE_now -= q;                                          /* return back to the beginning of the v vector */
        *etaVRE_now = 0.0;
        for (j = 0; j < qRE; j++){
          *etaVRE_now += (*vRE_now) * (*thetaREP);
          vRE_now++;
          thetaREP++;
        }

        h0        = h1;
        *eta_now  = *etaX_now + (*etaV_now) + (*etaXRE_now) + (*etaVRE_now);
        eta0      = *offset_now + (*eta_now);
        h1        = invlogit_AK(eta0); 
        *prob_now = h0 - h1;                                /* prob[s] */
        if (*prob_now < _AK_ZERO){
          *prob_now = 0;
          *anyZero  = 1;
        }
        prob_sum += (*prob_now);

        eta_now++;
        etaV_now++;
	etaVRE_now++;
        offset_now++;
        prob_now++;
      }

      /** prob[C] **/
      *prob_now = h1;
      if (*prob_now < _AK_ZERO){
        *prob_now = 0;
        *anyZero  = 1;
      }
      prob_sum += (*prob_now);

      /*** Scale probabilities to sum up to 1 ***/
      prob_now -= C;
      for (s = 0; s <= C; s++){
        *prob_now /= prob_sum;
        prob_now++;
      }

      etaX_now++;
      etaXRE_now++;
    }

    thetaRE_now = thetaREP + pRE;
    niP++;
  }

  return;
}


}    /*** end of the namespace util_cumlogit ***/
