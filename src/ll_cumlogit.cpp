/*** ll_cumlogit.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  04/08/2006
//     MAJOR MODIFICATIONS:  10/10/2006
//                           (order of 'p' and 'q' parameters changed such that the Hessian has a "better" structure
//                            with respect to decompositions, see p. 49 of Rue and Held, 2005)
//
//           ll_cumlogitNR:             10/08/2006
//           ll_cumlogitFS:             18/08/2006
//
//           ll_cumlogitNR2:            18/09/2006
//           ll_cumlogitFS2:            18/09/2006
//
//           scoreNR_cumlogit:          10/10/2006
//           scoreFS_cumlogit:          10/10/2006
//           infoMatNR_cumlogit:        10/10/2006
//           infoMatFS_cumlogit:        10/10/2006
//           infoMatSDE_SDS__cumlogit:  10/10/2006
//
// PURPOSE: Cumulative logit model, log-likelihood and its derivatives
//          * minus Hessian is computed as 
//            a) ll_cumlogitNR: observed information matrix
//            b) ll_cumlogitFS: expected information matrix
//
/* ********************************************************************************* */
//
// LT(a) = number of elements in the lower triangle of the matrix a x a,
//       = 0.5*a*(a+1)
// m = C*q + p (number of unknown regression parameters)
//

#include "ll_cumlogit.h"


/************ ll_cumlogitNR ****************************************************************************/
/***                                                                                                ***/
/*** Computes the value of the log-likelihood, score and the OBSERVED information matrix            ***/
/***                                                                                                ***/
/******************************************************************************************************/
//
// ll[1]:           value of the log-likelihood
// U[1,m]:          score vector
// I[LT(m)]:        lower triangle of the observed information matrix (minus Hessian) 
//
// etaX[1,n]:       etaX[i] = x[i]'beta, i=0,...,n-1, part of the linear predictor which is proportional w.r.t. odds
// etaV[C,n]:       etaV[c,i] = v[i]'gamma[c], i=0,...,n-1, c=1,...,C,
//                    part of the linear predictor which is not proportional w.r.t. odds
//                    * in this function, only 
//                      a) y[i] \in {1,...,C-1}:
//                        etaV[y[i]-1,i] = v[i]'gamma[y[i]] and 
//                        etaV[y[i],i]   = v[i]'gamma[y[i]+1] 
//                      b) y[i] = 0:
//                        etaV[-1,i] standing for Infty
//                        etaV[0,i]  = v[i]'gamma[1] 
//                      c) y[i] = C:
//                        etaV[C-1,i] = v[i]'gamma[C] and
//                        etaV[C,i] standing for -Infty
//                      are needed
//                    * the rest would be needed for Fisher-scoring, but here can be filled arbitrarily
//                    * only needed parts will be updated by this function
//
// prob[C+1,n]:     will be needed for Fisher scoring (on input filled by predictive category probabilities)
//                    * in this function not used and not updated
//
// work[0]:         will be needed for Fisher scoring
// anyZero[0]:      will be needed for Fisher scoring
//
// offset[1, n]:    possible offset term, must not depend on c
// 
// theta[1,m]:      current vector of regression parameters (both beta and gamma)
// y[1,n]:          observations taking values from {0,...,C}
// x[p,n]:          x[i,j], covariates whose effect is proportional w.r.t. odds
// v[q,n]:          v[i,j], covariates whose effect is not proportional w.r.t odds
// xx[n][LT(p)]:    array of lower triangles carrying x[i,j]*x[i,k]
// vv[n][LT(q)]:    array of lower triangles carrying v[i,j]*v[i,k]
// xv[n][p,q]:      array of rectangles carrying x[i,j]*v[i,k]
//
// order:           0: compute only the value of the log-likelihood
//                  1: compute also the score vector
//                  2: compute also the information matrix
//
/***********************************************************************************************************/
void
ll_cumlogitNR(double *ll,          
              MatrixRect<double> *U,
              MatrixLT<double> *I,
              MatrixRect<double> *etaX,  
              MatrixRect<double> *etaV,
              MatrixRect<double> *prob,
              double *work,
              int *anyZero,
              const MatrixRect<double> *offset,
              const MatrixRect<double> *theta,
              const MatrixRect<int> *y,        
              const MatrixRect<double> *x,     
              const MatrixRect<double> *v,  
              const MatrixLT<double> *xx,  
              const MatrixLT<double> *vv,  
              const MatrixRect<double> *xv,
              const int &order)
{
  static int i, j;
  static double offset_etaX, eta0, eta1, h0, h1, h0_h1, dh0, dh1, dh0_h0_h1, dh1_h0_h1;
  static double *U_X, *I_XX, *etaX_now, *etaV_now;
  static const double *thetaP, *thetaXP, *offset_now;
  static const int *y_now;
  static const double *x_now, *v_now;
  static const MatrixLT<double> *xx_now, *vv_now;
  static const MatrixRect<double> *xv_now;

  const int C = prob->nrow()-1;          /* number of response categories minus 1 */
  const int p = x->nrow();               /* number of beta regressors             */
  const int q = v->nrow();               /* number of gamma regressors            */
  const int Cq = C*q;

  *anyZero = 0;

  /* Reset log-likelihood, score and information matrix */
  /******************************************************/
  *ll = 0.0;
  if (order > 0){
    U->fillBy(0.0);
    if (order > 1){
      I->fillBy(0.0);
    }
  }
    
  /* Loop over observations */
  /**************************/
  thetaXP = theta->aconst() + Cq;     /* first beta parameter               */
  if (p){
    U_X = U->a() + Cq;                  /* score for the first beta parameter */
    I_XX = I->a() + I->diagI(Cq);
  }
  else{
    U_X = NULL;
    I_XX = NULL;
  }

  y_now = y->aconst();
  offset_now = offset->aconst();
  etaX_now = etaX->a();
  etaV_now = etaV->a();
  x_now = x->aconst();
  v_now = v->aconst();
  xx_now = xx;
  vv_now = vv;
  xv_now = xv;

  for (i = 0; i < y->length(); i++){    /* loop observations */    

    /* Update etaX */
    /***************/
    thetaP = thetaXP;
    *etaX_now = 0.0;
    for (j = 0; j < p; j++){
      *etaX_now += (*x_now) * (*thetaP);
      x_now++;
      thetaP++;
    }
    x_now -= p;             /* return back to the beginning of the x vector for this observation */

    offset_etaX = (*offset_now) + (*etaX_now);

    /* Update etaV (needed 2 elements) and compute h0 and h1 */
    /*********************************************************/
    thetaP = theta->aconst() + (*y_now-1)*q;                          /* beginning of the first needed gamma parameters */
    etaV_now += *y_now-1;                                             /* the first needed etaV                          */ 

      /* etaV[y-1] and h0 = exp(eta[y-1])/(1 + exp(eta[y-1])) */
    if (*y_now > 0){
      *etaV_now = 0.0;      
      for (j = 0; j < q; j++){
        *etaV_now += (*v_now) * (*thetaP);
        v_now++;
        thetaP++;
      }
      v_now -= q;          /* return back to the beginning of the v vector for this observation */
      eta0 = offset_etaX + (*etaV_now);
      h0 = invlogit_AK(eta0);      
    }
    else{
      thetaP += q;
      h0 = 1.0;
    }
    etaV_now++;

      /* etaV[y] and h1 = exp(eta[y])/(1 + exp(eta[y])) */
    if (*y_now < C){
      *etaV_now = 0.0;
      for (j = 0; j < q; j++){
        *etaV_now += (*v_now) * (*thetaP);
        v_now++;
        thetaP++;
      }
      v_now -= q;                   /* return back to the beginning of the v vector for this observation */
      eta1 = offset_etaX + (*etaV_now);
      h1 = invlogit_AK(eta1);
    }
    else{
      h1 = 0.0;
    }

    offset_now++;                   /* offset of the next observation */
    etaX_now++;                     /* etaX of the next observation */
    etaV_now += C - (*y_now);       /* etaV of the next observation */
    

    /* Log-likelihood contribution */
    /*******************************/
    h0_h1 = h0 - h1;
    if (h0_h1 < _AK_ZERO){             /* probability equal to zero */
      *ll = R_NegInf;                  /* minus infinity            */
      break;
    }
    *ll += log(h0_h1);


    /* Score vector */ 
    /****************/
    if (order > 0){
      scoreNR_cumlogit(U->a(), U_X, dh0, dh1, dh0_h0_h1, dh1_h0_h1, y_now, &x_now, &v_now, p, q, C, Cq, h0, h1, h0_h1);   
               /** performs also x_now += p; v_now += q  **/


    /* Information matrix (lower triangle filled by columns) */
    /*********************************************************/
      if (order > 1){
        infoMatNR_cumlogit(I->a(), I_XX, I->diagI(), y_now, xx_now->aconst(), vv_now->aconst(), xv_now->aconst(), 
                           p, q, C, Cq, vv_now->length(), h0_h1, dh0, dh1, dh0_h0_h1, dh1_h0_h1);

        xx_now++;
        vv_now++;
        xv_now++;
      }    /* end of if (order > 1) */
    }         /* end of if (order > 0) */
    else{
      x_now += p;
      v_now += q;
    }

    y_now++;
  }                  /* end of loop observations */

  /*** Fill the super-diagonal elements of the sub-diagonal submatrices in the information matrix ***/
  if (order > 1  && q){
    infoMatSDE_SDS_cumlogit(I->a(), p, q, C);
  }

  return;
}


/************ ll_cumlogitFS ****************************************************************************/
/***                                                                                                ***/
/*** Computes the value of the log-likelihood, score and the EXPECTED information matrix            ***/
/***                                                                                                ***/
/******************************************************************************************************/
//
// ll[1]:           value of the log-likelihood
// U[1,m]:          score vector
// I[LT(m)]:        lower triangle of the observed information matrix (minus Hessian) 
//
// etaX[1,n]:       etaX[i] = x[i]'beta, i=0,...,n-1, part of the linear predictor which is proportional w.r.t. odds
// etaV[C,n]:       etaV[c,i] = v[i]'gamma[c], i=0,...,n-1, c=1,...,C,
//                    part of the linear predictor which is not proportional w.r.t. odds
//
// prob[C+1,n]:     predictive probabilities of the response categories
// work[5*C]:       working array used to store 'h[1,...,C]' values, 'dh[1,...,C]' values, 'num_xv[1,...,C]' values, etc.
//
// anyZero[1]:      equal to 1 if at least one predictive probability is equal to 0 (less than _AK_ZERO)
//                  in that case, matrix I is not computed!
// 
// offset[1, n]:    possible offset term, must not depend on c
//
// theta[1,m]:      current vector of regression parameters (both beta and gamma)
// y[1,n]:          observations taking values from {0,...,C}
// x[p,n]:          x[i,j], covariates whose effect is proportional w.r.t. odds
// v[q,n]:          v[i,j], covariates whose effect is not proportional w.r.t odds
// xx[n][LT(p)]:    array of lower triangles carrying x[i,j]*x[i,k]
// vv[n][LT(q)]:    array of lower triangles carrying v[i,j]*v[i,k]
// xv[n][p,q]:      array of rectangles carrying x[i,j]*v[i,k]
//
// order:           0: compute only the value of the log-likelihood
//                  1: compute also the score vector
//                  2: compute also the information matrix
//
/***********************************************************************************************************/
void
ll_cumlogitFS(double *ll,          
              MatrixRect<double> *U,
              MatrixLT<double> *I,
              MatrixRect<double> *etaX,  
              MatrixRect<double> *etaV,
              MatrixRect<double> *prob,
              double *work,
              int *anyZero,
              const MatrixRect<double> *offset,
              const MatrixRect<double> *theta,
              const MatrixRect<int> *y,        
              const MatrixRect<double> *x,     
              const MatrixRect<double> *v,  
              const MatrixLT<double> *xx,  
              const MatrixLT<double> *vv,  
              const MatrixRect<double> *xv,
              const int &order)
{
  static int i;
  static double h0_h1;
  static double *num_xv, *num_vv, *num_vv2;
  static double *U_X, *etaX_now, *etaV_now;
  static const double *offset_now;
  static const int *y_now;
  static const double *x_now, *v_now;
  static const MatrixLT<double> *xx_now, *vv_now;
  static const MatrixRect<double> *xv_now;
  static double *prob_now, *h_now, *dh_now;

  const int C = prob->nrow()-1;          /* number of response categories minus 1 */
  const int p = x->nrow();               /* number of beta regressors             */
  const int q = v->nrow();               /* number of gamma regressors            */
  const int Cq = C*q;

  *anyZero = 0;

  /* Reset log-likelihood, score and information matrix */
  /******************************************************/
  *ll = 0.0;
  if (order > 0){
    U->fillBy(0.0);
    if (order > 1){
      I->fillBy(0.0);
    }
  }
    
  /* Loop over observations */
  /**************************/
  U_X = U->a() + Cq;                  /* score for the first beta parameter */

  y_now = y->aconst();
  prob_now = prob->a();
  offset_now = offset->aconst();
  etaX_now = etaX->a();
  etaV_now = etaV->a();
  x_now = x->aconst();
  v_now = v->aconst();
  xx_now = xx;
  vv_now = vv;
  xv_now = xv;

  h_now = work;
  dh_now = h_now + C;
  num_xv = dh_now + C;
  num_vv = num_xv + C;
  num_vv2 = num_vv + C;

  for (i = 0; i < y->length(); i++){    /* loop observations */

    /* Update etaX, etaV, prob, h, dh */
    /**********************************/
    util_cumlogit::prob01_cumlogit(anyZero, prob_now, etaX_now, etaV_now, h_now, dh_now, offset_now, theta->aconst(), x_now, v_now, p, q, C);

    offset_now++;                       /* offset of the next observation */
    etaX_now++;                         /* etaX of the next observation */
    etaV_now += C;                      /* etaV of the next observation */
    
    /* Log-likelihood contribution */
    /*******************************/
    h0_h1 = prob_now[*y_now];
    if (h0_h1 < _AK_ZERO){             /* probability equal to zero */
      *ll = R_NegInf;                  /* minus infinity            */
      break;
    }
    *ll += log(h0_h1);


    /* Score vector */ 
    /****************/
    if (order > 0){
      scoreFS_cumlogit(U->a(), U_X, y_now, &x_now, &v_now, p, q, C, Cq, h0_h1, dh_now);    /* inside performs also x_now += p; v_now += q */


    /* Information matrix (lower triangle filled by columns) */
    /*********************************************************/
      if (order > 1 && !(*anyZero)){
        infoMatFS_cumlogit(I->a(), num_xv, num_vv, num_vv2, y_now, xx_now->aconst(), vv_now->aconst(), xv_now->aconst(),
                           p, q, C, Cq, vv_now->length(), prob_now, dh_now);

        xx_now++;
        vv_now++;
        xv_now++;
      }    /* end of if (order > 1) */
    }         /* end of if (order > 0) */
    else{
      x_now += p;
      v_now += q;
    }

    y_now++;
    prob_now += C + 1;
  }                  /* end of loop observations */

  /*** Fill the super-diagonal elements of the sub-diagonal submatrices in the information matrix ***/
  if (order > 1 && q){
    infoMatSDE_SDS_cumlogit(I->a(), p, q, C);
  }

  return;
}


/************ ll_cumlogitNR2 ****************************************************************************/
/***                                                                                                  ***/
/*** Computes the value of the log-likelihood, score and the OBSERVED information matrix,             ***/
/***  also individual contributions to log-likelihood are returned (FINALLY NOT!)                     ***/
/***                                                                                                  ***/
/***  * version that uses on many places arrays directly and not matrices                             ***/
/***  * this will be useful to compute log-likelihood based on the subset of the data only            ***/
/***  * offset term may depend on c                                                                   ***/
/***  * function also returns needed values of the linear predictor (that is etaX + etaV)             ***/
/***                                                                                                  ***/
/***    eta:  linear predictor                                                                        ***/
/***          * array of length nObs                                                                  ***/
/***          * only needed parts are updated                                                         ***/
/***                                                                                                  ***/
/*** offset:  offset term, different for each logit                                                   ***/
/***          * matrix C x nObs stored in an array                                                    ***/
/***                                                                                                  ***/
/********************************************************************************************************/
void
ll_cumlogitNR2(double *ll,
               MatrixRect<double> *U,
               MatrixLT<double> *I,
               double *etaX,  
               double *etaV,
               double *eta,
               double *prob,
               double *work,
               int *anyZero,
               const double *offset,
               const double *theta,
               const int *y,        
               const double *x,     
               const double *v,  
               const MatrixLT<double> *xx,  
               const MatrixLT<double> *vv,  
               const MatrixRect<double> *xv,
               const int &nObs,
               const int &C,                         /* number of response categories minus 1 */
               const int &p,                         /* number of beta regressors             */
               const int &q,                         /* number of gamma regressors            */
               const int &update_etas,               /* IGNORED HERE  */
               const int &order)
{
  static int i, j;
  static double eta0, h0, h1, h0_h1, dh0, dh1, dh0_h0_h1, dh1_h0_h1;
  static const double *thetaP;

  const int Cq = C*q;

  *anyZero = 0;

  /* Reset log-likelihood, score and information matrix */
  /******************************************************/
  *ll = 0.0;
  if (order > 0){
    U->fillBy(0.0);
    if (order > 1){
      I->fillBy(0.0);
    }
  }
    
  /* Loop over observations */
  /**************************/
  const double *thetaXP = theta + Cq;     /* first beta parameter               */
  double *U_X, *I_XX;
  if (p){
    U_X = U->a() + Cq;                  /* score for the first beta parameter */
    I_XX = I->a() + I->diagI(Cq);
  }
  else{
    U_X = NULL;
    I_XX = NULL;
  }

  const int *y_now = y;
  const double *offset_now = offset;
  double *etaX_now = etaX;
  double *etaV_now = etaV;
  double *eta_now = eta;
  const double *x_now = x;
  const double *v_now = v;
  const MatrixLT<double> *xx_now = xx;
  const MatrixLT<double> *vv_now = vv;
  const MatrixRect<double> *xv_now = xv;

  for (i = 0; i < nObs; i++){    /* loop observations */    

    /* Update etaX */
    /***************/
    thetaP = thetaXP;
    *etaX_now = 0.0;
    for (j = 0; j < p; j++){
      *etaX_now += (*x_now) * (*thetaP);
      x_now++;
      thetaP++;
    }
    x_now -= p;             /* return back to the beginning of the x vector for this observation */    

    /* Update etaV (needed 2 elements) and compute h0 and h1 */
    /*********************************************************/
    thetaP = theta + (*y_now-1)*q;                          /* beginning of the first needed gamma parameters */
    etaV_now += *y_now-1;                                   /* the first needed etaV                          */ 
    offset_now += *y_now-1;
    eta_now += *y_now-1;

      /* etaV[y-1] and h0 = exp(eta[y-1])/(1 + exp(eta[y-1])) */
    if (*y_now > 0){
      *etaV_now = 0.0;      
      for (j = 0; j < q; j++){
        *etaV_now += (*v_now) * (*thetaP);
        v_now++;
        thetaP++;
      }
      v_now -= q;          /* return back to the beginning of the v vector for this observation */
      *eta_now = *etaX_now + (*etaV_now);
      eta0 = *offset_now + (*eta_now);
      h0 = invlogit_AK(eta0);      
    }
    else{
      thetaP += q;
      h0 = 1.0;
    }
    etaV_now++;
    offset_now++;
    eta_now++;

      /* etaV[y] and h1 = exp(eta[y])/(1 + exp(eta[y])) */
    if (*y_now < C){
      *etaV_now = 0.0;
      for (j = 0; j < q; j++){
        *etaV_now += (*v_now) * (*thetaP);
        v_now++;
        thetaP++;
      }
      v_now -= q;                   /* return back to the beginning of the v vector for this observation */
      *eta_now = *etaX_now + (*etaV_now);
      eta0 = *offset_now + (*eta_now);
      h1 = invlogit_AK(eta0);
    }
    else{
      h1 = 0.0;
    }
    offset_now += C - (*y_now);       /* offset of the next observation */
    etaX_now++;                       /* etaX of the next observation   */
    etaV_now += C - (*y_now);         /* etaV of the next observation   */
    eta_now += C - (*y_now);


    /* Log-likelihood contribution */
    /*******************************/
    h0_h1 = h0 - h1;
    if (h0_h1 < _AK_ZERO){             /* probability equal to zero */
      *ll = R_NegInf;                  /* minus infinity            */
      break;
    }
    *ll += log(h0_h1);


    /* Score vector */ 
    /****************/
    if (order > 0){
      scoreNR_cumlogit(U->a(), U_X, dh0, dh1, dh0_h0_h1, dh1_h0_h1, y_now, &x_now, &v_now, p, q, C, Cq, h0, h1, h0_h1);   
           /** performs also x_now += p; v_now += q  **/


    /* Information matrix (lower triangle filled by columns) */
    /*********************************************************/
    /* It contains potentially many zeros!                                                     */
    /* --------------------------------------------------------------------------------------- */
      if (order > 1){
        infoMatNR_cumlogit(I->a(), I_XX, I->diagI(), y_now, xx_now->aconst(), vv_now->aconst(), xv_now->aconst(), 
                           p, q, C, Cq, vv_now->length(), h0_h1, dh0, dh1, dh0_h0_h1, dh1_h0_h1);

        xx_now++;
        vv_now++;
        xv_now++;
      }    /* end of if (order > 1) */
    }         /* end of if (order > 0) */
    else{
      x_now += p;
      v_now += q;
    }

    y_now++;
  }                  /* end of loop observations */

  /*** Fill the super-diagonal elements of the sub-diagonal submatrices in the information matrix ***/
  if (order > 1 && q){
    infoMatSDE_SDS_cumlogit(I->a(), p, q, C);
  }

  return;
}


/************ ll_cumlogitFS2 ****************************************************************************/
/***                                                                                                  ***/
/*** Computes the value of the log-likelihood, score and the EXPECTED information matrix,             ***/
/***  also individual contributions to log-likelihood are returned (FINALLY NOT!)                     ***/
/***                                                                                                  ***/
/***  * version that uses on many places arrays directly and not matrices                             ***/
/***  * this will be useful to compute log-likelihood based on the subset of the data only            ***/
/***  * offset term may depend on c                                                                   ***/
/***  * function also returns the value of the linear predictor (that is etaX + etaV)                 ***/
/***                                                                                                  ***/
/***    eta:  linear predictor                                                                        ***/
/***          * array of length nObs                                                                  ***/
/***          * only needed parts are updated                                                         ***/
/***                                                                                                  ***/
/*** offset:  offset term, different for each logit                                                   ***/
/***          * matrix C x nObs stored in an array                                                    ***/
/***                                                                                                  ***/
/********************************************************************************************************/
void
ll_cumlogitFS2(double *ll,
               MatrixRect<double> *U,
               MatrixLT<double> *I,
               double *etaX,  
               double *etaV,
               double *eta,
               double *prob,
               double *work,
               int *anyZero,
               const double *offset,
               const double *theta,
               const int *y,        
               const double *x,     
               const double *v,  
               const MatrixLT<double> *xx,  
               const MatrixLT<double> *vv,  
               const MatrixRect<double> *xv,
               const int &nObs,
               const int &C,                         /* number of response categories minus 1                        */
               const int &p,                         /* number of beta regressors                                    */
               const int &q,                         /* number of gamma regressors                                   */
               const int &update_etas,               /* 0/1 indicating whether linear predictors are to be updated   */
               const int &order)
{
  static int i;
  static double h0_h1;

  const int Cq = C*q;

  *anyZero = 0;

  /* Reset log-likelihood, score and information matrix */
  /******************************************************/
  *ll = 0.0;
  if (order > 0){
    U->fillBy(0.0);
    if (order > 1){
      I->fillBy(0.0);
    }
  }
    
  /* Update of linear predictors if needed and update of probabilities */
  /*********************************************************************/
  if (update_etas){
    util_cumlogit::prob_etas_cumlogit2(anyZero, prob, etaX, etaV, eta, offset, theta, x, v, nObs, C, p, q);
  }
  else{
    util_cumlogit::prob_cumlogit2(anyZero, prob, offset, etaX, etaV, nObs, C);
  }

  /* Loop over observations */
  /**************************/
  double *U_X = U->a() + Cq;              /* score for the first beta parameter */

  const int *y_now = y;
  const double *prob_now = prob;
  const double *x_now = x;
  const double *v_now = v;
  const MatrixLT<double> *xx_now = xx;
  const MatrixLT<double> *vv_now = vv;
  const MatrixRect<double> *xv_now = xv;
  double *h_now = work;
  double *dh_now = h_now + C;
  double *num_xv = dh_now + C;
  double *num_vv = num_xv + C;
  double *num_vv2 = num_vv + C;

  for (i = 0; i < nObs; i++){    /* loop observations */

    /* Update h, dh */
    /****************/
    util_cumlogit::prob01_cumlogit2(h_now, dh_now, prob_now, C);
    
    /* Log-likelihood contribution */
    /*******************************/
    h0_h1 = prob_now[*y_now];
    if (h0_h1 < _AK_ZERO){             /* probability equal to zero */
      *ll = R_NegInf;                  /* minus infinity            */
      break;
    }
    *ll += log(h0_h1);


    /* Score vector */ 
    /****************/
    if (order > 0){
      scoreFS_cumlogit(U->a(), U_X, y_now, &x_now, &v_now, p, q, C, Cq, h0_h1, dh_now);   /* inside performs also x_now += p; v_now += q */


    /* Information matrix (lower triangle filled by columns) */
    /*********************************************************/
      if (order > 1 && !(*anyZero)){
        infoMatFS_cumlogit(I->a(), num_xv, num_vv, num_vv2, y_now, xx_now->aconst(), vv_now->aconst(), xv_now->aconst(),
                           p, q, C, Cq, vv_now->length(), prob_now, dh_now);
     
        xx_now++;
        vv_now++;
        xv_now++;
      }    /* end of if (order > 1) */
    }         /* end of if (order > 0) */
    else{
      x_now += p;
      v_now += q;
    }

    y_now++;
    prob_now += C + 1;
  }                  /* end of loop observations */

  /*** Fill the super-diagonal elements of the sub-diagonal submatrices in the information matrix ***/
  if (order > 1 && q){
    infoMatSDE_SDS_cumlogit(I->a(), p, q, C);
  }

  return;
}


/************ scoreNR_cumlogit ************************************************************************/
/***                                                                                                ***/
/*** Computes a score contribution and add it to current score vector.                              ***/
/*** This is just a subpart of above 'NR' functions.                                                ***/
/***                                                                                                ***/
/******************************************************************************************************/
//
// U[C*q+p]:  score vector
// U_X[p]:    pointer to U[C*q] (the first beta element of the score vector)
// x[p]:      pointer to X covariates, IT IS INCREMENTED!
// v[q]:      pointer to V covariates, IT IS INCREMENTED!
void
scoreNR_cumlogit(double *U,            double *U_X,           
                 double &dh0,          double &dh1,           double &dh0_h0_h1,   double &dh1_h0_h1,
                 const int *y,         const double **x,      const double **v,
                 const int &p,         const int &q,          const int &C,        const int &Cq,
                 const double &h0,     const double &h1,      const double &h0_h1)
{
  static int j;
  static double *U_now;

  static double dh0_h0_h1_dh1_h0_h1;

  if (*y == 0){
    dh1 = h1 * (1 - h1);
    dh1_h0_h1 = dh1/h0_h1;

    /** U(1)++ **/
    U_now = U;
    for (j = 0; j < q; j++){
      *U_now -= dh1_h0_h1 * (**v);
      (*v)++;
      U_now++;
    }

    /** U(star)++ **/
    U_now = U_X;
    for (j = 0; j < p; j++){
      *U_now -= dh1_h0_h1 * (**x);
      (*x)++;
      U_now++;
    }

  }
  else{
    if (*y == C){
      dh0 = h0 * (1 - h0);
      dh0_h0_h1 = dh0/h0_h1;

      /** U(C) **/
      U_now = U + (Cq - q);
      for (j = 0; j < q; j++){
        *U_now += dh0_h0_h1 * (**v);
        (*v)++;
        U_now++;
      }

      /** U(star) **/
      for (j = 0; j < p; j++){
        *U_now += dh0_h0_h1 * (**x);
        (*x)++;
        U_now++;
      }
    }
    else{
      dh0 = h0 * (1 - h0);
      dh1 = h1 * (1 - h1);
      dh0_h0_h1 = dh0/h0_h1;
      dh1_h0_h1 = dh1/h0_h1;
      dh0_h0_h1_dh1_h0_h1 = dh0_h0_h1 - dh1_h0_h1; 

      /** U(y) and U(y+1) **/
      U_now = U + (q * (*y-1));
      for (j = 0; j < q; j++){
        *U_now += dh0_h0_h1 * (**v);
        *(U_now+q) -= dh1_h0_h1 * (**v);
        (*v)++;
        U_now++;
      }

      /** U(star) **/
      U_now = U_X;
      for (j = 0; j < p; j++){
        *U_now += dh0_h0_h1_dh1_h0_h1 * (**x);
        (*x)++;
        U_now++;
      }
    }
  }

  return;
}


/************ scoreFS_cumlogit ************************************************************************/
/***                                                                                                ***/
/*** Computes a score contribution and add it to current score vector.                              ***/
/*** This is just a subpart of above 'FS' functions.                                                ***/
/***                                                                                                ***/
/******************************************************************************************************/
//
// U[C*q+p]:  score vector
// U_X[p]:    pointer to U[C*q] (the first beta element of the score vector)
// x[p]:      pointer to X covariates, IT IS INCREMENTED!
// v[q]:      pointer to V covariates, IT IS INCREMENTED!
void
scoreFS_cumlogit(double *U,            double *U_X,
                 const int *y,         const double **x,      const double **v,
                 const int &p,         const int &q,          const int &C,      const int &Cq,
                 const double &h0_h1,  const double *dh)
{
  static int j;
  static double *U_now;
  static const double *dh1;

  static double dh1_h0_h1, dh0_h0_h1, dh0_h0_h1_dh1_h0_h1;

  if (*y == 0){
    dh1_h0_h1 = dh[0]/h0_h1;

    /** U(1)++ **/
    U_now = U;
    for (j = 0; j < q; j++){
      *U_now -= dh1_h0_h1 * (**v);
      (*v)++;
      U_now++;
    }

    /** U(star) **/
    U_now = U_X;
    for (j = 0; j < p; j++){
      *U_now -= dh1_h0_h1 * (**x);
      (*x)++;
      U_now++;
    }
  }
  else{
    if (*y == C){
      dh0_h0_h1 = dh[C-1]/h0_h1;

      /** U(C)++ **/
      U_now = U + (Cq - q);
      for (j = 0; j < q; j++){
        *U_now += dh0_h0_h1 * (**v);
        (*v)++;
        U_now++;
      }

      /** U(star) **/
      for (j = 0; j < p; j++){
        *U_now += dh0_h0_h1 * (**x);
        (*x)++;
        U_now++;
      }
    }
    else{
      dh1 = dh + (*y);
      dh0_h0_h1 = *(dh1-1)/h0_h1;
      dh1_h0_h1 = *dh1/h0_h1;
      dh0_h0_h1_dh1_h0_h1 = dh0_h0_h1 - dh1_h0_h1; 
      
      /** U(y) and U(y+1) **/
      U_now = U + (q * (*y-1));
      for (j = 0; j < q; j++){
        *U_now += dh0_h0_h1 * (**v);
        *(U_now+q) -= dh1_h0_h1 * (**v);
        (*v)++;
        U_now++;
      }

      /** U(star) **/
      U_now = U_X;      
      for (j = 0; j < p; j++){
        *U_now += dh0_h0_h1_dh1_h0_h1 * (**x);
        (*x)++;
        U_now++;
      }
    }
  }
  
  return;
}


/************ infoMatNR_cumlogit **********************************************************************/
/***                                                                                                ***/
/*** Computes a contribution to the OBSERVED information matrix and add it to its current value.    ***/
/*** This is just a subpart of above 'NR' functions.                                                ***/
/***                                                                                                ***/
/*** !!! It does not fill super-diagonal elements of the sub-diagonal submatrices !!!               ***/
/*** !!! These must be filled at the end of the loop over the observations !!!                      ***/
/***                                                                                                ***/
/******************************************************************************************************/
//
// I_XX:  pointer to the first element of I^(*,*) (xx-block)
// lvv:   length of the array vv
//
void
infoMatNR_cumlogit(double *I,             double *I_XX,      const int *diagI,
	           const int *y,          const double *xx,  const double *vv,  const double *xv,
                   const int &p,          const int &q,      const int &C,      const int &Cq,           const int &lvv,
                   const double& h0_h1,   const double &dh0, const double &dh1, const double &dh0_h0_h1, const double &dh1_h0_h1)
{
  static int j, k;
  static const double *xx_obs, *vv_obs, *xv_obs;
  static double *I_now0, *I_now10;
  static double num0, num1, num10;

  const int pq = p*q;

  xx_obs = xx;
  vv_obs = vv;
  xv_obs = xv;

  if (*y == 0){

    /* I^(1,1) and I^(*,1) */
    I_now0 = I;
    for (j = 0; j < q; j++){

      for (k = j; k < q; k++){
        *I_now0 += dh1 * (*vv_obs);
        vv_obs++;
        I_now0++;
      }
      
      I_now0 += (Cq - q);
      for (k = 0; k < p; k++){
        *I_now0 += dh1 * (*xv_obs);
        xv_obs++;
        I_now0++;
      }
    }

    /* I^(*,*) */
    I_now0 = I_XX;
    for (j = 0; j < p; j++){
      for (k = j; k < p; k++){
        *I_now0 += dh1 * (*xx_obs);
        xx_obs++;
        I_now0++;
      }
    }    
  }        /* end of if (*y == 0) */
  else{
    if (*y == C){

      /* I^(C,C) and I^(*,C) */
      I_now0 = I + diagI[Cq - q];
      for (j = 0; j < q; j++){

        for (k = j; k < q; k++){
          *I_now0 += dh0 * (*vv_obs);
          vv_obs++;
          I_now0++;
        }
      
        for (k = 0; k < p; k++){
          *I_now0 += dh0 * (*xv_obs);
          xv_obs++;
          I_now0++;
        }
      }

      /* I^(*,*) */
      for (j = 0; j < p; j++){
        for (k = j; k < p; k++){
          *I_now0 += dh0 * (*xx_obs);
          xx_obs++;
          I_now0++;
        }
      }    
    }        /* end of if (*y == C) */
    else{    /* *y in {1,...,C-1} */
      num0 = dh0 * (1 + dh1_h0_h1/h0_h1);
      num1 = dh1 * (1 + dh0_h0_h1/h0_h1);
      num10 = -dh0_h0_h1 * dh1_h0_h1;

      /* I^(y,y), I^(y+1,y) and I^(*,y) */
      I_now0 = I + diagI[(*y - 1)*q];      
      for (j = 0; j < q; j++){

        I_now10 = I_now0 + q;
        for (k = j; k < q; k++){
          *I_now0 += num0 * (*vv_obs);
          *I_now10 += num10 * (*vv_obs);
          vv_obs++;
          I_now0++;
          I_now10++;
        }
  
        I_now0 = I_now10 + ((C - (*y) - 1)*q);    /* pointer to the 1st row of the j-th column of I^(*,y) */
        for (k = 0; k < p; k++){
          *I_now0 += dh0 * (*xv_obs);
          xv_obs++;
	  I_now0++;
        }
      }
      vv_obs -= lvv;
      xv_obs -= pq;


      /* I^(y+1,y+1) and I^(*,y+1) */
      for (j = 0; j < q; j++){
        for (k = j; k < q; k++){
          *I_now0 += num1 * (*vv_obs);
          vv_obs++;
          I_now0++;
        }
  
        I_now0 += ((C - (*y) - 1)*q);    /* pointer to the 1st row of the j-th column of I^(*,y) */
        for (k = 0; k < p; k++){
          *I_now0 += dh1 * (*xv_obs);
          xv_obs++;
	  I_now0++;
        }
      }

      /* I^(*,*) */
      I_now0 = I_XX;
      for (j = 0; j < p; j++){
        for (k = j; k < p; k++){
          *I_now0 += (dh0 + dh1) * (*xx_obs);
          xx_obs++;
          I_now0++;
        }
      }
    }
  }

  return;
}


/************ infoMatFS_cumlogit **********************************************************************/
/***                                                                                                ***/
/*** Computes a contribution to the EXPECTED information matrix and add it to its current value.    ***/
/*** This is just a subpart of above 'FS' functions.                                                ***/
/***                                                                                                ***/
/*** !!! It does not fill super-diagonal elements of the sub-diagonal submatrices !!!               ***/
/*** !!! These must be filled at the end of the loop over the observations !!!                      ***/
/***                                                                                                ***/
/******************************************************************************************************/
//
// lvv:   length of the array vv
//
void
infoMatFS_cumlogit(double *I, 
                   double *num_xv,        double *num_vv,    double *num_vv2,
	           const int *y,          const double *xx,  const double *vv,  const double *xv,
                   const int &p,          const int &q,      const int &C,      const int &Cq,     const int &lvv,
                   const double *prob,    const double *dh)
{
  static int j, k, s;
  static const double *xx_obs, *vv_obs, *xv_obs, *prob1, *dh1;
  static double *I_now0, *I_now10, *num_xv_now, *num_vv_now, *num_vv2_now;
  static double num_xx;

  const int pq = p*q;

  /*** Constants for xx-block, xv-block, vv-block, vv2-block ***/
  /*** ===================================================== ***/
  dh1 = dh;                                                                 
  prob1 = prob;                                                             
  num_xx = *dh1 * (*prob1);                                

  num_xv_now = num_xv;                                     
  num_vv_now = num_vv;                                     
  num_vv2_now = num_vv2;                                                        
  for (s = 1; s < C; s++){
    *num_xv_now = *dh1;                                    
    *num_vv_now = *dh1;                                    

    dh1++;                                                                      
    prob1++;                                                                    
    num_xx += (*(dh1-1) + *dh1) * (*prob1);                

    *num_xv_now *= (*(prob1-1) + *prob1);                  
    num_xv_now++;                                          

    if (s == 1){
      *num_vv_now *= (*(prob1-1) + *prob1) + (*dh1)/(*prob1);
    }
    else{
      *num_vv_now *= *(prob1-1) + *prob1 + (*dh1)/(*prob1) + (*(dh1-2))/(*(prob1-1));
    }
    num_vv_now++;                                          

    *num_vv2_now = -(*(dh1-1))*(*dh1)/(*prob1);                                 
    num_vv2_now++;                                                              
  }

    /** s = C **/
  prob1++;                                                                      
  num_xx += *dh1 * (*prob1);                               
  *num_xv_now = *dh1 * (*(prob1-1) + *prob1);              
  *num_vv_now = *dh1 * (*(prob1-1) + *prob1 + (*(dh1-1))/(*(prob1-1)));

  
  /*** Computation of blocks of the matrix ***/
  /*** =================================== ***/

  /* I^(s,s), I^(s+1,s) and I^(*,s) */
  I_now0 = I;                      /* pointer to I^(1,1) */
  xv_obs = xv;
  vv_obs = vv;

  num_vv_now = num_vv;                                                            
  num_vv2_now = num_vv2;
  num_xv_now = num_xv;

  for (s = 1; s < C; s++){
    for (j = 0; j < q; j++){

      /* I^(s,s) and I^(s+1,s) */
      I_now10 = I_now0 + q;            /* pointer to the j-th diagonal element of I^(s+1,s) */
      for (k = j; k < q; k++){
        *I_now0 += *num_vv_now * (*vv_obs);                             
        *I_now10 += *num_vv2_now * (*vv_obs);
        vv_obs++;
        I_now0++;
        I_now10++;
      }     

      /* I^(*,s) */
      I_now0 = I_now10 + ((C - s - 1)*q);   /* pointer to the 1st row of the j-th column of I^(*,s) */
      for (k = 0; k < p; k++){
        *I_now0 += *num_xv_now * (*xv_obs);                                       
        xv_obs++;
        I_now0++;
      }
    }  
    num_vv_now++;
    num_vv2_now++;
    num_xv_now++;
    vv_obs -= lvv;
    xv_obs -= pq;
  }

  /* s = C */
  for (j = 0; j < q; j++){

    /* I^(C,C) */
    for (k = j; k < q; k++){
      *I_now0 += *num_vv_now * (*vv_obs);                             
      vv_obs++;
      I_now0++;
    }     

    /* I^(*,C) */
    for (k = 0; k < p; k++){
      *I_now0 += *num_xv_now * (*xv_obs);                                       
      xv_obs++;
      I_now0++;
    }
  }  

  /* I^(*,*) */  
  xx_obs = xx;
  for (j = 0; j < p; j++){
    for (k = j; k < p; k++){           
      *I_now0 += num_xx * (*xx_obs);                                     
      xx_obs++;
      I_now0++;                                                                      
    }
  }

  return;
}

/************ infoMatSDE_SDS_cumlogit *****************************************************************/
/***                                                                                                ***/
/*** Fill in super-diagonal elements of the sub-diagonal submatrices                                ***/
/***   by copying sub-diagonal elements of these submatrices.                                       ***/
/***                                                                                                ***/
/*** To be performed usually after the loop over the observations.                                  ***/
/***                                                                                                ***/
/*** This is just a subpart of above functions.                                                     ***/
/***                                                                                                ***/
/******************************************************************************************************/
void
infoMatSDE_SDS_cumlogit(double *I,  const int &p,  const int &q,  const int &C)
{
  static int s, j, k, skip0, skip1, skip10;
  static double *I_now10;

  I_now10 = I + q + 1;
  for (s = 1; s < C; s++){
    skip0 = q * (C - s) + p;
    skip1 = q * (C - (s+1)) + p;

    for (j = 0; j < q; j++){
      for (k = j+1; k < q; k++){
        skip10 = (q - k) + (k - j)*skip1 + k*(2*q - 2*j - k - 1)/2 + (k - j - 1)*q + j;    /* skip to the diagonal-opposite element */
        *(I_now10 + skip10) = *I_now10;
        I_now10++;
      }
      I_now10 += skip0 + 1;
    }
  }

  return;
}


