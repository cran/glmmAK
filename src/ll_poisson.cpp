/*** ll_poisson.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  07/02/2007
//
//         ll_poisson:  07/02/2007
//
// PURPOSE: Poisson regresion model, log-likelihood and its derivatives
//
/* ********************************************************************************* */
//
// LT(a) = number of elements in the lower triangle of the matrix a x a,
//       = 0.5*a*(a+1)

#include "ll_poisson.h"

namespace Ll_poisson{

/************ ll_poisson ******************************************************************************/
/***                                                                                                ***/
/*** Computes the value of the log-likelihood, score and the information matrix (OBSERVED=EXPECTED) ***/
/***                                                                                                ***/
/******************************************************************************************************/
//
// ll[1]:                Value of the log-likelihood
// U[p]:                 Value of the score vector
// I[LT(p)]:             Value of the information matrix
// eta[nObs]:            Linear predictor = log(mu)
// mu[nObs]:             Expected counts  = exp(eta)
// offset[nObs]:         Offset
// theta[p]:             Regression coefficients
// y[nObs]:              Observed counts
// log_y_factor[nObs]:   log(y!) for each observation
// x[p,nObs]:            Covariates
// xx[nObs, LT(p)]:      Array of matrices x*x'
// nObs:                 Number of observations
// p:                    Number of regression coefficients
// order:                Order of derivatives to be computed
//
void
ll_poisson(double *ll,
           double *U,
           double *I,
           double *eta,
           double *mu,
           const double *offset,
           const double *theta,
           const int *y,
           const double *log_y_factor,
           const double *x,
           const MatrixLT<double> *xx,
           const int &nObs,
           const int &p,
           const int &order)
{
  static int i, j;
  static double ll_now;
  static double *U_now, *I_now;
  static const double *thetaP, *xx_obs;

  const double _ZERO_ = 0.0;
  int LTp = (p*(p+1))/2;

  /* Reset log-likelihood, score and information matrix */
  /******************************************************/
  *ll = 0.0;
  if (order > 0){
    AK_BLAS_LAPACK::fillArray(U, &_ZERO_, p);
    if (order > 1){
      AK_BLAS_LAPACK::fillArray(I, &_ZERO_, LTp);
    }
  }


  /* Loop over observations */
  /**************************/
  const int *y_now               = y;
  const double *log_y_factor_now = log_y_factor;
  const double *offset_now       = offset;
  const double *x_now            = x;
  const MatrixLT<double> *xx_now = xx;
  double *eta_now = eta;
  double *mu_now  = mu;

  for (i = 0; i < nObs; i++){    /* loop observations */    

    /* Update eta and mu */
    /*********************/
    thetaP   = theta;
    *eta_now = 0.0;
    for (j = 0; j < p; j++){
      *eta_now += (*x_now) * (*thetaP);
      x_now++;
      thetaP++;
    }

    *mu_now  = exp_AK(*eta_now + (*offset_now));


    /* Log-likelihood contribution */
    /*******************************/
    ll_now = (*y_now)*(*eta_now) - (*mu_now) - (*log_y_factor_now);
    //Rprintf("i=%d,  y=%d,  eta=%g,  mu=%g,  ll=%g\n", i, *y_now, *eta_now, *mu_now, ll_now);
    if (ll_now <= R_NegInf){
      *ll = R_NegInf;
      break;
    }
    *ll += ll_now;


    /* Score vector */ 
    /****************/
    if (order > 0){
      x_now -= p;             /* return back to the beginning of the x vector for this observation */    
      U_now = U;

      for (j = 0; j < p; j++){
        *U_now += ((*y_now) - (*mu_now)) * (*x_now);
        U_now++;
        x_now++;
      }

    /* Information matrix (lower triangle filled by columns) */
    /*********************************************************/
      if (order > 1){
        xx_obs = xx_now->aconst();
        I_now  = I;

        for (j = 0; j < LTp; j++){
          *I_now += (*mu_now) * (*xx_obs);
          I_now++;
          xx_obs++;
	}

        xx_now++;
      }    /* end of if(order > 1) */
    }    /* end of if(order > 0) */

    eta_now++;
    mu_now++;
    offset_now++;
    y_now++;
    log_y_factor_now++;
  }    /* end of loop observations */

  return;
}

}  /** end of the namespace Ll_poisson **/
