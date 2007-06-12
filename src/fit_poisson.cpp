/*** fit_poisson.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  07/02/2007
//
// PURPOSE: Poisson regression model, main fitting function
//
/* ********************************************************************************* */

#include "fit_poisson.h"

namespace Fit_poisson{

extern "C"{

// theta[p]:
// eta[n]:     OUTPUT:  values of the linear predictor at maximum
// mu[n]:      OUTPUT:  values of the expected counts at maximum
// ll[1]:
// U[p]:
// I[LT(p)]:
// Y[n]:
// offset[n]:
// X[p,n]:
// p:          Number of regression parameters
// n:          Number of observations
// niter:
// toler:
// trace:
// err:        Error flag
//

void
fit_poisson(double *theta,       double *eta,           double *mu,
            double *ll,          double *U,             double *I,
            const int *Y,        const double *offset,  const double *X,
            const int *p,        const int *n,
            int* niter,          const double *toler,   const int *trace,
            int *err)
{
  try{
    const char *fname = "fit_poisson.cpp";
    *err = 0;

    int k;
    int LTp = ((*p) * (*p + 1))/2;
    double *dP1;
    const double *cdP1;
    const int *iP1;

    /*** Matrices derived from the design matrices (needed to compute derivatives) ***/
    MatrixLT<double> *xx;
    xx = new MatrixLT<double>[*n];
    if (!xx) throw returnR(fname, "Out of memory (object 'xx')", 99);
    cdP1 = X;
    for (k = 0; k < *n; k++){
      xx[k] = MatrixLT<double>(*p, cdP1);
      cdP1 += (*p);
    }

    /*** Space for Newton-Raphson step ***/
    double *NR_step = (double*) calloc(*p, sizeof(double));
    if (!NR_step) throw returnR(fname, "Out of memory (object 'NR_step')", 99);

    /*** Working space for the information matrix ***/
    double *I_temp = (double*) calloc(LTp, sizeof(double));
    if (!I_temp) throw returnR(fname, "Out of memory (object 'I_temp')", 99);

    /*** Correction to the log-likelihood for the factorials ***/
    double *log_y_factor = (double*) calloc(*n, sizeof(double));
    if (!log_y_factor) throw returnR(fname, "Out of memory (object 'log_y_factor')", 99);    
    dP1 = log_y_factor;
    iP1 = Y;
    for (k = 0; k < *n; k++){
      *dP1 = lgamma1p(double(*iP1));
      dP1++;
      iP1++;
    }

    /*** Initial values of the log-likelihood, score and observed information matrix ***/
    Ll_poisson::ll_poisson(ll, U, I, eta, mu, offset, theta, Y, log_y_factor, X, xx, *n, *p, 2);
    
    /*** Optimization (Newton-Raphson) ***/
    AK_Optim::optim_newton_raphson03(theta, ll, U, I, I_temp, NR_step, Ll_poisson::ll_poisson, 
                                     eta, mu, offset, Y, log_y_factor, X, xx, *n, *p, 
                                     *niter, Fit_poisson::_AK_MAX_STEPHALF, *toler, *trace,
                                     Fit_poisson::_AK_MAX_PD_ATTEMPTS, Fit_poisson::_AK_EPS_PD_ATTEMPT, *err);


    /*** Clean ***/
    free(log_y_factor);
    free(I_temp);
    free(NR_step);
    delete [] xx;

  }                      /* end of try */
  catch(returnR rr){
    *err = rr.errflag() + 1000;
    return;
  }
}  /* end of fit_poisson */

}  /* end of extern "C" */

}  /* end of the namespace Fit_poisson */
