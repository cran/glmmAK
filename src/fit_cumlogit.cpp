/*** fit_cumlogit.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  07/08/2006
//      FIRST WORKING VERSION: 18/08/2006
//      MAJOR MODIFICATIONS:   11/10/2006
//                            (order of 'p' and 'q' parameters changed such that the Hessian has a "better" structure
//                            with respect to decompositions, see p. 49 of Rue and Held, 2005)
//
//
// PURPOSE: Cumulative logit model, main fitting function
//
/* ********************************************************************************* */

#include "fit_cumlogit.h"

namespace Fit_cumlogit{

extern "C"{

// theta[C*q+p]:        regression coefficients
// ll[1]:               value of the log-likelihood
// U[C*q+p]:            score vector
// I_obs[(C*q+p)^2]:    observed information matrix (minus Hessian) 
//                        * full matrix in column major order
//                      OUTPUT: inverse of the observed information matrix
//
// I_exp[(C*q+p)^2]:    expected information matrix
//                        * full matrix in column major order
//                      OUTPUT: inverse of the expected information matrix
//
// Y[n]:                observations taking values from {0,...,C}
// X[p*n]:              x[i,j], covariates whose effect is proportional w.r.t. odds
//                        * X matrix from R must be first TRANSPOSED (in R) before passing to C++!!!
// V[q*n]:              v[i,j], covariates whose effect is not proportional w.r.t odds
//                        * V matrix from R must be first TRANSPOSED (in R) before passing to C++!!!
//
// C[1]:                number of response categories minus 1, i.e., response comes from {0,...,C}
// p[1]:                number of covariates proportional w.r.t. odds
// q[1]:                number of covariates not proportional w.r.t. odds (including intercepts)
// n[1]:                number of observations
//
// niter[1]:            INPUT: maximal number of iterations of the optimizer
//                        * if equal to 0, only log-likelihood, its derivatives and inverse of the information matrix are computed 
//                      OUTPUT: number of iterations performed
// toler[1]:            tolerance to detect convergence
//
// err[1]:              error flag

void
fit_cumlogit(double *theta,
             double *ll,          double *U,           double *I_obs,      double *I_exp,
             const int *Y,        const double *X,     const double *V,  
             const int *C,        const int *p,        const int *q,       const int *n,
             int* niter,          const double *toler, const int *trace,
             int *err)
{
  try{
    const char *fname = "fit_cumlogit.cpp";
    *err = 0;

    int k;
    int itemp;
    const double *cdP1, *cdP2;

    const int Cq = (*C) * (*q);
    const int nBG = Cq + (*p);        /* number of unknown regression parameters */

    /*** Function to evaluate log-likelihood and compute derivatives ***/
    void (*eval_ll)(double*, MatrixRect<double>*, MatrixLT<double>*, 
                    MatrixRect<double>*, MatrixRect<double>*, MatrixRect<double>*, double*, int*,
                    const MatrixRect<double>*, const MatrixRect<double>*, const MatrixRect<int>*, 
                    const MatrixRect<double>*, const MatrixRect<double>*,
                    const MatrixLT<double>*, const MatrixLT<double>*,
                    const MatrixRect<double>*, const int&);

    const int algor = Fit_cumlogit::NR;
    switch (algor){
    case Fit_cumlogit::NR:
      eval_ll = ll_cumlogitNR;
      break;
    case Fit_cumlogit::FS:
      eval_ll = ll_cumlogitFS;
      break;
    default:
      throw returnR(fname, "Unknown argument 'algor'", 98);
    }
        
    /*** Response matrix (row vector) ***/
    MatrixRect<int> *y;
    y = new MatrixRect<int>(1, *n, Y);
    if (!y) throw returnR(fname, "Out of memory (object 'y')", 99);

    /**** Design matrices         ****/
    /* * they are transposed!        */
    MatrixRect<double> *x, *v;
    x = new MatrixRect<double>(*p, *n, X);
    if (!x) throw returnR(fname, "Out of memory (object 'x')", 99);

    v = new MatrixRect<double>(*q, *n, V);
    if (!v) throw returnR(fname, "Out of memory (object 'v')", 99);

    /*** Matrices derived from the design matrices (needed to compute derivatives) ***/
    MatrixLT<double> *xx, *vv;
    MatrixRect<double> *xv;
    xx = new MatrixLT<double>[*n];
    if (!xx) throw returnR(fname, "Out of memory (object 'xx')", 99);
    cdP1 = X;
    for (k = 0; k < *n; k++){
      xx[k] = MatrixLT<double>(*p, cdP1);
      cdP1 += (*p);
    }

    vv = new MatrixLT<double>[*n];
    if (!vv) throw returnR(fname, "Out of memory (object 'vv')", 99);
    cdP1 = V;
    for (k = 0; k < *n; k++){
      vv[k] = MatrixLT<double>(*q, cdP1);
      cdP1 += (*q);
    }

    xv = new MatrixRect<double>[*n];
    if (!xv) throw returnR(fname, "Out of memory (object 'xv')", 99);
    cdP1 = X;
    cdP2 = V;
    for (k = 0; k < *n; k++){
      xv[k] = MatrixRect<double>(*p, *q, cdP1, cdP2);
      cdP1 += (*p);
      cdP2 += (*q);
    }
   
    /*** Regression coefficients (both beta and gamma) ***/
    MatrixRect<double> *theta_vec;
    theta_vec = new MatrixRect<double>(1, nBG, theta);
    if (!theta_vec) throw returnR(fname, "Out of memory (object 'theta_vec')", 99);

    /*** Offset term (zeros at this moment) ***/
    MatrixRect<double> *offset_vec;
    offset_vec = new MatrixRect<double>(1, *n);
    if (!offset_vec) throw returnR(fname, "Out of memory (object 'offset_vec')", 99);

    /*** Linear predictors ***/
    MatrixRect<double> *etaX, *etaV;
    etaX = new MatrixRect<double>(1, *n);
    if (!etaX) throw returnR(fname, "Out of memory (object 'etaX')", 99);
    if (*p) etaX->bA(theta_vec, x, Cq);
   
    etaV = new MatrixRect<double>(*C, *n);
    if (!etaV) throw returnR(fname, "Out of memory (object 'etaV')", 99);
    if (*q) etaV->bA2(theta_vec, v, 0);

    /*** Category probabilities              ***/
    /* (will be needed by Fisher scoring)      */
    int anyZero = 0;
    MatrixRect<double> *prob;
    prob = new MatrixRect<double>(*C+1, *n);
    if (!prob) throw returnR(fname, "Out of memory (object 'prob')", 99);
    util_cumlogit::prob_cumlogit(prob, offset_vec, etaX, etaV);

    /*** Working array, needed by Fisher scoring ***/
    double *h_dh;
    h_dh = (double*) calloc(5*(*C), sizeof(double));
    if (!h_dh) throw returnR(fname, "Out of memory (object 'h_dh')", 99);    

    /*** Derivative matrices ***/
    MatrixRect<double> *U_vec;
    U_vec = new MatrixRect<double>(1, nBG);
    if (!U_vec) throw returnR(fname, "Out of memory (object 'U_vec')", 99);

    MatrixRect<double> *NR_step;
    NR_step = new MatrixRect<double>(1, nBG);
    if (!NR_step) throw returnR(fname, "Out of memory (object 'NR_step')", 99);    

    MatrixLT<double> *I_mat;
    I_mat = new MatrixLT<double>(nBG);
    if (!I_mat) throw returnR(fname, "Out of memory (object 'I_mat')", 99);

    /*** Initial values of the log-likelihood, score and observed information matrix ***/
    eval_ll(ll, U_vec, I_mat, etaX, etaV, prob, h_dh, &anyZero, offset_vec, theta_vec, y, x, v, xx, vv, xv, 2);

    /*** Optimization (Newton-Raphson) ***/
    if (anyZero){
      *err = 6;
    }
    else{
      AK_Optim::optim_newton_raphson01(theta_vec, ll, U_vec, I_mat, NR_step, eval_ll,
                                       etaX, etaV, prob, h_dh, &anyZero, offset_vec, y, x, v, xx, vv, xv, 
                                       *niter, Fit_cumlogit::_AK_MAX_STEPHALF, *toler, *trace,
                                       Fit_cumlogit::_AK_MAX_PD_ATTEMPTS, Fit_cumlogit::_AK_EPS_PD_ATTEMPT, *err);

      if (!(*err)){

        /*** Copy values to return ***/
        theta_vec->mat2array(theta);
        U_vec->mat2array(U);
        I_mat->mat2array(I_obs, 1);

        /*** Expected information matrix and its inverse ***/
        ll_cumlogitFS(ll, U_vec, I_mat, etaX, etaV, prob, h_dh, &anyZero, offset_vec, theta_vec, y, x, v, xx, vv, xv, 2);

        if (anyZero){
          *err = 106;
        }
        else{
          itemp = I_mat->chol_inv(0, NULL);
          if (itemp < I_mat->nrow()){
            *err = 107;
          }
          else{
            I_mat->mat2array(I_exp, 1);
          }
        }
      }
    }

    /*** Clean ***/
    delete I_mat;
    delete NR_step;
    delete U_vec;
    free(h_dh);
    delete prob;
    delete etaV;
    delete etaX;
    delete offset_vec;
    delete theta_vec;
    delete [] xv;
    delete [] vv;
    delete [] xx;
    delete v;
    delete x;
    delete y;

    return;
  }     /* end of try */
  catch(returnR rr){
    *err = rr.errflag() + 1000;
    return;
  }
}  /* end of fit_cumlogit */

}  /* end of extern "C" */

}  /* end of the namespace Fit_cumlogit */
