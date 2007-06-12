/*** AK_Optim.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  13/08/2006
//
//     optim_newton_raphson01:  13/08/2006
//     optim_newton_raphson02:  15/01/2007
//     optim_newton_raphson03:  08/02/2007
//      solver_newton_raphson:  13/08/2006
//
// PURPOSE: Various optimization routines
//
/* ********************************************************************************* */

#include "AK_Optim.h"

namespace AK_Optim{

/************************************************************************************************/
/***                                                                                          ***/
/*** optim_newton_raphson01:  Optimization using the Newton-Raphson algorithm                 ***/
/***                          (written mainly for logit models)                               ***/
/***                                                                                          ***/
/***                                                                                          ***/
/************************************************************************************************/
//
// If niter = 0, the function only checks whether the information matrix is positive definite
//   and if so, computes its inverse.
//
// theta:    argument of the function we optimize
//           INPUT:   starting value
//           OUTPUT:  point of maximum
//
// ll:       value of the function we optimize
//           INPUT:   at starting point
//           OUTPUT:  at maximum
//
// d_ll:     value of the score vector
//           INPUT:  at starting point
//           OUTPUT: at maximum
//
// dd_ll:    value of the information matrix or its inverse
//           INPUT:  information matrix at starting point
//           OUTPUT: inverse of the information matrix at theta if (err=0, 1, 2)
//                   the information matrix at theta if (err=100, 101, 102)
// 
// NR_step:  matrix used to store the Newton-Raphson step
//
// eval_ll:  function used to compute the value of the function we maximize, its first and minus second derivatives
//           (e.g., function 'll_cumlogitNR' from ll_cumlogit.cpp)
// 
// etaX:
// etaV:
// prob:
// work:
// anyZero:
//
/// y:
// x:
// v:
// xx:
// vv:
// xv:             arguments needed to evaluate (*eval_ll), see, e.g. 'll_cumlogitNR' in ll_cumlogit.cpp
//
// niter:          INPUT:  maximal number of Newton-Raphson iterations
//                 OUTPUT: number of Newton-Raphson iterations performed
//
// max_stephalf:   maximum number of step-halving steps
// toler:          relative toleration to detect convergence
// trace:          if !=0 some information is printed during the iterations
//
// max_nAttempt:   maximum number of attempts to make the second derivative matrix positive definite if it is not
// eps_Attempt:    epsilon for MatrixLT::chol_solvePD function
//
// err:            error flag
//                 0 = everything O.K.
//                 100 = convergence reached but the final minus Hessian is not PD
//                 1 = maximum number of step-halving steps was performed without increase of the objective function
//                 101 = 1 and the final minus Hessian is not PD    
//                 2 = maximum number of iterations was used without convergence
//                 102 = 2 and the final minus Hessian is not PD
//                 3 = minus second derivative matrix not positive definite and not possible to make it PD
//                     by inflation of the diagonal
//                 4 = NaN appeared in initial ll(x), ll'(x) or ll''(x)
//                 5 = NaN appeared in ll(x), ll'(x) or ll''(x)
//                 6 = +-Inf appeared in the Hessian (can happen only with Fisher scoring)
//
//
void
optim_newton_raphson01(MatrixRect<double> *theta,
                       double *ll,
                       MatrixRect<double> *d_ll,
                       MatrixLT<double> *dd_ll,
                       MatrixRect<double> *NR_step,
                       void (*eval_ll)(double*, MatrixRect<double>*, MatrixLT<double>*, 
                                       MatrixRect<double>*, MatrixRect<double>*, MatrixRect<double>*, double*, int*,
                                       const MatrixRect<double>*, const MatrixRect<double>*, const MatrixRect<int>*, 
                                       const MatrixRect<double>*, const MatrixRect<double>*,
                                       const MatrixLT<double>*, const MatrixLT<double>*,
                                       const MatrixRect<double>*, const int&),
                       MatrixRect<double> *etaX,  
                       MatrixRect<double> *etaV,
                       MatrixRect<double> *prob,
                       double *work,
                       int *anyZero,
                       const MatrixRect<double> *offset,
                       const MatrixRect<int> *y,        
                       const MatrixRect<double> *x,     
                       const MatrixRect<double> *v,  
                       const MatrixLT<double> *xx,  
                       const MatrixLT<double> *vv,  
                       const MatrixRect<double> *xv,
                       int &niter,
                       const int &max_stephalf,
                       const double &toler,
                       const int &trace,
                       const int &max_nAttempt,
                       const double &eps_Attempt,
                       int &err)
{
  err = 0;

  int i, rank, Attempt, halfstep;
  double *aP, *bP, old_ll, relat_diff;
  const int maxiter = niter;

  /*** Check initials for NaN ***/
  if (!R_finite(*ll)){
    if (trace){
      Rprintf("\nInitial values lead to the log-likelihood of -Inf.\n");
    }
    err = 4;
    return;
  }
  // if (!R_finite(*ll) || !d_ll->is_finite() || !dd_ll->is_finite()){
  //   err = 4;
  //   return;
  // }
  
  /*** Iterate ***/
  for (niter = 0; niter < maxiter; niter++){

    if (trace){
      Rprintf("\n*** Iteration %d (log-likelihood = %g):\n", niter, *ll);
      Rprintf("   Theta: ");
      theta->print(0);
      Rprintf("   Score: ");
      d_ll->print(0);
      // Rprintf("   Minus Hessian:\n");
      // dd_ll->print(0);
    }
    
    /*** Newton-Raphson step, try to make the minus second derivative matrix positive definite if it is not ***/
    rank = dd_ll->chol_solvePD(Attempt, d_ll->a(), NR_step->a(), max_nAttempt, eps_Attempt);
    if (trace){
      Rprintf("   Number of PD attempts = %d\n", Attempt);
      Rprintf("   NR-step: ");
      NR_step->print(0);
    }
    if (rank < dd_ll->nrow()){
      err = 3;
      return;
    }

    /*** Compute new value of theta ***/
    theta->add_A(NR_step);    

    /*** Update derivatives ***/
    old_ll = *ll;
    eval_ll(ll, d_ll, dd_ll, etaX, etaV, prob, work, anyZero, offset, theta, y, x, v, xx, vv, xv, 2);
    if (*anyZero){
      err = 6;
      return;
    }
    
    /*** Check convergence ***/
    relat_diff = R_finite(*ll) ? fabs(1 - old_ll/(*ll)) : R_PosInf;
    if (relat_diff <= toler){
      if (trace){
        Rprintf("\nConvergence reached with log-likelihood = %g\n", *ll);
        Rprintf("   Theta: ");
        theta->print(0);
        Rprintf("   Score: ");
        d_ll->print(0);
        // Rprintf("   Minus Hessian:\n");
        // dd_ll->print(0);
      }
      break;
    }

    /*** If not yet convergence, check whether the objective function increases ***/
    /***   if no increase perform step-halving                                  ***/
    if (!R_finite(*ll) || *ll < old_ll){      
      if (trace){
        Rprintf("   Step-halfing (0) -> log-likelihood = %g\n", *ll);
      }
      for (halfstep = 0; halfstep < max_stephalf; halfstep++){
        aP = theta->a();
        bP = NR_step->a();
        for (i = 0; i < theta->nrow(); i++){
          *bP *= 0.5;
          *aP -= (*bP);
          aP++;
          bP++;
        }
        eval_ll(ll, d_ll, dd_ll, etaX, etaV, prob, work, anyZero, offset, theta, y, x, v, xx, vv, xv, 0);
        if (trace){
          Rprintf("   Step-halfing (%d) -> log-likelihood = %g\n", halfstep+1, *ll);
        }
        if (*ll >= old_ll){
          eval_ll(ll, d_ll, dd_ll, etaX, etaV, prob, work, anyZero, offset, theta, y, x, v, xx, vv, xv, 2);
          if (*anyZero){
            err = 6;
            return;
          }
          break;
        }
      }
      if (halfstep == max_stephalf){
        eval_ll(ll, d_ll, dd_ll, etaX, etaV, prob, work, anyZero, offset, theta, y, x, v, xx, vv, xv, 2);
        if (*anyZero){
          err = 6;
          return;
        }
        err = 1;
        break;
      }
    }

    /*** Check values of ll and its derivatives for NaN ***/
    // if (!R_finite(*ll) || !d_ll->is_finite() || !dd_ll->is_finite()){
    //   err = 5;
    //   return;
    // }

  }  /*** end of iterations ***/
  
  /*** Check values of ll and its derivatives for NaN ***/
  // if (!R_finite(*ll) || !d_ll->is_finite() || !dd_ll->is_finite()){
  //   err = 5;
  //   return;
  // }

  /*** Check niter ***/
  if (maxiter && niter == maxiter) err = 2;
  else                             niter++;     /* to get the number of iterations really performed */

  /*** Invert the minus second derivative matrix and check whether it is PD ***/
  rank = dd_ll->chol_inv(0, NULL);
  if (rank < dd_ll->nrow()) err += 100;
  return;    
}


/************************************************************************************************/
/***                                                                                          ***/
/*** optim_newton_raphson02:  Optimization using the Newton-Raphson algorithm                 ***/
/***                More or less the same as 'newton_raphson' function from                   ***/
/***                newton_raphson.cpp of bayesSurv package                                   ***/
/***                                                                                          ***/
/***                                                                                          ***/
/************************************************************************************************/
//
// x ............. on INPUT: starting point for the Newton-Raphson
//                 on OUTPUT: point where the maximum is reached
// gx ............ on INPUT : value of g at the starting point
//                 on OUTPUT: value of g at the maximum
// dgx............ on INPUT : value of g' at the starting point
//                 on OUTPUT: value of g' at the maximum
// ddgx........... on INPUT : value of -g'' at the starting point
//                 on OUTPUT: value of -g'' at the maximum
// parmD ......... additional parameters to evaluate g
// parmI ......... additional parameters to evaluate g
// eval2 ......... routine to compute the function we want to maximize and its first and second derivatives
//                 const double* ..... point x where to evalulate the function g
//                 double* ........... g(x)
//                 double* ........... g'(x)
//                 double* ........... -g''(x)
//                 const double* ..... additional parameters to evaluate g
//                 const int* ........ additional parameters to evaluate g
//                 const int& ........ integer indicating what should be computed
//                                   0 = all g(x), g'(x), g''(x)
//                                   1 = only g(x)
//                                   2 = only g'(x) and g''(x)
// iter .......... on OUTPUT: number of iterations needed
// maxiter ....... maximum number of iterations
// max_stephalf... maximum number of step-halving steps
// toler ......... relative toleration to detect convergence
// zero .......... small number to detect singular second derivative
// err ........... error flag
//                 1 = maximum number of step-halving steps was performed without increase of the objective function
//                 2 = maximum number of iterations was used without convergence
//                 3 = bad initials
//                 4 = NaN appeared in g'(x) or g''(x)
//
void
optim_newton_raphson02(double* x,            double* gx,          double* dgx,  double* ddgx,  
                       const double* parmD,  const int* parmI,
	               void (*eval2)(const double*,  double*, double*, double*, const double*, const int*, const int&),
                       int* iter,            const int* maxiter,  const int* max_stephalf,  
                       const double* toler,  const double* zero,  int* err)
{
  *err = 0;
  static double newgx, newx, NRstep, relat_diff;
  static int halfstep;

  if (!R_finite(*gx) || !R_finite(*dgx) || !R_finite(*ddgx)){
    *err = 3;
    return;
  }

  for (*iter = 0; *iter < *maxiter; (*iter)++){
    if (fabs(*ddgx) <= *zero) *ddgx = (*zero);
    NRstep = (*dgx)/(*ddgx);
    newx = (*x) + NRstep;
    for (halfstep = 0; halfstep < *max_stephalf; halfstep++){
      eval2(&newx, &newgx, dgx, ddgx, parmD, parmI, 1);
      relat_diff = fabs(1 - ((*gx)/newgx));
      if (newgx >= *gx || relat_diff <= *toler) break;
      newx = 0.5*((*x) + newx);
    }
    if (halfstep == *max_stephalf){
      *err = 1;
      return;
    }
    *x = newx;
    *gx = newgx;
    eval2(x, gx, dgx, ddgx, parmD, parmI, 2);
    if (!R_finite(*dgx) || !R_finite(*ddgx)){
      *err = 4;
      return;
    }    
    if (relat_diff <= *toler) break;
  }
  
  if (*iter == *maxiter) *err = 2;
  return;
}


/************************************************************************************************/
/***                                                                                          ***/
/*** optim_newton_raphson03:  Optimization using the Newton-Raphson algorithm                 ***/
/***                          (written mainly for poisson models)                             ***/
/***                                                                                          ***/
/***                                                                                          ***/
/************************************************************************************************/
//
// If niter = 0, the function only checks whether the information matrix is positive definite
//   and if so, computes its inverse.
//
// theta:    argument of the function we optimize
//           INPUT:   starting value
//           OUTPUT:  point of maximum
//
// ll:       value of the function we optimize
//           INPUT:   at starting point
//           OUTPUT:  at maximum
//
// d_ll:     value of the score vector
//           INPUT:  at starting point
//           OUTPUT: at maximum
//
// dd_ll[LT(p)]:    value of the information matrix or its inverse
//                  INPUT:  information matrix at starting point
//                  OUTPUT: inverse of the information matrix at theta if (err=0, 1, 2)
//                          the information matrix at theta if (err=100, 101, 102)
//
// dd_ll_temp[LT(p)]:  working space
// 
// NR_step:  matrix used to store the Newton-Raphson step
//
// eval_ll:  function used to compute the value of the function we maximize, its first and minus second derivatives
//           (e.g., function 'll_poisson' from ll_poisson.cpp)
// 
// eta:
// mu:
// offset:
//
// y:
// log_y_factor:
// x:
// xx:            
// nObs:
// p:              arguments needed to evaluate (*eval_ll), see, e.g. 'll_poisson' in ll_poisson.cpp
//              
// niter:          INPUT:  maximal number of Newton-Raphson iterations
//                 OUTPUT: number of Newton-Raphson iterations performed
//
// max_stephalf:   maximum number of step-halving steps
// toler:          relative toleration to detect convergence
// trace:          if !=0 some information is printed during the iterations
//
// max_nAttempt:   maximum number of attempts to make the second derivative matrix positive definite if it is not
// eps_Attempt:    epsilon for MatrixLT::chol_solvePD function
//
// err:            error flag
//                 0 = everything O.K.
//                 100 = convergence reached but the final minus Hessian is not PD
//                 1 = maximum number of step-halving steps was performed without increase of the objective function
//                 101 = 1 and the final minus Hessian is not PD    
//                 2 = maximum number of iterations was used without convergence
//                 102 = 2 and the final minus Hessian is not PD
//                 3 = minus second derivative matrix not positive definite and not possible to make it PD
//                     by inflation of the diagonal
//                 4 = NaN appeared in initial ll(x), ll'(x) or ll''(x)
//                 5 = NaN appeared in ll(x), ll'(x) or ll''(x)
//
//
void
optim_newton_raphson03(double *theta,
                       double *ll,
                       double *d_ll,
                       double *dd_ll,
                       double *dd_ll_temp,
                       double *NR_step,
                       void (*eval_ll)(double*, double*, double*, double*, double*, const double*, const double*, 
                                       const int*, const double*, const double*, const MatrixLT<double>*, const int&, const int&, const int&),
                       double *eta,  
                       double *mu,
                       const double *offset,
                       const int *y,        
                       const double *log_y_factor,
                       const double *x,     
                       const MatrixLT<double> *xx,  
                       const int &nObs,
                       const int &p,
                       int &niter,
                       const int &max_stephalf,
                       const double &toler,
                       const int &trace,
                       const int &max_nAttempt,
                       const double &eps_Attempt,
                       int &err)
{
  err = 0;

  const int _ONE_ = 1;

  static int i, info, Attempt, halfstep;
  static double *aP, *bP, old_ll, relat_diff;
  const int maxiter = niter;

  /*** Check initials for NaN ***/
  if (!R_finite(*ll)){
    if (trace){
      Rprintf("\nInitial values lead to the log-likelihood of -Inf.\n");
    }
    err = 4;
    return;
  }
  // if (!R_finite(*ll) || !d_ll->is_finite() || !dd_ll->is_finite()){
  //   err = 4;
  //   return;
  // }
  
  /*** Iterate ***/
  for (niter = 0; niter < maxiter; niter++){

    if (trace){
      Rprintf("\n*** Iteration %d (log-likelihood = %g):\n", niter, *ll);
      Rprintf("   Theta: ");
      AK_BLAS_LAPACK::printArray(theta, p);
      Rprintf("   Score: ");
      AK_BLAS_LAPACK::printArray(d_ll, p);
      // Rprintf("   Minus Hessian:\n");
      // AK_BLAS_LAPACK::printLT(dd_ll, p);
    }
    
    /*** Decomposition of the Hessian, try to make the minus second derivative matrix positive definite if it is not ***/
    AK_BLAS_LAPACK::chol_dpptrfPD(dd_ll, dd_ll_temp, &p, &Attempt, &max_nAttempt, &eps_Attempt, &info);
    if (info){
      err = 3;
      return;
    }

    /*** Compute the Newton-Raphson step ***/
    AK_BLAS_LAPACK::copyArray(NR_step, d_ll, p);
    AK_BLAS_LAPACK::chol_solve_system(NR_step, dd_ll, &p, &_ONE_);
    if (trace){
      Rprintf("   Number of PD attempts = %d\n", Attempt);
      Rprintf("   NR-step: ");
      AK_BLAS_LAPACK::printArray(NR_step, p);
    }

    /*** Compute new value of theta ***/
    AK_BLAS_LAPACK::a_aPlusb(theta, NR_step, p);


    /*** Update derivatives ***/
    old_ll = *ll;
    eval_ll(ll, d_ll, dd_ll, eta, mu, offset, theta, y, log_y_factor, x, xx, nObs, p, 2);
    
    /*** Check convergence ***/
    relat_diff = R_finite(*ll) ? fabs(1 - old_ll/(*ll)) : R_PosInf;
    if (relat_diff <= toler){
      if (trace){
        Rprintf("\nConvergence reached with log-likelihood = %g\n", *ll);
        Rprintf("   Theta: ");
        AK_BLAS_LAPACK::printArray(theta, p);
        Rprintf("   Score: ");
        AK_BLAS_LAPACK::printArray(d_ll, p);
        // Rprintf("   Minus Hessian:\n");
        // AK_BLAS_LAPACK::printLT(dd_ll, p);
      }
      break;
    }

    /*** If not yet convergence, check whether the objective function increases ***/
    /***   if no increase perform step-halving                                  ***/
    if (!R_finite(*ll) || *ll < old_ll){      
      if (trace){
        Rprintf("   Step-halfing (0) -> log-likelihood = %g\n", *ll);
      }
      for (halfstep = 0; halfstep < max_stephalf; halfstep++){
        aP = theta;
        bP = NR_step;
        for (i = 0; i < p; i++){
          *bP *= 0.5;
          *aP -= (*bP);
          aP++;
          bP++;
        }
        eval_ll(ll, d_ll, dd_ll, eta, mu, offset, theta, y, log_y_factor, x, xx, nObs, p, 0);
        if (trace){
          Rprintf("   Step-halfing (%d) -> log-likelihood = %g\n", halfstep+1, *ll);
        }
        if (*ll >= old_ll){
          eval_ll(ll, d_ll, dd_ll, eta, mu, offset, theta, y, log_y_factor, x, xx, nObs, p, 2);
          break;
        }
      }
      if (halfstep == max_stephalf){
        eval_ll(ll, d_ll, dd_ll, eta, mu, offset, theta, y, log_y_factor, x, xx, nObs, p, 2);
        err = 1;
        break;
      }
    }

    /*** Check values of ll and its derivatives for NaN ***/
    // if (!R_finite(*ll) || !d_ll->is_finite() || !dd_ll->is_finite()){
    //   err = 5;
    //   return;
    // }

  }  /*** end of iterations ***/
  
  /*** Check values of ll and its derivatives for NaN ***/
  // if (!R_finite(*ll) || !d_ll->is_finite() || !dd_ll->is_finite()){
  //   err = 5;
  //   return;
  // }

  /*** Check niter ***/
  if (maxiter && niter == maxiter) err = 2;
  else                             niter++;     /* to get the number of iterations really performed */

  /*** Invert the minus second derivative matrix and check whether it is PD ***/
  AK_BLAS_LAPACK::chol_dpptrf(dd_ll, &p, &info);
  if (info){
    err += 100;
    return;
  }
  AK_BLAS_LAPACK::chol_dpptri(dd_ll, &p, &info);
  if (info){
    err += 100;
    return;
  }
  return;    
}


/**************************************************************************************************************/
/*****                                                                                                    *****/
/***** solver_newton_raphson:  Find a solution of an equation g(x) = b using a Newton-Raphson algorithm   *****/
/*****                                                                                                    *****/
/**************************************************************************************************************/
//
// Taken from newton_raphson.cpp of bayesSurv package
// 
// ASSUMPTIONS: * starting point is in correct region (i.e. only one solution is found)
//              * there exists a solution
//         typical usage: g(x) is concave and unimodal with no limits imposed on the range of the function 'g'
//
// x ............. on INPUT: starting point for the Newton-Raphson
//                 on OUTPUT: point where the equation is solved
// gx ............ on INPUT : value of g at the starting point
//                 on OUTPUT: value of g at the solution (this should be close to b)
// dgx............ on INPUT : value of g' at the starting point
//                 on OUTPUT: value of g' at the solution
// b ............. right hand side
// parmD ......... additional parameters to evaluate g
// parmI ......... additional parameters to evaluate g
// eval2 ......... routine to compute the function g(x)and its first derivative
//                 const double* ..... point x where to evalulate the function g
//                 double* ........... g(x)
//                 double* ........... g'(x)
//                 double* ........... -g''(x) (nowhere needed by this routine)
//                 const double* ..... additional parameters to evaluate g
//                 const int* ........ additional parameters to evaluate g
//                 const int& ........ integer indicating what should be computed
//                                   0 = all g(x), g'(x), g''(x)
//                                   1 = only g(x)
//                                   2 = only g'(x) and g''(x)
//                                   3 = only g(x) and g'(x)
// iter .......... on OUTPUT: number of iterations needed
// maxiter ....... maximum number of iterations
// toler ......... relative toleration to detect convergence
// zero .......... small number to detect singular first derivative
// err ........... error flag
//                 (1) = maximum number of step-halving steps was performed without increase of the objective function (not applicable here)
//                 2   = maximum number of iterations was used without convergence
//                 3   = bad initials
//                 4   = NaN appeared in g(x) or g'(x)
//
void
solver_newton_raphson(double* x,            double* gx,          double* dgx,  const double* b,
                      const double* parmD,  const int* parmI,
	              void (*eval2)(const double*,  double*, double*, double*, const double*, const int*, const int&),
                      int* iter,            const int* maxiter,
                      const double* toler,  const double* zero,  int* err)
{
  *err = 0;
  static double ddgx, NRstep, _diff;

  if (!R_finite(*gx) || !R_finite(*dgx) || !R_finite(*b)){
    REprintf("x=%g,  gx=%g,  dgx=%g,  b=%g\n", *x, *gx, *dgx, *b);
    *err = 3;
    return;
  }

  _diff = *b - (*gx);
  for (*iter = 0; *iter < *maxiter; (*iter)++){
    if (fabs(*dgx) <= *zero) *dgx = (*zero);
    NRstep = _diff/(*dgx);
    (*x) += NRstep;
    eval2(x, gx, dgx, &ddgx, parmD, parmI, 3);
    if (!R_finite(*gx) || !R_finite(*dgx)){
      REprintf("x=%g,  gx=%g,  dgx=%g,  b=%g\n", *x, *gx, *dgx, *b);
      *err = 4;
      return;
    }    
    _diff = *b - (*gx);
    if (fabs(_diff/(*b)) <= *toler) break;
  }
  
  if (*iter == *maxiter) *err = 2;
  return;
}

}  /*** end of the namespace AK_Optim ***/
