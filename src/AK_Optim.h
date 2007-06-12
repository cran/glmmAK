/*** AK_Optim.h ***/

#ifndef _AK_OPTIM_H_
#define _AK_OPTIM_H_

#include <cmath>

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

#include "MatrixLT.h"
#include "MatrixRect.h"
#include "AK_BLAS_LAPACK.h"

namespace AK_Optim{

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
                       int &err);

void
optim_newton_raphson02(double* x,            double* gx,          double* dgx,  double* ddgx,  
                       const double* parmD,  const int* parmI,
	               void (*eval2)(const double*,  double*, double*, double*, const double*, const int*, const int&),
                       int* iter,            const int* maxiter,  const int* max_stephalf,  
                       const double* toler,  const double* zero,  int* err);

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
                       int &err);

void
solver_newton_raphson(double* x,            double* gx,          double* dgx,  const double* b,
                      const double* parmD,  const int* parmI,
	              void (*eval2)(const double*,  double*, double*, double*, const double*, const int*, const int&),
                      int* iter,            const int* maxiter,
                      const double* toler,  const double* zero,  int* err);

}

#endif
