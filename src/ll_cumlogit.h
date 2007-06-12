/*** ll_cumlogit.h ***/

#ifndef _LL_CUMLOGIT_H_
#define _LL_CUMLOGIT_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

#include "MatrixLT.h"
#include "MatrixRect.h"

#include "util_cumlogit.h"

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
              const MatrixRect<double> *vx,
              const int &order=2);

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
              const MatrixRect<double> *vx,
              const int &order=2);

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
               const MatrixRect<double> *vx,
               const int &nObs,
               const int &C,
               const int &p,
               const int &q,
               const int &update_etas=1,               /* IGNORED HERE  */
               const int &order=2);

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
               const MatrixRect<double> *vx,
               const int &nObs,
               const int &C,                         /* number of response categories minus 1 */
               const int &p,                         /* number of beta regressors             */
               const int &q,                         /* number of gamma regressors            */
               const int &update_etas=1,             /* 0/1 indicating whether linear predictors are to be updated   */
               const int &order=2);

void
scoreNR_cumlogit(double *U,            double *U_X,           
                 double &dh0,          double &dh1,           double &dh0_h0_h1,   double &dh1_h0_h1,
                 const int *y,         const double **x,      const double **v,
                 const int &p,         const int &q,          const int &C,        const int &Cq,
                 const double &h0,     const double &h1,      const double &h0_h1);

void
scoreFS_cumlogit(double *U,            double *U_X,
                 const int *y,         const double **x,      const double **v,
                 const int &p,         const int &q,          const int &C,      const int &Cq,
                 const double &h0_h1,  const double *dh);

void
infoMatNR_cumlogit(double *I,             double *I_XX,      const int *diagI,
	           const int *y,          const double *xx,  const double *vv,  const double *xv,
                   const int &p,          const int &q,      const int &C,      const int &Cq,           const int &lvv,
                   const double& h0_h1,   const double &dh0, const double &dh1, const double &dh0_h0_h1, const double &dh1_h0_h1);

void
infoMatFS_cumlogit(double *I, 
                   double *num_xv,        double *num_vv,    double *num_vv2,
	           const int *y,          const double *xx,  const double *vv,  const double *xv,
                   const int &p,          const int &q,      const int &C,      const int &Cq,     const int &lvv,
                   const double *prob,    const double *dh);

void
infoMatSDE_SDS_cumlogit(double *I,  const int &p,  const int &q,  const int &C);

#endif

