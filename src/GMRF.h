/*** GMRF.h ***/

#ifndef _GAUSSIAN_MARKOV_RANDOM_FIELD_H_
#define _GAUSSIAN_MARKOV_RANDOM_FIELD_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_BLAS_LAPACK.h"

namespace GMRF {

const double _AK_TOL_CHOL2 = 1e-300;        /*** Tolerance for the Cholesky decomposition used here ***/
const double _AK_SQRT_TOL_CHOL2 = 1e-150;

extern "C"{
    void
  diff(double* Da,  const int *order,  const int *na);
}

extern "C"{
  void
  tdiff(double *Qa,  const double *Da,  const int *diffOper,  const int *order,  const int *na);
}

extern "C"{
  void
  diff_operator(int *diffOper,  const int *order);
}

extern "C"{
  void
  Q_matrix(double *Q,  const int *order,  const int *na);
}

void
log_density_Ax_x(double *VALUE, const double *A,  const int *nx,  const int *nc,  double *work);

void
rGMRF(double *x,                    double *log_dens,
      const double *mu,             const double *Li,  const double *log_dets,     
      const double *A,              const double *e,   const double *U,  
      const double *log_dens_Ax_x, 
      const int *nx,                const int *nc,     
      const int *mu_nonZERO,        const int *e_nonZERO,    
      double *work);

void
rGMRF_inputArgs(double *log_dets,       double *U,        
                const double *mu,       const double *Li,  const double *A,  const double *e,
                const int *nx,          const int *nc,  
                const int *mu_nonZERO,  const int *e_nonZERO,
                double *work);

void
dGMRF(double *VALUE,       
      const double *x,              const int *unlog,
      const double *mu,             const double *Li,   const double *log_dets,   
      const double *log_dens_Ax_x,  
      const int *nx,                const int *nc,              
      const int *mu_nonZERO, 
      double *work);

void
dGMRF_inputArgs(double *log_dets,       double *LW,        double *V,         
                const double *mu,       const double *Li,  const double *A,  const double *e,
                const int *nx,          const int *nc,  
                const int *mu_nonZERO,  const int *e_nonZERO,
                double *work);

void
dGMRF2(double *VALUE,       
       const double *x,              const int *unlog,
       const double *mu,             const double *Li,           const double *log_dets,   
       const double *mustar,         const double *LiSigmastar,
       const int *nx,                const int *nc,              
       const int *mu_nonZERO, 
       double *work);

void
dGMRF2_inputArgs(double *log_dets,
                 double *mustar,         double *LiSigmastar,
                 const double *mu,       const double *Li,        const double *A,   const double *e,
                 const int *nx,          const int *nc,  
                 const int *mu_nonZERO,  const int *e_nonZERO,    
                 double *work);

void
dscale_norm_const(const double *F, double *RES);

void
dscale(const double *x,  double *gx,  double *dgx,  double *ddgx,  const double *parD,  const int *parI,  const int &what);

void
rscale(double *x,  const double *parD);


extern "C"{
  void
  rGMRFR(double *x,              double *log_dens,  
         const double *mu,       const double *Li,      const double *A,     const double *e,
         const int *nx,          const int *nc,         const int *nrandom,  
         const int *mu_nonZERO,  const int *e_nonZERO);
}

extern "C"{
  void
  dGMRFR(double *VALUE,  
         const double *x,        const int *unlog,
         const double *mu,       const double *Li,     const double *A,     const double *e,  
         const int *nx,          const int *nc,        const int *npoints,
         const int *mu_nonZERO,  const int *e_nonZERO);
}

extern "C"{
  void
  dGMRF2R(double *VALUE,  
          double *mustar,         double *LiSigmastar,
          const double *x,        const int *unlog,
          const double *mu,       const double *Li,     const double *A,     const double *e,  
          const int *nx,          const int *nc,        const int *npoints,
          const int *mu_nonZERO,  const int *e_nonZERO);
}

extern "C"{
  void
  rscaleR(double *x,  const int *n,  const double *F);
}

}  /*** end of namespace GMRF ***/

#endif
