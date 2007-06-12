/*** GMRF.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  07/11/2006
//
//    PURPOSE:  Functions that implement (intrinsic) Gaussian random fields
//              and some useful utilities
//
//                          diff:  08/11/2006
//                         tdiff:  09/11/2006
//                 diff_operator:  08/11/2006
//                      Q_matrix:  18/11/2006
//
//              log_density_Ax_x:  17/11/2006
//                         rGMRF:  13/11/2006
//               rGMRF_inputArgs:  13/11/2006
//                         dGMRF:  17/11/2006
//               dGMRF_inputArgs:  17/11/2006
//                        dGMRF2:  16/11/2006
//              dGMRF2_inputArgs:  16/11/2006
//
//             dscale_norm_const:  07/11/2006
//                        dscale:  07/11/2006
//                        rscale:  07/11/2006
//
//                        rGMRFR:  13/11/2006
//                        dGMRFR:  17/11/2006
//                       dGMRF2R:  16/11/2006
//                       rscaleR:  07/11/2006
//                         
/* ********************************************************* */
 
#include "GMRF.h"

namespace GMRF {

/*** ============================================================================= ***/
/*** Difference operator applied to a vector                                       ***/
/***  = D * a, where D is (na-order) x na difference operator matrix               ***/
/***           and a is vector of length na                                        ***/
/***                                                                               ***/
/*** ============================================================================= ***/
//
// order:  order of the differences (at least 1)
// na:    length of the whole vector 'a'
// Da:     INPUT: vector for which differences should be computed
//        OUTPUT: computed differences (at places 0, ..., na-1-order)
extern "C"{
  void
  diff(double *Da,  const int *order,  const int *na)
  {
    static int j, d;
    static double *DaP;

    if (*order < 0 || *order > *na-1){
      REprintf("diff:  order=%d,  na=%d\n", *order, *na);
      throw returnR("Error in GMRF.cpp: diff, order must be >= 0 & <= na-1", 1);
    }

    for (d = 1; d <= *order; d++){
      DaP = Da;
      for (j = 0; j < *na-d; j++){
        *DaP = *(DaP+1) - *DaP;
        DaP++;
      }
    }

    return;
  }  
}  /*** end of extern "C" ***/


/*** ================================================================================== ***/
/*** Transposed difference operator applied to a vector                                 ***/
/***  = t(D) * b, where t(D) is na x (na-order) transposed difference operator matrix   ***/
/***              and b is a vector of length (na-order)                                ***/
/***                                                                                    ***/
/*** TYPICAL USAGE:        b = D * a,                                                   ***/
/***                t(D) * b = t(D) * D * a = Q * a                                     ***/
/***                                                                                    ***/
/*** ================================================================================== ***/
//
// Qa[na]:           computed t(D) * b
// Da[na-order]:     input vector b
// diffOper[order]:  vector defining the difference operator (see below)
// order:
// na:
//
extern "C"{
  void
  tdiff(double *Qa,  const double *Da,  const int *diffOper,  const int *order,  const int *na)
  {
    //Rprintf("\ntdiff: Da (order=%d, na=%d): ", *order, *na);
    //AK_BLAS_LAPACK::printArray(Da, *na);
    //Rprintf("tdiff: diffOper:  ");
    //AK_BLAS_LAPACK::printIArray(diffOper, *order+1);

    static int i, j, end, nDa_1;
    static double *QaP;
    static const int *diffOperEndP, *diffOperP;
    static const double *DaStartP, *DaP;

    if (*order < 0 || *order > *na-1){
      REprintf("tdiff:  order=%d,  na=%d\n", *order, *na);
      throw returnR("Error in GMRF.cpp: tdiff, order must be >= 0 & <= na-1", 1);
    }

    diffOperEndP = diffOper + (*order);
    nDa_1 = *na - (*order) - 1;                  /** length of Da minus 1 **/

    /** Qa[0,...,order-1] **/
    QaP      = Qa;
    DaStartP = Da; 
    for (i = 0; i < *order; i++){
      diffOperP = diffOper + i;            /** t(D)[i, 0]   **/
      DaP = DaStartP;                      /** Da[0]        **/
      *QaP = 0.0;                          /** Qa[i] = 0.0  **/
      end = (i <= nDa_1) ? i : nDa_1;
      for (j = 0; j <= end; j++){
        *QaP += (*diffOperP) * (*DaP);
        diffOperP--;                     /** t(D)[i, j+1] **/
        DaP++;                           /** Da[j+1]      **/
      }
      QaP++;
    }

    /** Qa[order,...,na-1] **/
    for (i = *order; i < *na; i++){
      diffOperP = diffOperEndP;
      DaP = DaStartP;
      *QaP = 0.0;                          /** Qa[i] = 0.0  **/
      end = (i <= nDa_1) ? i : nDa_1;
      for (j = i-(*order); j <= end; j++){
        *QaP += (*diffOperP) * (*DaP);
        diffOperP--;                     /** t(D)[i, j+1] **/
        DaP++;                           /** Da[j+1]      **/
      }
      QaP++;
      DaStartP++;
    }
 
    return;
  }
}


/*** ============================================================================= ***/
/*** Vector defining the difference operator                                       ***/
/***                                                                               ***/
/*** ============================================================================= ***/
//
// diffOper[order+1]:  vector defining the difference operator, e.g.,
//                        order = 0:   1
//                        order = 1:   -1, 1
//                        order = 2:   1, -2, 1
//                        order = 3:   -1, 3, -3, 1
//                         etc.
extern "C"{
  void
  diff_operator(int *diffOper,  const int *order)
  {   
    int j;

    if (*order < 0){
      REprintf("diff_operator:  order=%d\n", *order);
      throw returnR("Error in GMRF.cpp: diff_operator, order must be >= 0", 1);
    }

    int *work = (int*) calloc(*order+1, sizeof(int));    
    if (!work) throw returnR("Error in GMRF.cpp: diff_operator, out of memory", 99);
    
    *diffOper = 1;
    for (int d=1; d <= *order; d++){
      work[0] = 0;
      for (j = 0; j < d; j++){
        work[j+1] = diffOper[j];
        diffOper[j] *= (-1);
      }
      for (j = 0; j <= d; j++){      
        diffOper[j] += work[j];
      }
    }

    free(work);

    return;
  }
}


/*** ============================================================================= ***/
/*** Matrix Q = t(D) * D                                                           ***/
/***     (in packed form as lower triangle)                                        ***/
/***                                                                               ***/
/*** ============================================================================= ***/
//
// Q[LT(na)]: Computed matrix Q = t(D)*D (its lower triangle)
//
extern "C"{
  void
  Q_matrix(double *Q,  const int *order,  const int *na)
  {
    int j, i;
    double *QP, *tDP;
    const int *diffOperP;
    const double * QtempP;

    if (*order < 0 || *order > *na-1){
      REprintf("Q_matrix:  order=%d,  na=%d\n", *order, *na);
      throw returnR("Error in GMRF.cpp: Q_matrix, order must be >= 0 & <= na-1", 1);
    }

    if (*order == 0){    /** Q is identity matrix **/
      QP = Q;
      for (j = 0; j < *na; j++){
        *QP = 1;
        QP++;
        for (i = j+1; i < *na; i++){
          *QP = 0;
          QP++;
        }
      }
      return;
    }

    /***** order > 0 *****/
    /***** ========= *****/
    int na_order = *na - (*order);

    /** Difference operator **/
    int *diffOper = (int*) calloc(*order+1, sizeof(int));    
    if (!diffOper) throw returnR("Error in GMRF.cpp: Q_matrix, out of memory", 99);
    diff_operator(diffOper, order);

    /** Matrix D' **/
    double *tD = (double*) calloc((*na) * na_order, sizeof(double));
    if (!tD) throw returnR("Error in GMRF.cpp: Q_matrix, out of memory", 99);
    tDP =tD; 
    for (j = 0; j < na_order; j++){     
      for (i = 0; i < j; i++){
        *tDP = 0;
        tDP++;
      }
      diffOperP = diffOper;
      for (i = j; i <= j + *order; i++){
        *tDP = *diffOperP;
        tDP++;
        diffOperP++;
      }
      for (i = j + *order + 1; i < *na; i++){
        *tDP = 0;
        tDP++;
      }
    }
    free(diffOper);

    /** Matrix Q = D'*D **/
    double *Qtemp = (double*) calloc((*na) * (*na), sizeof(double));
    if (!Qtemp) throw returnR("Error in GMRF.cpp: Q_matrix, out of memory", 99);
    AK_BLAS_LAPACK::C_AtB(Qtemp, tD, tD, na, &na_order, na);
    QP     = Q;
    QtempP = Qtemp;
    for (j = 0; j < *na; j++){
      for (i = 0; i < j; i++){
        QtempP++;
      }
      for (i = j; i < *na; i++){
        *QP = *QtempP;
        QP++;
        QtempP++;
      }
    }

    free(tD);
    free(Qtemp);
    return;
  }
}


/*** ======================================================================================================= ***/
/*** Logarithm of pi(Ax|x) = -0.5*log|A*t(A)|                                                                ***/
/***                                                                                                         ***/
/*** ======================================================================================================= ***/
//
//  A[nc x nx]
//  work[nc*nc]:   Working array
//                 * length = nc*nc
//                 * used to store A * t(A)
//
void
log_density_Ax_x(double *VALUE, const double *A,  const int *nx,  const int *nc,  double *work)
{
  int info;

  if (*nc == 0){
    *VALUE = 0;
    return;
  }

  /** Compute A*t(A) **/
  AK_BLAS_LAPACK::C_AtB(work, A, A, nc, nx, nc);

  /** Cholesky decomposition of A*t(A) **/
  AK_BLAS_LAPACK::chol_dpotrf(work, nc, &info);
  if (info) throw returnR("Error in GMRF.cpp: log_density_Ax_x(). Constraints in A are linearly dependent.", 1);

  /** Log-density = -sum(log(diagonal elements of the decomposition)) **/
  const double *LP = work;
  *VALUE = 0;
  for (int i = 0; i < *nc; i++){
    *VALUE += log_AK(*LP);
    LP += *nc + 1;
  }
  *VALUE *= -1;

  return;
}


/*** ======================================================================================================= ***/
/*** Sampling from GMRF, possibly under a hard constraint                                                    ***/
/***   Algorithm 2.6 (p. 38), Rue and Held (2005). Gaussian Markov Random Fields: Theiry and Applications    ***/
/***      (conditioning by kriging)                                                                          ***/
/***                                                                                                         ***/
/*** Return also the value of the log-density evaluated in the sampled point                                 ***/
/***                                                                                                         ***/
/***   x ~ N(mu, Q^{-1})                                                                                     ***/
/***       under A*x = e                                                                                     ***/
/***                                                                                                         ***/
/***     Q = Li*t(Li)                                                                                        ***/
/***                                                                                                         ***/
/*** ======================================================================================================= ***/
//
// OUTPUT ARGUMENTS
// ================
// x[nx]:               Sampled value
// log_dens[1]:         Value of the log-density evaluated at sampled point
//
// PARAMETERS OF THE GMRF
// ======================
// mu[nx]:                  Mean of the unconstrained GMRF 
//                            (not needed when mu_nonZERO = 0 in which case it can be set to NULL)
// Li[LT(nx)]:              Cholesky decomposition of the unconstrained precision matrix Q,
//                            array of length nx*(nx+1)/2
// log_dets[4]:             log_dets[0] = log(|Q|^{1/2})  = sum(log(Li[j,j])),    Q = Li*t(Li)
//                          log_dets[1] = -(nx - nc)*log(sqrt(2pi))
//                          log_dets[2] = log(|W|^{-1/2}) = -sum(log(LW[j,j])),   W = LW*t(LW)
//                          log_dets[3] = -0.5*t(A*mu - e)*W^{-1}*(A*mu - e)
//
// A[nc x nx]:              Matrix defining the constraints (stored in column major order)
// e[nc]:                   Right-hand side of the constraint (not needed when e_nonZERO = 0 in which case it can be set to NULL)
// U[nc x nx]:              Matrix U = (A*Q^{-1}*t(A))^{-1} * A * Q^{-1}
// 
// log_dens_Ax_x[1]:        Value of log(pi(Ax|x)) which is equal to -0.5*log(det(A*t(A)))
//
// DIMENSIONS etc.
// ===============
// nx:                       Dimension of the GMRF
// nc:                       Number of constraints (can be 0 as well)
// mu_nonZERO:               If true then the argument 'mu' is used as the unconstrained mean of the GMRF
// e_nonZERO:                If true then the argument 'e' is used as the right-hand side of the constraint equation
// work[max(nc, nx)]:        Working array (needed only if nc > 0)
//                            * length = max(nc, nx)
//                            * used to store c = A*x and 
//                                            x - mu when evaluating the log-density
//
void
rGMRF(double *x,                    double *log_dens,
      const double *mu,             const double *Li,  const double *log_dets,     
      const double *A,              const double *e,   const double *U,  
      const double *log_dens_Ax_x, 
      const int *nx,                const int *nc,     
      const int *mu_nonZERO,        const int *e_nonZERO,    
      double *work)
{
  static int j;
  static double *xP;
  static const double *log_detsP;

  /*** Sample z ~ N(0, I) ***/
  xP = x;
  for (j = 0; j < *nx; j++){
    *xP = norm_rand();
    xP++;
  }

  /***** No constraints *****/
  /***** ============== *****/  
  if (*nc == 0){

    /*** Compute -0.5*t(z)*z to get the -0.5*(x - mu)'Q(x - mu) part of the unconstrained log-density ***/
    AK_BLAS_LAPACK::ddot2(log_dens, x, *nx);
    *log_dens *= -0.5;

    /*** Solve t(L)*v = z,  then v = t(L)^{-1}*z ~ N(0, Q^{-1}) ***/
    AK_BLAS_LAPACK::chol_solve_backward(x, Li, nx);

    /*** Compute x = mu + v, then x = mu + L^{-1}z ~ N(mu, Q^{-1}) ***/
    if (*mu_nonZERO) AK_BLAS_LAPACK::a_aPlusb(x, mu, *nx);                                     /** x = x + mu **/

    /*** Add -(n-nc)/2)*log(2*pi) + sum(log(Li[j,j])) to the log of the unconstrained density ***/
    log_detsP = log_dets;
    *log_dens += *log_detsP;
    log_detsP++;
    *log_dens += *log_detsP;
    log_detsP++;

    return;
  }

  /***** At least 1 constraint *****/
  /***** ===================== *****/

  /***** Sample from unconstrained GMRF *****/
  /***** ------------------------------ *****/

  /*** Compute -0.5*t(z)*z to get the -0.5*(x - mu)'Q(x - mu) part of the unconstrained log-density ***/
  AK_BLAS_LAPACK::ddot2(log_dens, x, *nx);
  *log_dens *= -0.5;

  /*** Solve t(L)*v = z,  then v = t(L)^{-1}*z ~ N(0, Q^{-1}) ***/
  AK_BLAS_LAPACK::chol_solve_backward(x, Li, nx);

  /*** Compute x = mu + v, then x = mu + L^{-1}z ~ N(mu, Q^{-1}) ***/
  if (*mu_nonZERO){
    AK_BLAS_LAPACK::a_aPlusb(x, mu, *nx);      /** x = x + mu **/
  }

  /***** Correction for constraints (conditioning by kriging) *****/
  /***** ---------------------------------------------------- *****/

  if (*nc == 1){
    /*** Compute c = A*x ***/
    AK_BLAS_LAPACK::ddot(work, A, x, *nx);              /** c = A*x   **/

    /*** Subtract e from c ***/
    if (*e_nonZERO){ 
      *work -= *e;            /** c = c - e **/
    }

    /*** Compute x^* = x - t(U)*c ***/
    AK_BLAS_LAPACK::a_aMinus_betabConst(x, U, work, *nx);
  }
  else{
    /*** Compute c = A*x ***/
    AK_BLAS_LAPACK::c_Ab(work, A, x, nc, nx);

    /*** Subtract e from c ***/
    if (*e_nonZERO){ 
      AK_BLAS_LAPACK::a_aMinusb(work, e, *nc);          /** c = c - e **/
    }

    /*** Compute x^* = x - t(U)*c ***/
    AK_BLAS_LAPACK::a_aMinustAb(x, U,  work, nc,  nx);
  }

  /***** Evaluate log-density *****/
  /***** -------------------- *****/

  /*** -0.5*(x-mu)'*Q*(x-mu) = -0.5*(x-mu)'*Li*t(Li)*(x-mu) part of pi(x) ***/
  if (*mu_nonZERO) AK_BLAS_LAPACK::c_aMinusb(work, x, mu, *nx);
  else             AK_BLAS_LAPACK::copyArray(work, x, *nx);
  AK_BLAS_LAPACK::a_tLa(work, Li, nx);  
  AK_BLAS_LAPACK::ddot2(log_dens, work, *nx);
  *log_dens *= -0.5;

  /*** -(n-nc)/2)*log(2*pi) + sum(log(Li[j,j])) in the log of the unconstrained density ***/
  log_detsP = log_dets;
  *log_dens += *log_detsP;
  log_detsP++;
  *log_dens += *log_detsP;
  log_detsP++;

  /*** Correct log-dens by log(pi(A*x|x)) ***/
  *log_dens += *log_dens_Ax_x;

  /*** Correct log-dens by log(determinant(W)) part of log(pi(A*x)) ***/  
  *log_dens -= *log_detsP;
  log_detsP++;

  /*** Correct log-dens by -0.5*(Amu-e)'*W^{-1}*(Amu-e) part of log(pi(A*x))  ***/
  *log_dens -= *log_detsP;  

  return;
}


/*** ======================================================================================================= ***/
/*** Compute some input arguments for 'rGMRF'                                                                ***/
/***                                                                                                         ***/
/*** ======================================================================================================= ***/
//
// OUTPUT arguments:
// =================
//   Filled always:           log_dets
//   Filled only if (nc > 0): U
//
//   log_dets[4]:    Parts needed to evaluate log-density 
//                      log_dets[0] = log(|Q|^{1/2})  = sum(log(Li[j,j])),    Q = Li*t(Li)
//                      log_dets[1] = -(nx - nc)*log(sqrt(2pi))
//                      log_dets[2] = log(|W|^{-1/2}) = -sum(log(LW[j,j])),   W = LW*t(LW)
//                      log_dets[3] = -0.5*t(A*mu - e)*W^{-1}*(A*mu - e)
//   U[nc x nx]:     Matrix U = (A*Q^{-1}*t(A))^{-1} * A * Q^{-1}
//
// INPUT arguments:
// ================
//   mu[nx]:         Mean of the unconstrained GMRF (not needed when mu_nonZERO = 0, in which case it can be set to NULL)
//   Li[LT(nx)]:     Cholesky decomposition of the unconstrained precision matrix Q,
//   A[nc x nx]:     Matrix defining the constraints (stored in column major order)
//   e[nc]:          Right-hand side of the constraint (not needed when e_nonZERO = 0, in which case it can be set to NULL)
//
// DIMENSIONS etc:
// ===============
//   nx:             Dimension of the GMRF
//   nc:             Number of constraints (can be 0 as well)
//   mu_nonZERO:     If true then the argument 'mu' is used as the unconstrained mean of the GMRF
//   e_nonZERO:      If true then the argument 'e' is used as the right-hand side of the constraint equation
//   work[]:         Working array (needed only if nc > 0)
//                          * length = nc*(nc+1)/2 + nx*nc + work(dGMRF_inputArgs)
//                                   = nc*(nc+1)/2 + nx*nc + nc
//
//                          * used to store LW = Cholesky decomposition of matrix W = A*Q^{-1}*t(A), lower triangle of matrix nc x nc
//                                          V = Q^{-1}*t(A), matrix nx x nc
//                                          working array for dGMRF_inputArgs 
//
void
rGMRF_inputArgs(double *log_dets,       double *U,        
                const double *mu,       const double *Li,  const double *A,  const double *e,
                const int *nx,          const int *nc,  
                const int *mu_nonZERO,  const int *e_nonZERO,
                double *work)
{

  static double *LW, *V, *work2;
  LW    = work;                    // LW[LT(nc)]: Cholesky decomposition of the matrix W = A*Q^{-1}*t(A) = covariance matrix of A*x, 
  V     = LW + (*nc)*(*nc+1)/2;    // V[nx x nc]: matrix Q^{-1}*t(A)
  work2 = V  + (*nx)*(*nc);        // working array for dGMRF_inputArgs

  /*** Compute U = W^{-1}*t(V) ***/
  if (*nc == 0){
    /*** Compute log_dets ***/
    dGMRF_inputArgs(log_dets, LW, U, mu, Li, A, e, nx, nc, mu_nonZERO, e_nonZERO, work2);
    return;
  }

  if (*nc == 1){
    /*** Compute log_dets and matrix V = Q{-1}*t(A) (matrix nx x nc stored in U) ***/
    dGMRF_inputArgs(log_dets, LW, U, mu, Li, A, e, nx, nc, mu_nonZERO, e_nonZERO, work2);

    /*** Solve W*U = t(V),  U is matrix 1 x nx,  t(V) is matrix 1 x nx ***/
    AK_BLAS_LAPACK::chol_solve_system(U, LW, nc, nx);    

    return;
  }

  /** nc > 1 **/
  /** ====== **/

  /*** Compute log_dets and matrix V = Q{-1}*t(A) (matrix nx x nc stored in V = work + LT(nc)) ***/
  dGMRF_inputArgs(log_dets, LW, V, mu, Li, A, e, nx, nc, mu_nonZERO, e_nonZERO, work2);

  /*** Solve W*U = t(V),  U is matrix nc x nx,  t(V) is matrix 1 x nx ***/
  AK_BLAS_LAPACK::transposition(U, V, nx, nc);
  AK_BLAS_LAPACK::chol_solve_system(U, LW, nc, nx);    

  return;
}  


/*** ======================================================================================================= ***/
/*** (Log-)density of the GMRF, possibly under a hard constraint                                             ***/
/***   See p. 38, Rue and Held (2005). Gaussian Markov Random Fields: Theory and Applications                ***/
/***                                                                                                         ***/
/***   x ~ N(mu, Q^{-1})                                                                                     ***/
/***       under A*x = e                                                                                     ***/
/***                                                                                                         ***/
/***   Factorization pi(x|Ax) = pi(x)*pi(Ax|x)/pi(Ax) is used                                                ***/
/***                                                                                                         ***/
/*** ======================================================================================================= ***/
//
// OUTPUT ARGUMENTS
// ================
// VALUE[1]:         Value of the (log-)density evaluated at x
//
// INPUT ARGUMENTS:
// ================
// x[nx]:             Point where the density is to be evaluated
// unlog[1]:          If not zero, then the value of the density is computed and not log-density
// 
// PARAMETERS OF THE GMRF
// ======================
// mu[nx]:                  Mean of the unconstrained GMRF 
//                            (not needed when mu_nonZERO = 0 in which case it can be set to NULL)
// Li[LT(nx)]:              Cholesky decomposition of the unconstrained precision matrix Q,
//                            array of length nx*(nx+1)/2
// log_dets[4]:             log_dets[0] = log(|Q|^{1/2})  = sum(log(Li[j,j])),    Q = Li*t(Li)
//                          log_dets[1] = -(nx - nc)*log(sqrt(2pi))
//                          log_dets[2] = log(|W|^{-1/2}) = -sum(log(LW[j,j])),   W = LW*t(LW)
//                          log_dets[3] = -0.5*t(A*mu - e)*W^{-1}*(A*mu - e)
//
// log_dens_Ax_x[1]:        Value of log(pi(Ax|x)) which is equal to -0.5*log(det(A*t(A)))
//
// DIMENSIONS etc.
// ===============
// nx:             Dimension of the GMRF
// nc:             Number of constraints (can be 0 as well)
// mu_nonZERO:     If true then the argument 'mu' is used as the unconstrained mean of the GMRF
// work[nx]:       Working array (needed only if nc > 0)
//                        * length = nx
//                        * used to store  x - mu when evaluating the log-density
//
void
dGMRF(double *VALUE,       
      const double *x,              const int *unlog,
      const double *mu,             const double *Li,   const double *log_dets,   
      const double *log_dens_Ax_x,  
      const int *nx,                const int *nc,              
      const int *mu_nonZERO, 
      double *work)
{
  static const double *log_detsP;

  if (*nx <= *nc){
    if (*unlog) *VALUE = 1.0;
    else        *VALUE = 0.0;
    return;
  }

  /*** -0.5*(x-mu)'*Q*(x-mu) = -0.5*(x-mu)'*Li*t(Li)*(x-mu) part of pi(x) ***/
  if (*mu_nonZERO) AK_BLAS_LAPACK::c_aMinusb(work, x, mu, *nx);
  else             AK_BLAS_LAPACK::copyArray(work, x, *nx);
  AK_BLAS_LAPACK::a_tLa(work, Li, nx);  
  AK_BLAS_LAPACK::ddot2(VALUE, work, *nx);
  *VALUE *= -0.5;

  /*** -(n-nc)/2)*log(2*pi) + sum(log(Li[j,j])) in the log of the unconstrained density ***/
  log_detsP = log_dets;
  *VALUE += *log_detsP;
  log_detsP++;
  *VALUE += *log_detsP;
  log_detsP++;

  if (*nc > 0){
    /*** Correct log-dens by log(pi(A*x|x)) ***/
    *VALUE += *log_dens_Ax_x;

    /*** Correct log-dens by log(determinant(W)) part of log(pi(A*x)) ***/  
    *VALUE -= *log_detsP;
    log_detsP++;

    /*** Correct log-dens by -0.5*(Amu-e)'*W^{-1}*(Amu-e) part of log(pi(A*x))  ***/
    *VALUE -= *log_detsP;  
  }

  if (*unlog) *VALUE = exp_AK(*VALUE);
  return;
}


/*** ======================================================================================================= ***/
/*** Compute some input arguments for 'dGMRF'                                                                ***/
/***                                                                                                         ***/
/*** ======================================================================================================= ***/
//
// OUTPUT arguments:
// =================
//   Filled always:           log_dets
//   Filled only if (nc > 0): LW, V
//
//   log_dets[4]:    Parts needed to evaluate log-density 
//                      log_dets[0] = log(|Q|^{1/2})  = sum(log(Li[j,j])),    Q = Li*t(Li)
//                      log_dets[1] = -(nx - nc)*log(sqrt(2pi))
//                      log_dets[2] = log(|W|^{-1/2}) = -sum(log(LW[j,j])),   W = LW*t(LW)
//                      log_dets[3] = -0.5*t(A*mu - e)*W^{-1}*(A*mu - e)
//   LW[LT(nc)]:     Cholesky decomposition of the matrix W = A*Q^{-1}*t(A) = covariance matrix of A*x
//   V[nx x nc]:     Matrix V = Q^{-1}*t(A)
//
// INPUT arguments:
// ================
//   mu[nx]:         Mean of the unconstrained GMRF (not needed when mu_nonZERO = 0, in which case it can be set to NULL)
//   Li[LT(nx)]:     Cholesky decomposition of the unconstrained precision matrix Q,
//   A[nc x nx]:     Matrix defining the constraints (stored in column major order)
//   e[nc]:          Right-hand side of the constraint (not needed when e_nonZERO = 0, in which case it can be set to NULL)
//
// DIMENSIONS etc:
// ===============
//   nx:             Dimension of the GMRF
//   nc:             Number of constraints (can be 0 as well)
//   mu_nonZERO:     If true then the argument 'mu' is used as the unconstrained mean of the GMRF
//   e_nonZERO:      If true then the argument 'e' is used as the right-hand side of the constraint equation
//   work[nc]:       Working array (needed only if nc > 0)
//                         * length = nc
//                         * used to store vector A*mu - e
//
void
dGMRF_inputArgs(double *log_dets,       double *LW,        double *V,
                const double *mu,       const double *Li,  const double *A,  const double *e,
                const int *nx,          const int *nc,  
                const int *mu_nonZERO,  const int *e_nonZERO,
                double *work)
{
  static int j, nx_nc, info;
  static double *log_detsP;
  static const double *LP;
  static const double _ZERO_ = 0.0;

  static double *Amu_e;

  /*** log(|Q|^{1/2}) ***/
  log_detsP = log_dets;
  LP = Li;
  *log_detsP = 0;
  for (j = *nx; j > 0; j--){
    *log_detsP += log_AK(*LP);
    LP += j;
  }
  log_detsP++;

  if (*nc == 0){
    *log_detsP = -(*nx) * _AK_LOG_SQRT_2PI;
    log_detsP++;
    *log_detsP = 0.0;
    log_detsP++;
    *log_detsP = 0.0;
    return;
  }

  Amu_e = work;                  // Amu_e[nc]:      Vector A*mu - e = mean of A*x minus right hand side of the constraint
  if (*nc == 1){
    /*** Compute A*mu ***/
    if (*mu_nonZERO) AK_BLAS_LAPACK::ddot(Amu_e, A, mu, *nx);
    else             *Amu_e = 0.0;

    /*** Solve Q*V = t(A), V is matrix nx x 1 ***/    
    AK_BLAS_LAPACK::copyArray(V, A, *nx);
    AK_BLAS_LAPACK::chol_solve_system(V, Li, nx, nc);

    /*** Compute W = A*V,  A is matrix  1 x nx, V is matrix nx x 1, store W in LW ***/
    AK_BLAS_LAPACK::ddot(LW, A, V, *nx);

    /*** Compute Cholesky decomposition of the matrix W, matrix 1 x 1 ***/
    if (*LW <= _AK_TOL_CHOL2) *LW = _AK_SQRT_TOL_CHOL2;
    else                      *LW = sqrt(*LW);

    /***  -(nx - nc)*log(sqrt(2pi)), log(|W|^{-1/2})  ***/
    *log_detsP = -((*nx) - (*nc)) * _AK_LOG_SQRT_2PI;
    log_detsP++;
    *log_detsP = -log_AK(*LW);
    log_detsP++;

    /*** Compute (A*mu - e) ***/
    if (*e_nonZERO) *Amu_e -= *e;      /** Amu - e **/

    /*** Solve LW*z2 = (A*mu - e), store z2 in Amu_e  ***/
    *Amu_e /= (*LW);

    /** -0.5*t(z2)*z2 to log_dens  **/
    *log_detsP = -0.5*(*Amu_e)*(*Amu_e);

    return;
  }

  /**** nc > 1 ****/
  /**** ====== ****/
  nx_nc = *nx * (*nc);

  /*** Compute A*mu ***/
  if (*mu_nonZERO) AK_BLAS_LAPACK::c_Ab(Amu_e, A, mu, nc, nx);      
  else             AK_BLAS_LAPACK::fillArray(Amu_e, &_ZERO_, *nc);

  /*** Solve Q*V = t(A), V is matrix nx x nc ***/    
  AK_BLAS_LAPACK::transposition(V, A, nc, nx);
  AK_BLAS_LAPACK::chol_solve_system(V, Li, nx, nc);

  /*** Compute W = A*V,  A is matrix  nc x nx, V is matrix nx x nc, store W in LW ***/
  /*** we know that A*V = A*Q^{-1}*t(A) must be symmetric                         ***/
  REprintf("Number of constraints is %d.\n", *nc);
  throw returnR("ERROR in GMRF.cpp: dGMRF_inputArgs. MULTIPLICATION A*V MUST BE IMPLEMENTED FIRST.", 1);

  /*** Compute Cholesky decomposition of the matrix W, matrix nc x nc ***/
  AK_BLAS_LAPACK::chol_dpptrf(LW, nc, &info);
  if (info) throw returnR("Error in GMRF.cpp: dGMRF_inputArgs(). Matrix A*Q^{-1}*t(A) is not positive definite.", 1);    

  /***  -(nx - nc)*log(sqrt(2pi)), log(|W|^{-1/2}) ***/
  *log_detsP = -((*nx) - (*nc)) * _AK_LOG_SQRT_2PI;
  log_detsP++;
  LP = LW;
  for (j = *nc; j > 0; j--){
    *log_detsP -= log_AK(*LP);
    LP += j;
  }
  log_detsP++;

  /*** Compute (A*mu - e) ***/
  if (*e_nonZERO) AK_BLAS_LAPACK::a_aMinusb(Amu_e, e, *nc);

  /*** Solve LW*z2 = (A*mu - e), store z2 in Amu_e ***/
  AK_BLAS_LAPACK::chol_solve_forward(Amu_e, LW, nc);

  /** -0.5*t(z2)*z2 to log_dens  **/
  AK_BLAS_LAPACK::ddot2(log_detsP, Amu_e, *nc);
  *log_detsP *= -0.5;

  return;
}


/*** ======================================================================================================= ***/
/*** (Log-)density of the GMRF, possibly under a hard constraint                                             ***/
/***   See p. 38, Rue and Held (2005). Gaussian Markov Random Fields: Theory and Applications                ***/
/***                                                                                                         ***/
/***   x ~ N(mu, Q^{-1})                                                                                     ***/
/***       under A*x = e                                                                                     ***/
/***                                                                                                         ***/
/***   Decomposition of Sigma^* = V*Lambda*t(V) is used, where                                               ***/
/***       Sigma^* = covariance matrix of x | A*x=e                                                          ***/
/***        Lambda = diag(eigenvalues of Sigma^*)                                                            ***/
/***             V = eigenvectors of Sigma^* in columns                                                      ***/
/***                                                                                                         ***/
/*** ======================================================================================================= ***/
//
// OUTPUT ARGUMENTS
// ================
// VALUE[1]:         Value of the (log-)density evaluated at x
//
// INPUT ARGUMENTS:
// ================
// x[nx]:             Point where the density is to be evaluated
// unlog[1]:          If not zero, then the value of the density is computed and not log-density
// 
// PARAMETERS OF THE GMRF
// ======================
// mu[nx]:                  Mean of the unconstrained GMRF 
//                            (not needed when mu_nonZERO = 0 in which case it can be set to NULL)
// Li[LT(nx)]:              Cholesky decomposition of the unconstrained precision matrix Q,
//                            array of length nx*(nx+1)/2
// log_dets[4]:             log_dets[0] = log(|Q|^{1/2})  = sum(log(Li[j,j])),    Q = Li*t(Li)
//                          log_dets[1] = -(nx - nc)*log(sqrt(2pi))
//                          log_dets[2] = log(|W|^{-1/2}) = -sum(log(LW[j,j])),   W = LW*t(LW)
//                          log_dets[3] = -0.5*t(A*mu - e)*W^{-1}*(A*mu - e)
//                          log_dets[4] = log(|Sigmastar|^{-1/2})
//   mustar[nx]:            Mean of the constrained GMRF
//   LiSigmastar[nx x nx]:  Decomposition of Sigmastar^{-1} (generalized inverse) such that Sigmastar^{-1} = LiSigmastar * t(LiSigmastar)
//                            * Sigmastar = Covariance matrix of the constrained GMRF (rank = nx - nc) 
//                            * LiSigmastar is NEITHER SYMMETRIC NOR LOWER TRIANGULAR!
//                            * LiSigmastar = Lambda^{-1/2} * t(V), where V = matrix with eigenvectors in columns
//                                                            Lambda^{-1/2} = diagonal matrix with sqrt(1/eigenvalues) on diagonal
//                                                                       and zeros on diagonal on places of zero eigenvalues
//                            * LiSigmastar has first nc columns equal to zero
//
// DIMENSIONS etc.
// ===============
// nx:             Dimension of the GMRF
// nc:             Number of constraints (can be 0 as well)
// mu_nonZERO:     If true then the argument 'mu' is used as the unconstrained mean of the GMRF
// work[2*nx]:     Working array (needed only if nc > 0)
//                        * length = 2*nx
//                        * used to store  x - mustar 
//                                         zscore = t(Sigma^{-1/2})*(x - mustar)'
//
void
dGMRF2(double *VALUE,       
       const double *x,              const int *unlog,
       const double *mu,             const double *Li,           const double *log_dets,   
       const double *mustar,         const double *LiSigmastar,
       const int *nx,                const int *nc,              
       const int *mu_nonZERO, 
       double *work)
{
  static const double *log_detsP;

  static double *x_mustar, *zscore; 
  x_mustar    = work;
  zscore      = x_mustar + (*nx);

  *VALUE = 0.0;
  if (*nx <= *nc){
    if (*unlog) *VALUE = 1.0;
    return;
  }
 
  if (*nc == 0){                                         /** Density of proper normal distribution **/
    if (*mu_nonZERO) AK_BLAS_LAPACK::c_aMinusb(x_mustar, x, mu, *nx);    /** x - mu                           **/
    else             AK_BLAS_LAPACK::copyArray(x_mustar, x, *nx);
    AK_BLAS_LAPACK::a_tLa(x_mustar, Li, nx);                             /** t(L)*(x - mu), where Q = L*t(L)  **/
    AK_BLAS_LAPACK::ddot2(VALUE, x_mustar, *nx);                         /** (x - mu)'*Q*(x - mu)             **/
    *VALUE *= -0.5;    

    log_detsP = log_dets;
    *VALUE += *log_detsP;
    log_detsP++;
    *VALUE += *log_detsP;

    if (*unlog) *VALUE = exp_AK(*VALUE);
    return;
  }
  else{                                                  /** Density of improper normal distribution **/

    /** Compute x - mustar **/
    AK_BLAS_LAPACK::c_aMinusb(x_mustar, x, mustar, *nx);

    /** Compute t(LiSigmastar)*(x - mustar) = t(Sigma^{-1/2})*(x - mustar)', store it in zscore              **/
    AK_BLAS_LAPACK::c_tAb(zscore, LiSigmastar, x_mustar, nx, nx);        /** This part could be improved since the first nc columns of Sigma^{-1/2} are equal to zero  **/
                                                                         /** and at this moment I multiply some zeros                                                  **/

    /** Evaluate log-density **/
    AK_BLAS_LAPACK::ddot2(VALUE, zscore, *nx);        /** (x - mustar)'*Sigma^{-}*(x - mustar)             **/
    *VALUE *= -0.5;    
    *VALUE += log_dets[1];            /** += -(nx - nc)*log(sqrt(2pi))                     **/

    /** add log(|Sigmastar|^{-1/2}) **/
    *VALUE += log_dets[4];
  }

  if (*unlog) *VALUE = exp_AK(*VALUE);
  return;
}


/*** ======================================================================================================= ***/
/*** Compute some input arguments for 'dGMRF2'                                                               ***/
/***                                                                                                         ***/
/*** ======================================================================================================= ***/
//
// OUTPUT arguments:
// =================
//   Filled always:           log_dets
//   Filled only if (nc > 0): Amu_e, LW, U, mustar, LiSigmastar
//
//   log_dets[5]:            Parts needed to evaluate log-density 
//                              log_dets[0] = log(|Q|^{1/2})  = sum(log(Li[j,j])),    Q = Li*t(Li)
//                              log_dets[1] = -(nx - nc)*log(sqrt(2pi))
//                              log_dets[2] = log(|W|^{-1/2}) = -sum(log(LW[j,j])),   W = LW*t(LW)
//                              log_dets[3] = -0.5*t(A*mu - e)*W^{-1}*(A*mu - e)
//                              log_dets[4] = log(|Sigmastar|^{-1/2})
//   mustar[nx]:             Mean of the constrained GMRF
//   LiSigmastar[nx x nx]:   Decomposition of Sigmastar^{-1} (generalized inverse) such that Sigmastar^{-1} = LiSigmastar * t(LiSigmastar)
//                            * Sigmastar = Covariance matrix of the constrained GMRF (rank = nx - nc) 
//                            * LiSigmastar is NEITHER SYMMETRIC NOR LOWER TRIANGULAR!
//                            * LiSigmastar = Lambda^{-1/2} * t(V), where V = matrix with eigenvectors in columns
//                                                            Lambda^{-1/2} = diagonal matrix with sqrt(1/eigenvalues) on diagonal
//                                                                       and zeros on diagonal on places of zero eigenvalues
//                            * LiSigmastar has first nc columns equal to zero
//
// INPUT arguments: 
// ================
//   mu[nx]:         Mean of the unconstrained GMRF (not needed when mu_nonZERO = 0, in which case it can be set to NULL)
//   Li[LT(nx)]:     Cholesky decomposition of the unconstrained precision matrix Q,
//   A[nc x nx]:     Matrix defining the constraints (stored in column major order)
//   e[nc]:          Right-hand side of the constraint (not needed when e_nonZERO = 0, in which case it can be set to NULL)
//
// DIMENSIONS etc:
// ===============
//   nx:             Dimension of the GMRF
//   nc:             Number of constraints (can be 0 as well)
//   mu_nonZERO:     If true then the argument 'mu' is used as the unconstrained mean of the GMRF
//   e_nonZERO:      If true then the argument 'e' is used as the right-hand side of the constraint equation
//   work[]:         Working array (needed only if nc > 0)
//                      * length = nc + nx*nc + nx + nx*nx + max{work(rGMRF_inputArgs), work(eigen_dspev)}
//                               = nc + nx*nc + nx + nx*nx + max{nc*(nc+1)/2 + nx*nc + nc,  3*nx}
//
//                      * used to store vector A*mu - e, vector of length nc
//                                      matrix U = (A * Q^{-1} * t(A))^{-1} * A * Q^{-1}, matrix nc x nx
//                                      Lambda = vector of eigenvalues of the matrix Sigmastar, vector of length nx
//                                      V = matrix of eigenvectors of the matrix Sigmastar, matrix nx x nx
//                                      working array for rGMRF_inputArgs or eigen_dspev
//
void
dGMRF2_inputArgs(double *log_dets,
                 double *mustar,         double *LiSigmastar,
                 const double *mu,       const double *Li,        const double *A,   const double *e,
                 const int *nx,          const int *nc,  
                 const int *mu_nonZERO,  const int *e_nonZERO,    
                 double *work)
{
  static int i, j, info;
  static double *LambdaP, *VP, *LiSigmastarP, *log_detsP;

  static const double _ZERO_      = 0;
  static const double _MINUS_ONE_ = -1;

  static double *Amu_e, *U, *Lambda, *V, *work2;
  Amu_e      = work;
  U          = Amu_e   + (*nc);
  Lambda     = U       + (*nx) * (*nc);
  V          = Lambda  + (*nx);
  work2      = V       + (*nx) * (*nx);

  if (*nx <= *nc){
    return;
  }

  /** Compute log_dets, U **/
  rGMRF_inputArgs(log_dets, U, mu, Li, A, e, nx, nc, mu_nonZERO, e_nonZERO, work2);

  /*** No constraints ***/
  /*** ============== ***/
  if (*nc == 0){
    return;
  }

  /*** Some constraints ***/
  /*** ================ ***/

  /** Compute A*mu - e **/
  if (*nc == 1){
    if (*mu_nonZERO) AK_BLAS_LAPACK::ddot(Amu_e, A, mu, *nx);
    else             *Amu_e = 0.0;
    if (*e_nonZERO) *Amu_e -= *e;
  }
  else{
    if (*mu_nonZERO) AK_BLAS_LAPACK::c_Ab(Amu_e, A, mu, nc, nx);      
    else             AK_BLAS_LAPACK::fillArray(Amu_e, &_ZERO_, *nc);
    if (*e_nonZERO) AK_BLAS_LAPACK::a_aMinusb(Amu_e, e, *nc);
  }

  /** Mean of the constrained GMRF **/
  AK_BLAS_LAPACK::c_tAb(mustar, U, Amu_e, nc, nx);                                /* t(U) * (A*mu - e)                 */
  AK_BLAS_LAPACK::a_alphaaPlusb(mustar, &_MINUS_ONE_, mu, *nx);                   /* mu - t(U) * (A*mu - e)   = mustar */

  /** Variance of the constrained GMRF **/
  AK_BLAS_LAPACK::copyLT(LiSigmastar, Li, *nx);
  AK_BLAS_LAPACK::chol_dpptri(LiSigmastar, nx, &info);                            /* Compute Q^{-1} */
  if (info) throw returnR("Error in GMRF.cpp: dGMRF2_inputArgs(). Supplied precision matrix is not positive definite.", 1);

  AK_BLAS_LAPACK::C_tAB(V, A, U, nc, nx, nx);                    /* t(A) * U                                */
  AK_BLAS_LAPACK::chol_solve_system(V, Li, nx, nx);              /* Q^{-1} * t(A) * U                       */
  AK_BLAS_LAPACK::ALT_ALTminusB(LiSigmastar, V, *nx);            /* Q^{-1} - Q^{-1} * t(A) * U  = Sigmastar */

  //Rprintf("\nConstrained mean is equal to: ");                     /***** TEMPORAR *****/
  //AK_BLAS_LAPACK::printArray(mustar, *nx);                         /***** TEMPORAR *****/
  //Rprintf("\nConstrained covariance matrix is equal to:\n");       /***** TEMPORAR *****/
  //AK_BLAS_LAPACK::printLT(LiSigmastar, *nx);                       /***** TEMPORAR *****/
  //Rprintf("\n\n");                                                 /***** TEMPORAR *****/

  /** Eigenvalue decomposition of Sigmastar, first nc elements of the computed Lambda should be zero **/
  AK_BLAS_LAPACK::eigen_dspev(LiSigmastar, Lambda, V, nx, work2, &info);
  if (info) throw returnR("Error in GMRF.cpp: dGMRF2_inputArgs(). Eigenvalue decomposition of Sigmastar failed.", 1);  

  /** Invert matrix Lambda and compute square roots of the diagonal elements (first nx - nc diagonal elements must be inverted and sqrted) **/
  //Rprintf("\nEigen values of the covariance matrix: ");  /***** TEMPORAR *****/
  //AK_BLAS_LAPACK::printArray(Lambda, *nx);               /***** TEMPORAR *****/
  //Rprintf("Eigen vectors of the covariance matrix:\n");  /***** TEMPORAR *****/
  //AK_BLAS_LAPACK::printMatrix(V, *nx, *nx);              /***** TEMPORAR *****/
  LambdaP = Lambda;
  for (j = 0; j < *nc; j++){
    *LambdaP = 0;
    LambdaP++;
  }
  for (j = *nc; j < *nx; j++){
    *LambdaP = sqrt(1/(*LambdaP));
    LambdaP++;
  }
  //Rprintf("\nDiagonal of the matrix Lambda^{-1/2}: ");    /***** TEMPORAR *****/
  //printArray(Lambda, *nx);                                /***** TEMPORAR *****/

  /** Compute Sigmastar^{-1/2} = V * Lambda^{-1/2} (matrix with first nc columns filled by zeros) **/
  LiSigmastarP = LiSigmastar;
  VP = V;
  LambdaP = Lambda;
  for (j = 0; j < *nc; j++){
    for (i = 0; i < *nx; i++){
      *LiSigmastarP = 0;
      LiSigmastarP++;
      VP++;
    }
    LambdaP++;
  }
  for (j = *nc; j < *nx; j++){
    for (i = 0; i < *nx; i++){
      *LiSigmastarP = *LambdaP * (*VP);
      LiSigmastarP++;
      VP++;
    }
    LambdaP++;
  }
  //Rprintf("Sigmastar^{-1/2}:\n");                          /***** TEMPORAR *****/
  //AK_BLAS_LAPACK::printMatrix(LiSigmastar, *nx, *nx);      /***** TEMPORAR *****/

  /** Evaluate sum(log(Lambda^{-1/2})) **/
  log_detsP = log_dets + 4;
  *log_detsP = 0;
  LambdaP = Lambda + (*nc);  
  for (j = *nc; j < *nx; j++){   
    *log_detsP += log_AK(*LambdaP);
    LambdaP++;
  }

  return;
}



/*** ============================================================================= ***/
/*** Functions allowing to sample from g(x) \propto 1 + 1/x, x \in [1/F, F], F>1   ***/
/***                                                                               ***/
/***   FUNCTIONS: dscale_norm_const, dscale, rscale, rscaleR                       ***/
/***                                                                               ***/
/*** ============================================================================= ***/

// Normalizing constant and other constants needed to get quantile function etc.
//
// OUTPUT: 
//         parD[0] = 1/F
//         parD[1] = F
//         parD[2] = log(F)
//         parD[3] = c(F) = F - 1/F + 2*log(F)
//         parD[4] = 1/F - log(F)                  
//         parD[5] = F + log(F)
//       
void
dscale_norm_const(const double *F, double *parD)
{  
  if (*F <= 1) throw returnR("Error in GMRF.cpp: dscale_norm_const(F). F must be > 1", 1);

  parD[0] = 1/(*F);
  parD[1] = *F;
  parD[2] = log(*F);
  parD[3] = *F - parD[0] + 2*parD[2];
  parD[4] = parD[0] - parD[2];
  parD[5] = *F + parD[2];

  return;
}

// Evaluate g(x) = x + log(x), its first and -second derivatives
//
//   x:     const double* ..... point x where to evalulate the function g
//   gx:    double* ........... g(x)
//   dgx:   double* ........... g'(x)
//   ddgx:  double* ........... -g''(x) (nowhere needed by this routine)
//   parD:  const double* ..... additional parameters to evaluate g
//                       parD[0] = 1/F
//                       parD[1] = F
//                       parD[2] = log(F)
//                       parD[3] = c(F) = F - 1/F + 2*log(F)
//                       parD[4] = 1/F - log(F)
//                       parD[5] = F + log(F)
//   parI:  const int* ........ additional parameters to evaluate g (not needed here)
//   what:  const int& ........ integer indicating what should be computed
//                     0 = all g(x), g'(x), g''(x)
//                     1 = only g(x)
//                     2 = only g'(x) and g''(x)
//                     3 = only g(x) and g'(x)

void
dscale(const double *x,  double *gx,  double *dgx,  double *ddgx,  const double *parD,  const int *parI,  const int &what)
{
  const double *parP = parD;
  if (*x < *parP){
    *gx = R_NegInf;
    *dgx = *ddgx = 0.0;
    return;
  }
  parP++;

  if (*x > *parP){
    *gx = R_NegInf;
    *dgx = *ddgx = 0.0;
    return;
  }

  switch (what){
  case 0:
    *gx   = *x + log_AK(*x);
    *dgx  = 1/(*x);
    *ddgx = (*dgx) * (*dgx);
    *dgx += 1; 
    return;

  case 1:
    *gx   = *x + log_AK(*x);
    return;

  case 2:
    *dgx  = 1/(*x);
    *ddgx = (*dgx) * (*dgx);
    *dgx += 1; 
    return;

  case 3:
    *gx   = *x + log_AK(*x);
    *dgx  = 1 + 1/(*x);
    return;

  default:
    throw returnR("Error in GMRF.cpp: dscale(...). Unknown what argument.", 1);
  } 
}

// Sample from g(x) \propto 1 + 1/x, x \in [1/F, F]
//
void
rscale(double *x,  const double *parD)
{
  static const int maxiter_rscale = 10;
  static const double toler_rscale = 1e-3;      /** tolerance to detect convergence of the solver_newton_raphson                      **/
  static const double zero_rscale = 1e-10;      /** tolerance to detect zero first derivatives, not needed here since g'(x) = 1 + 1/x **/

  static int iter;
  static double u, gx, dgx, ddgx, _diff, NRstep;
  static const double *parP;

  /*** Right-hand side of the equation we have to solve ***/
  u = unif_rand();
  parP = parD + 3;
  u *= (*parP);        
  parP++;
  u += (*parP);

  /*** Solve the equation using Newton-Raphson                                                       ***/
  /*** !!! We must take care of boundaries since the NRstep can lead to points outside [1/F, F]  !!! ***/
  *x = 1;                                           /** starting point **/
  dscale(x, &gx, &dgx, &ddgx, parD, NULL, 3);    

  /*** The following part of the code is based on 'solver_newton_raphson' from AK_Optim.cpp. ***/
  _diff = u - gx;
  for (iter = 0; iter < maxiter_rscale; iter++){
    if (fabs(dgx) <= zero_rscale) dgx = zero_rscale;
    NRstep = _diff/dgx;
    (*x) += NRstep;
    dscale(x, &gx, &dgx, &ddgx, parD, NULL, 3);
    if (!R_finite(gx)){    /*** we are outside [1/F, F] ***/
      if (*x < parD[0]){
        *x = parD[0];
        gx = parD[4];
        dgx = 1 + parD[1];
      }
      else{
        if (*x > parD[1]){
          *x = parD[1];
          gx = parD[5];
          dgx = 1 + parD[0];
        }
      }
    }    
    _diff = u - gx;
    if (fabs(_diff/u) <= toler_rscale) break;
  }

  return;
}



/*** ============================================================================= ***/
/*** Interface to R                                                                ***/
/*** ============================================================================= ***/

/***********************************************************************************/
// Sample n points from (constrained) GMRF
//
/***********************************************************************************/
//
// nrandom:      Number of points to be sampled
//
extern "C"{
  void
  rGMRFR(double *x,              double *log_dens,  
         const double *mu,       const double *Li,      const double *A,     const double *e,
         const int *nx,          const int *nc,         const int *nrandom,  
         const int *mu_nonZERO,  const int *e_nonZERO)
  {
    try{
      GetRNGstate(); 

      const int nwork_log_dens_Ax_x   = (*nc) * (*nc);
      const int nwork_rGMRF           = (*nx > *nc ? *nx : *nc);
      const int nwork_rGMRF_inputArgs = (*nc)*(*nc+1)/2 + (*nc)*(*nx) + (*nc); 

      int nwork_max = nwork_log_dens_Ax_x > nwork_rGMRF ? nwork_log_dens_Ax_x : nwork_rGMRF;
      if (nwork_max < nwork_rGMRF_inputArgs) nwork_max = nwork_rGMRF_inputArgs;

      const int nwork = 1 + 4 + (*nc)*(*nx) + nwork_max;
      double *work = (double*) calloc(nwork, sizeof(double));
      if (!work) throw returnR("Out of memory in GMRF.cpp: rGMRFR().", 99);

      double *log_dens_Ax_x, *log_dets, *U, *work2;
      log_dens_Ax_x = work;
      log_dets      = log_dens_Ax_x + 1;
      U             = log_dets      + 4;
      work2         = U             + (*nc)*(*nx);

      log_density_Ax_x(log_dens_Ax_x, A, nx, nc, work2);      
      rGMRF_inputArgs(log_dets, U, mu, Li, A, e, nx, nc, mu_nonZERO, e_nonZERO, work2);

      double *x_NOW        = x;
      double *log_dens_NOW = log_dens;
      for (int nr = 0; nr < *nrandom; nr++){
        rGMRF(x_NOW, log_dens_NOW, mu, Li, log_dets, A, e, U, log_dens_Ax_x, nx, nc, mu_nonZERO, e_nonZERO, work2);
        x_NOW += *nx;
        log_dens_NOW++;
      }

      free(work);
      PutRNGstate();
      return;
    }
    catch(returnR rr){
      PutRNGstate();
      return;
    }
  }
}  /*** end of extern "C" ***/


/********************************************************************************************/
// Evaluate (log-)density of (constrained) GMRF, version using decomposition of pi(x|Ax)
//
/*******************************************************************************************/
//
// nrandom:      Number of points where the (log-)density is to be evaluated
//
extern "C"{
  void
  dGMRFR(double *VALUE,  
         const double *x,        const int *unlog,
         const double *mu,       const double *Li,     const double *A,     const double *e,  
         const int *nx,          const int *nc,        const int *npoints,
         const int *mu_nonZERO,  const int *e_nonZERO)
  {
    try{
      const int nwork_log_dens_Ax_x   = (*nc) * (*nc);
      const int nwork_dGMRF           = *nx;
      const int nwork_dGMRF_inputArgs = *nc; 

      int nwork_max = nwork_log_dens_Ax_x > nwork_dGMRF ? nwork_log_dens_Ax_x : nwork_dGMRF;
      if (nwork_max < nwork_dGMRF_inputArgs) nwork_max = nwork_dGMRF_inputArgs;

      const int nwork = 1 + 4 + (*nc)*(*nc+1)/2 + (*nc)*(*nx) + nwork_max;
      double *work = (double*) calloc(nwork, sizeof(double));
      if (!work) throw returnR("Out of memory in GMRF.cpp: dGMRFR().", 99);

      double *log_dens_Ax_x, *log_dets, *LW, *V, *work2;
      log_dens_Ax_x = work;
      log_dets      = log_dens_Ax_x + 1;
      LW            = log_dets      + 4;
      V             = LW            + (*nc)*(*nc+1)/2;
      work2         = V             + (*nc)*(*nx);

      log_density_Ax_x(log_dens_Ax_x, A, nx, nc, work2);
      dGMRF_inputArgs(log_dets, LW, V, mu, Li, A, e, nx, nc, mu_nonZERO, e_nonZERO, work2);

      const double *x_NOW = x;
      double *VALUE_NOW   = VALUE;
      for (int nr = 0; nr < *npoints; nr++){
        dGMRF(VALUE_NOW, x_NOW, unlog, mu, Li, log_dets, log_dens_Ax_x, nx, nc, mu_nonZERO, work2);
        x_NOW += *nx;
        VALUE_NOW++;
      }

      free(work);
      return;
    }
    catch(returnR rr){
      return;
    }
  }
}


/********************************************************************************************/
// Evaluate (log-)density of (constrained) GMRF, version using eigenvalue decomposition
//
/*******************************************************************************************/
//
// nrandom:      Number of points where the (log-)density is to be evaluated
//
extern "C"{
  void
  dGMRF2R(double *VALUE,  
          double *mustar,         double *LiSigmastar,
          const double *x,        const int *unlog,
          const double *mu,       const double *Li,     const double *A,     const double *e,  
          const int *nx,          const int *nc,        const int *npoints,
          const int *mu_nonZERO,  const int *e_nonZERO)
  {
    try{
      const int nwork_dGMRF2           = 2 * (*nx);
      const int nwork_rGMRF_inputArgs  = (*nc)*(*nc+1)/2 + (*nc)*(*nx) + (*nc); 
      const int nwork_eigen_dspev      = 3 * (*nx);
      const int nwork_dGMRF2_inputArgs = *nc + (*nx)*(*nc) + (*nx) + (*nx)*(*nx) + (nwork_rGMRF_inputArgs > nwork_eigen_dspev ? nwork_rGMRF_inputArgs : nwork_eigen_dspev); 

      const int nwork = 5 + (nwork_dGMRF2 > nwork_dGMRF2_inputArgs ? nwork_dGMRF2 : nwork_dGMRF2_inputArgs);
      double *work = (double*) calloc(nwork, sizeof(double));
      if (!work) throw returnR("Out of memory in GMRF.cpp: dGMRF2R().", 99);

      double *log_dets, *work2;
      log_dets = work;
      work2    = log_dets + 5;

      dGMRF2_inputArgs(log_dets, mustar, LiSigmastar, mu, Li, A, e, nx, nc, mu_nonZERO, e_nonZERO, work2);

      const double *x_NOW = x;
      double *VALUE_NOW   = VALUE;
      for (int nr = 0; nr < *npoints; nr++){
        dGMRF2(VALUE_NOW, x_NOW, unlog, mu, Li, log_dets, mustar, LiSigmastar, nx, nc, mu_nonZERO, work2);
        x_NOW += *nx;
        VALUE_NOW++;
      }

      free(work);
      return;
    }
    catch(returnR rr){
      return;
    }
  }
}



/***********************************************************************************/
// Sample n points from from g(x) \propto 1 + 1/x, x \in [1/F, F]
//
/***********************************************************************************/
extern "C"{
  void
  rscaleR(double *x,  const int *n,  const double *F)
  {
    try{
      double parD[5];
      dscale_norm_const(F, parD);
  
      GetRNGstate();
      double *xP = x;
      for (int i = 0; i < *n; i++){
        rscale(xP, parD);
        xP++;
      }
      PutRNGstate();
      return;
    }
    catch(returnR rr){
      PutRNGstate();
      return;
    }
  }

}  /*** end of extern "C" ***/

}  /*** end of namespace GMRF ***/
