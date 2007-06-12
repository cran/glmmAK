/*** AK_BLAS_LAPACK.h ***/

#ifndef _AK_BLAS_LAPACK_H_
#define _AK_BLAS_LAPACK_H_

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

namespace AK_BLAS_LAPACK {

extern "C"{
  void
  transposition(double *tA,  const double *A,  const int *nrowA,  const int *ncolA);
}

void
LT2Rect(double *Rect,  const double *LT,  const int &nrow);

void
Rect2LT(double *LT,  const double *Rect,  const int &nrow);

void
fillArray(double *a,  const double *value,  const int &length);

void
fillArrayI(int *a,  const int *value,  const int &length);

void
copyArray(double *to,  const double *from,  const int &length);

void
copyLT(double *to,  const double *from,  const int &nrow);

void
copyLT_add2diag(double *to,  const double *from,  const double *eps, const int &nrow);

void
add2array(double *Q,  const double *eps,  const int &length);

void
add2diag(double *Q,  const double *eps,  const int &nrow);

void
a_aPlusb(double *a,  const double *b,  const int &length);

void
c_aPlusb(double *c,  const double *a,  const double *b,  const int &length);

void
a_aMinusb(double *a,  const double *b,  const int &length);

void
c_aMinusb(double *c,  const double *a,  const double *b,  const int &length);

void
a_alphaaPlusb(double *a,  const double *alpha, const double *b,  const int &length);

void
a_aPlus_betab(double *a,  double *b,  const double *beta,  const int &length);

void
a_aMinus_betab(double *a,  double *b,  const double *beta,  const int &length);

void
a_aMinus_betabConst(double *a,  const double *b,  const double *beta,  const int &length);

void
a_alphaa(double *a,  const double *alpha,  const int &length);

void
ddot(double *RES,  const double *a,  const double *b,  const int &length);

void
ddot(double *RES,  const double *a,  const int *b,  const int &length);

void
ddot2(double *RES,  const double *a, const int &length);

extern "C"{
  void
  a_tLa(double *a,  const double *L,  const int *na);
}

extern "C"{
  void
  a_La(double *a,  const double *L,  const int *na);
}

extern "C"{
  void
  c_ALTb(double *c,  const double *A,  const double *b,  const int *nb);
}

extern "C"{
  void
  c_Ab(double *c,  const double *A,  const double *b,  const int *nrA,  const int *ncA);
}

extern "C"{
  void
  c_tAb(double *c,  const double *A,  const double *b,  const int *nrA,  const int *ncA);
}

extern "C"{
  void
  a_aMinustAb(double *a,  const double *A,  const double *b,  const int *nrA,  const int *ncA);
}

extern "C"{
  void
  C_AB(double *C,  const double *A,  const double *B,  const int *nrA,  const int *ncA,  const int *ncB);
}

extern "C"{
  void
  C_tAB(double *C,  const double *A,  const double *B,  const int *nrA,  const int *ncA,  const int *ncB);
}

extern "C"{
  void
  C_AtB(double *C,  const double *A,  const double *B,  const int *nrA,  const int *ncA,  const int *nrB);
}

void
ALT_ALTminusB(double *A,  const double *B,  const int &nrow);

void
ALT_addb2diag(double *A,  const double *b,  const int &nrow);

void
ALT_BLTremoveRowCol(double *A,  double *a,  const double *B,  const int &nrow,  const int &iremove);

void
ALT_pp_BLTremoveRowCol(double *A,  double *a,  const double *B,  const int &nrow,  const int &iremove);

extern "C"{
  void
  ALT_BLT_min1b_minb1_plusb(double *A,  double *a,  const double *B,  const int &nrow,  const int &iremove);

  void
  ALT_pp_BLT_min1b_minb1_plusb(double *A,  double *a,  const double *B,  const int &nrow,  const int &iremove);
}

extern "C"{
  void
  chol_dpotrf(double *Q,  const int *nrow,  int *info);
}

extern "C"{
  void
  chol_dpptrf(double *Q,  const int *nrow,  int *info);
}

extern "C"{
  void
  chol_dpptri(double *Q,  const int *nrow,  int *info);
}

void
chol_dpptrfPD(double *Q,  double *Qtemp,  const int *nrow,  int *Attempt,  const int *nAttempt,  const double *eps,  int *info);

extern "C"{
  void
  chol_solve_system(double *x,  const double *L,  const int *nx,  const int *neq);
}

extern "C"{
  void
  chol_solve_forward_system(double *x,  const double *L,  const int *nx,  const int *neq);
}

extern "C"{
  void
  chol_solve_backward_system(double *x,  const double *L,  const int *nx,  const int *neq);
}

extern "C"{
  void
  chol_solve_forward(double *x,  const double *L,  const int *nx);
}

extern "C"{
  void
  chol_solve_backward(double *x,  const double *L,  const int *nx);
}

extern "C"{
  void
  eigen_dspev(double *Q,  double *evalues,  double *evectors,  const int *nrow,  double *work,  int *info);
}

void
printArray(const double *a,  const int &length);

void
printIArray(const int *a,  const int &length);

void
printLT(const double *Q,  const int &nrow);

void
printLT4R(const double *Q,  const int &nrow);

void
printMatrix(const double *Q,  const int &nrow,  const int &ncol);

}  /*** end of namespace AK_BLAS_LAPACK ***/

#endif
