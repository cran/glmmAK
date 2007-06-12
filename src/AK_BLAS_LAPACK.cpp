/*** AK_BLAS_LAPACK.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  10/11/2006
//
//    PURPOSE: Wrappers to LAPACK and BLAS  + some related routines written directly by me
//
//                 transposition:  13/11/2006
//                       LT2Rect:  11/12/2006
//                       Rect2LT:  11/12/2006
//
//                     fillArray:  14/11/2006
//                     copyArray:  10/11/2006
//                        copyLT:  10/11/2006
//               copyLT_add2diag:  10/11/2006
//
//                     add2array:  10/11/2006
//                      add2diag:  10/11/2006
//
//                      a_aPlusb:  10/11/2006
//                      c_aPlusb:  14/02/2007
//                     a_aMinusb:  13/11/2006
//                     c_aMinusb:  13/11/2006
//                 a_alphaaPlusb:  14/11/2006
//                 a_aPlus_betab:  10/11/2006
//                a_aMinus_betab:  10/11/2006
//           a_aMinus_betabConst:  13/11/2006
//                      a_alphaa:  13/11/2006
//                          ddot:  13/11/2006
//                         ddot2:  13/11/2006
//
//                         a_tLa:  13/11/2006  (wrapper to BLAS 'dtpmv')
//                          a_La:  13/11/2006  (wrapper to BLAS 'dtpmv')
//                        c_ALTb:  14/02/2007  (wrapper to BLAS 'dspmv')
//                          c_Ab:  14/11/2006  (wrapper to BLAS 'dgemv')
//                         c_tAb:  14/11/2006  (wrapper to BLAS 'dgemv')
//                    a_MinustAb:  16/11/2006  (wrapper to BLAS 'dgemv')
//
//                          C_AB:  10/11/2006  (wrapper to BLAS 'dgemm')
//                         C_tAB:  14/11/2006  (wrapper to BLAS 'dgemm')
//                         C_AtB:  17/11/2006  (wrapper to BLAS 'dgemm')
//                 ALT_ALTminusB:  14/11/2006
//                 ALT_addb2diag:  13/02/2007
//           ALT_BLTremoveRowCol:  21/11/2006
//        ALT_pp_BLTremoveRowCol:  21/11/2006
//     ALT_BLT_min1b_minb1_plusb:  21/11/2006
//  ALT_pp_BLT_min1b_minb1_plusb:  21/11/2006
//
//                   chol_dpotrf:  17/11/2006  (wrapper to LAPACK 'dpotrf')
//                   chol_dpptrf:  10/11/2006  (wrapper to LAPACK 'dpptrf')
//                   chol_dpptri:  14/11/2006  (wrapper to LAPACK 'dpptri')
//                   eigen_dspev:  14/11/2006  (wrapper to LAPACK 'dspev')
//             chol_solve_system:  10/11/2006
//     chol_solve_forward_system:  11/12/2006
//    chol_solve_backward_system:  11/12/2006
//            chol_solve_forward:  10/11/2006
//           chol_solve_backward:  10/11/2006
//
//                    printArray:  16/11/2006
//                   printIArray:  22/11/2006
//                       printLT:  16/11/2006
//                     printLT4R:  22/11/2006
//                   printMatrix:  16/11/2006
//             
//
/* ************************************************************************************************* */

#include "AK_BLAS_LAPACK.h"

namespace AK_BLAS_LAPACK {

/* ******************************************************************************** */
/* transposition:  Transposition of a general matrix                                */
/*                                                                                  */
/* ******************************************************************************** */
//
// tA[ncolA*nrowA]:  output matrix
// A[nrowA*ncolA]:   input matrix
// nrowA:            number of rows in the original matrix
// ncolA:            number of columns in the original matrix
//
extern "C"{
  void
  transposition(double *tA,  const double *A,  const int *nrowA,  const int *ncolA)
  {
    static int i, j;
    static double *tAP;
    static const double *AP;

    tAP = tA;
    for (j = 0; j < *nrowA; j++){       /** Loop over columns of t(A) **/
      AP = A + j;                       /** A[j, 0]                   **/
      for (i = 0; i < *ncolA; i++){     /** Loop over rows of t(A)    **/
        *tAP = *AP;
        AP += *nrowA;                   /** go to A[j, i+1]           **/
        tAP++;
      }
    }
 
    return;
  }
}

/* *********************************************************************************************** */
/* LT2Rect:   Copy a symmetric matrix stored in packed format into the rectangular array           */
/*                                                                                                 */
/* *********************************************************************************************** */
//
// Rect[nrow*nrow]
// LT[nrow*(nrow+1)/2]
//
void
LT2Rect(double *Rect,  const double *LT,  const int &nrow)
{
  static int j, i;
  static double *RectColP, *RectRowP, *RectDiagP;
  static const double *LTP;

  LTP       = LT;
  RectDiagP = Rect;
  for (j = 0; j < nrow; j++){
    *RectDiagP = *LTP;

    RectColP = RectDiagP + 1;
    RectRowP = RectDiagP + nrow;

    RectDiagP += nrow + 1;
    LTP++;

    for (i = j+1; i < nrow; i++){
      *RectColP = *LTP;
      *RectRowP = *LTP;
      RectColP++;
      RectRowP += nrow;
      LTP++;
    }
  }

  return;
}


/* *********************************************************************************************** */
/* Rect2LT:   Copy a symmetric matrix stored in the rectangular array into a packed array          */
/*                                                                                                 */
/* *********************************************************************************************** */
//
// LT[nrow*(nrow+1)/2]
// Rect[nrow*nrow]
//
void
Rect2LT(double *LT,  const double *Rect,  const int &nrow)
{
  static int j, i;
  static const double *RectP;
  static double *LTP;

  LTP   = LT;
  RectP = Rect;
  for (j = 0; j < nrow; j++){
    for (i = 0; i < j; i++){
      RectP++;
    }

    for (i = j; i < nrow; i++){
      *LTP = *RectP;
      LTP++;
      RectP++;
    }
  }

  return;
}


/* ******************************************************************************** */
/* fillArray:   Fill array by some value                                            */
/*                                                                                  */
/* ******************************************************************************** */
void
fillArray(double *a,  const double *value,  const int &length)
{
  static int j;
  static double *aP;

  aP = a;
  for (j = 0; j < length; j++){
    *aP = *value;
    aP++;
  }

  return;
}


/* ******************************************************************************** */
/* fillArrayI:   Fill integer array by some value                                   */
/*                                                                                  */
/* ******************************************************************************** */
void
fillArrayI(int *a,  const int *value,  const int &length)
{
  static int j;
  static int *aP;

  aP = a;
  for (j = 0; j < length; j++){
    *aP = *value;
    aP++;
  }

  return;
}


/* ******************************************************************************** */
/* copyArray:   Copy one array to another one                                       */
/*                                                                                  */
/* ******************************************************************************** */
void
copyArray(double *to,  const double *from,  const int &length)
{
  static int j;
  static const double *fromP;
  static double *toP;

  fromP = from;
  toP = to;
  for (j = 0; j < length; j++){
    *toP = *fromP;
    toP++;
    fromP++;
  }

  return;
}


/* ******************************************************************************** */
/* copyLT:   Copy lower triangle of the symmetric matrix                            */
/*                                                                                  */
/* ******************************************************************************** */
void
copyLT(double *to,  const double *from,  const int &nrow)
{
  static int i, j;
  static const double *fromP;
  static double *toP;

  fromP = from;
  toP   = to;
  for (j = 0; j < nrow; j++){
    for (i = j; i < nrow; i++){
      *toP = *fromP;
      toP++;
      fromP++;
    }
  }

  return;
}


/* ******************************************************************************** */
/* copyLT_add2diag:   Copy lower triangle of the symmetric matrix                   */
/*                    and add eps to the diagonal                                   */
/*                                                                                  */
/* ******************************************************************************** */
void
copyLT_add2diag(double *to,  const double *from,  const double *eps, const int &nrow)
{
  static int i, j;
  static const double *fromP;
  static double *toP;

  fromP = from;
  toP   = to;
  for (j = 0; j < nrow; j++){
    *toP = *fromP + (*eps);
    toP++;
    fromP++;
    for (i = j+1; i < nrow; i++){
      *toP = *fromP;
      toP++;
      fromP++;
    }
  }

  return;
}


/* ******************************************************************************** */
/* add2array:  Add a value to all elements of an array                              */
/*                                                                                  */
/* ******************************************************************************** */
void
add2array(double *Q,  const double *eps,  const int &length)
{
  static int j;
  static double *QP;  

  QP = Q;
  for (j = 0; j < length; j++){
    *QP += *eps;
    QP++;
  }

  return;
}


/* ******************************************************************************** */
/* add2diag:   Add a value to the diagonal of the symmetric matrix whose            */
/*             lower triangle is stored in Q                                        */
/*                                                                                  */
/* ******************************************************************************** */
void
add2diag(double *Q,  const double *eps,  const int &nrow)
{
  static int j;
  static double *QP;  

  QP = Q;
  for (j = nrow; j > 0; j--){
    *QP += *eps;
    QP += j;
  }

  return;
}


/*********************************************************************************************************************************/
/*********************************************************************************************************************************/
/*** BLAS, level 1                                                                                                             ***/
/*********************************************************************************************************************************/
/*********************************************************************************************************************************/

/* ******************************************************************************** */
/* a_aPlusb:  a = a + b                                                             */
/*                                                                                  */
/* ******************************************************************************** */
void
a_aPlusb(double *a,  const double *b,  const int &length)
{
  static int j;
  static double *aP;
  static const double *bP;

  aP = a;
  bP = b;
  for (j = 0; j < length; j++){
    *aP += (*bP);
    aP++;
    bP++;
  }

  return;
}


/* ******************************************************************************** */
/* c_aPlusb:  c = a + b                                                             */
/*                                                                                  */
/* ******************************************************************************** */
void
c_aPlusb(double *c,  const double *a,  const double *b,  const int &length)
{
  static int j;
  static double *cP;
  static const double *aP, *bP;

  aP = a;
  bP = b;
  cP = c;
  for (j = 0; j < length; j++){
    *cP = (*aP) + (*bP);
    aP++;
    bP++;
    cP++;
  }

  return;
}


/* ******************************************************************************** */
/* a_aMinusb:  a = a - b                                                            */
/*                                                                                  */
/* ******************************************************************************** */
void
a_aMinusb(double *a,  const double *b,  const int &length)
{
  static int j;
  static double *aP;
  static const double *bP;

  aP = a;
  bP = b;
  for (j = 0; j < length; j++){
    *aP -= (*bP);
    aP++;
    bP++;
  }

  return;
}



/* ******************************************************************************** */
/* c_aMinusb:  c = a - b                                                            */
/*                                                                                  */
/* ******************************************************************************** */
void
c_aMinusb(double *c,  const double *a,  const double *b,  const int &length)
{
  static int j;
  static double *cP;
  static const double *aP, *bP;

  aP = a;
  bP = b;
  cP = c;
  for (j = 0; j < length; j++){
    *cP = (*aP) - (*bP);
    aP++;
    bP++;
    cP++;
  }

  return;
}


/* ******************************************************************************** */
/* a_alphaaPlusb:  a = alpha*a + b                                                  */
/*                                                                                  */
/* ******************************************************************************** */
void
a_alphaaPlusb(double *a,  const double *alpha, const double *b,  const int &length)
{
  static int j;
  static double *aP;
  static const double *bP;

  aP = a;
  bP = b;
  for (j = 0; j < length; j++){
    *aP *= (*alpha);
    *aP += *bP;
    aP++;
    bP++;
  }

  return;
}


/* ******************************************************************************** */
/* a_aPlus_betab:  a = a + beta*b                                                   */
/*                 b = beta*b                                                       */
/*                                                                                  */
/* ******************************************************************************** */
void
a_aPlus_betab(double *a,  double *b,  const double *beta,  const int &length)
{
  static int j;
  static double *aP, *bP;

  aP = a;
  bP = b;
  for (j = 0; j < length; j++){
    *bP *= (*beta);
    *aP += (*bP);
    aP++;
    bP++;
  }

  return;
}


/* ******************************************************************************** */
/* a_aMinus_betab:  a = a - beta*b                                                  */
/*                  b = beta*b                                                      */
/*                                                                                  */
/* ******************************************************************************** */
void
a_aMinus_betab(double *a,  double *b,  const double *beta,  const int &length)
{
  static int j;
  static double *aP, *bP;

  aP = a;
  bP = b;
  for (j = 0; j < length; j++){
    *bP *= (*beta);
    *aP -= (*bP);
    aP++;
    bP++;
  }

  return;
}


/* ******************************************************************************** */
/* a_aMinus_betabConst:  a = a - beta*b                                             */
/*                                                                                  */
/* ******************************************************************************** */
void
a_aMinus_betabConst(double *a,  const double *b,  const double *beta,  const int &length)
{
  static int j;
  static double *aP;
  static const double *bP;

  aP = a;
  bP = b;
  for (j = 0; j < length; j++){
    *aP -= (*beta)*(*bP);
    aP++;
    bP++;
  }

  return;
}


/* ******************************************************************************** */
/* a_alphaa:  a = alpha*a                                                           */
/*                                                                                  */
/* ******************************************************************************** */
void
a_alphaa(double *a,  const double *alpha,  const int &length)
{
  static int j;
  static double *aP = a;

  for (j = 0; j < length; j++){
    *aP *= (*alpha);
    aP++;
  }

  return;
}


/* ******************************************************************************** */
/* ddot:  Compute t(a) * b                                                          */
/* IT IS OVERLOADED                                                                 */
/* ******************************************************************************** */
void
ddot(double *RES,  const double *a,  const double *b,  const int &length)
{
  static int j;
  static const double *aP, *bP;

  *RES = 0.0;
  aP = a;
  bP = b;
  for (j = 0; j < length; j++){
    *RES += (*aP)*(*bP);
    aP++;
    bP++;
  }

  return;
}


/* ******************************************************************************** */
/* ddot:  Compute t(a) * b                                                          */
/*        version with integer b                                                    */
/* IT IS OVERLOADED                                                                 */
/* ******************************************************************************** */
void
ddot(double *RES,  const double *a,  const int *b,  const int &length)
{
  static int j;
  static const double *aP;
  static const int *bP;

  *RES = 0.0;
  aP = a;
  bP = b;
  for (j = 0; j < length; j++){
    *RES += (*aP)*(*bP);
    aP++;
    bP++;
  }

  return;
}



/* ******************************************************************************** */
/* ddot2:  Compute t(a) * a                                                          */
/*                                                                                  */
/* ******************************************************************************** */
void
ddot2(double *RES,  const double *a, const int &length)
{
  static int j;
  static const double *aP;

  *RES = 0.0;
  aP = a;
  for (j = 0; j < length; j++){
    *RES += (*aP)*(*aP);
    aP++;
  }

  return;
}


/*********************************************************************************************************************************/
/*********************************************************************************************************************************/
/*** BLAS, level 2                                                                                                             ***/
/*********************************************************************************************************************************/
/*********************************************************************************************************************************/

/* ********************************************************************************************** */
/* a_tLa:  Compute a = t(L)*a,                                                                    */
/*         where L is lower triangular matrix (stored in packed form in column major order)       */
/*         wrapper to BLAS 'dtpmv'                                                                */
/*                                                                                                */
/* ********************************************************************************************** */
extern "C"{
  void
  a_tLa(double *a,  const double *L,  const int *na)
  {
    static const char *UPLO  = "L";
    static const char *TRANS = "T";
    static const char *DIAG  = "N";
    static const int INCX = 1;
    F77_CALL(dtpmv)(UPLO, TRANS, DIAG, na, L, a, &INCX);

    return;
  }
}
      

/* ********************************************************************************************** */
/* a_La:  Compute a = L*a,                                                                        */
/*        where L is lower triangular matrix (stored in packed form in column major order)        */
/*        wrapper to BLAS 'dtpmv'                                                                 */
/*                                                                                                */
/* ********************************************************************************************** */
extern "C"{
  void
  a_La(double *a,  const double *L,  const int *na)
  {
    static const char *UPLO  = "L";
    static const char *TRANS = "N";
    static const char *DIAG  = "N";
    static const int INCX = 1;
    F77_CALL(dtpmv)(UPLO, TRANS, DIAG, na, L, a, &INCX);

    return;
  }
}

/* ******************************************************************************************************** */
/* c_ALTb:  Compute c = A*b,                                                                                */
/*        where A is a symmetric matrix from which only a lower triangular is stored in column major order  */
/*        wrapper to BLAS 'dspmv'                                                                           */
/*                                                                                                          */
/* ******************************************************************************************************** */
//
// c[nb]:
// A[LT(nb)]:
// b[nb]:
//
extern "C"{
  void
  c_ALTb(double *c,  const double *A,  const double *b,  const int *nb)
  {
    static const char   *UPLO = "L";
    static const double ALPHA = 1;
    static const int    INCX  = 1;
    static const double BETA  = 0;
    static const int    INCY  = 1;
    F77_CALL(dspmv)(UPLO, nb, &ALPHA, A, b, &INCX, &BETA, c, &INCY);

    return;
  }
}


/* ******************************************************************************** */
/* c_Ab:   Compute c = A * b                                                        */
/*         wrapper to BLAS 'dgemv'                                                  */
/*                                                                                  */
/* ******************************************************************************** */
//
// All matrices stored in column major order
// c[nrA]:
// A[nrA x ncA]:
// b[ncA]:
//
extern "C"{
  void
  c_Ab(double *c,  const double *A,  const double *b,  const int *nrA,  const int *ncA)
  {
    static const char *TRANS = "N";
    static const double ALPHA = 1;
    static const int INCX = 1;
    static const int INCY = 1;
    static const double BETA = 0;
    F77_CALL(dgemv)(TRANS, nrA, ncA, &ALPHA, A, nrA, b, &INCX, &BETA, c, &INCY);

    return;    
  }
}


/* ******************************************************************************** */
/* c_tAb:   Compute c = t(A) * b                                                    */
/*         wrapper to BLAS 'dgemv'                                                  */
/*                                                                                  */
/* ******************************************************************************** */
//
// All matrices stored in column major order
// c[ncA]:
// A[nrA x ncA]:
// b[nrA]:
//
extern "C"{
  void
  c_tAb(double *c,  const double *A,  const double *b,  const int *nrA,  const int *ncA)
  {
    static const char *TRANS = "T";
    static const double ALPHA = 1;
    static const int INCX = 1;
    static const int INCY = 1;
    static const double BETA = 0;
    F77_CALL(dgemv)(TRANS, nrA, ncA, &ALPHA, A, nrA, b, &INCX, &BETA, c, &INCY);

    return;    
  }
}


/* ******************************************************************************** */
/* a_aMinustAb:   Compute a = a - t(A) * b                                          */
/*         wrapper to BLAS 'dgemv'                                                  */
/*                                                                                  */
/* ******************************************************************************** */
//
// All matrices stored in column major order
// c[ncA]:
// A[nrA x ncA]:
// b[nrA]:
//
extern "C"{
  void
  a_aMinustAb(double *a,  const double *A,  const double *b,  const int *nrA,  const int *ncA)
  {
    static const char *TRANS = "T";
    static const double ALPHA = -1;
    static const int INCX = 1;
    static const int INCY = 1;
    static const double BETA = 1;
    F77_CALL(dgemv)(TRANS, nrA, ncA, &ALPHA, A, nrA, b, &INCX, &BETA, a, &INCY);

    return;    
  }
}


/*********************************************************************************************************************************/
/*********************************************************************************************************************************/
/*** BLAS, level 3                                                                                                             ***/
/*********************************************************************************************************************************/
/*********************************************************************************************************************************/

/* ******************************************************************************** */
/* C_AB:   Compute C = A * B                                                        */
/*         wrapper to BLAS 'dgemm'                                                  */
/*                                                                                  */
/* ******************************************************************************** */
//
// All matrices stored in column major order
// C[nrA x ncB]:
// A[nrA x ncA]:
// B[ncA x ncB]:
//
extern "C"{
  void
  C_AB(double *C,  const double *A,  const double *B,  const int *nrA,  const int *ncA,  const int *ncB)
  {
    static const char *TRANS = "N";
    static const double ALPHA = 1;
    static const double BETA = 0;
    F77_CALL(dgemm)(TRANS, TRANS, nrA, ncB, ncA, &ALPHA, A, nrA, B, ncA, &BETA, C, nrA);

    return;    
  }
}


/* ******************************************************************************** */
/* C_tAB:   Compute C = t(A) * B                                                    */
/*          wrapper to BLAS 'dgemm'                                                 */
/*                                                                                  */
/* ******************************************************************************** */
//
// All matrices stored in column major order
// C[ncA x ncB]:
// A[nrA x ncA]:
// B[nrA x ncB]:
//
extern "C"{
  void
  C_tAB(double *C,  const double *A,  const double *B,  const int *nrA,  const int *ncA,  const int *ncB)
  {
    static const char *TRANS = "T";
    static const char *TRANSB = "N";
    static const double ALPHA = 1;
    static const double BETA = 0;
    F77_CALL(dgemm)(TRANS, TRANSB, ncA, ncB, nrA, &ALPHA, A, nrA, B, nrA, &BETA, C, ncA);

    return;    
  }
}


/* ******************************************************************************** */
/* C_AtB:   Compute C = A * t(B)                                                    */
/*          wrapper to BLAS 'dgemm'                                                 */
/*                                                                                  */
/* ******************************************************************************** */
//
// All matrices stored in column major order
// C[ncA x ncB]:
// A[nrA x ncA]:
// B[nrB x ncA]:
//
extern "C"{
  void
  C_AtB(double *C,  const double *A,  const double *B,  const int *nrA,  const int *ncA,  const int *nrB)
  {
    static const char *TRANS = "N";
    static const char *TRANSB = "T";
    static const double ALPHA = 1;
    static const double BETA = 0;
    F77_CALL(dgemm)(TRANS, TRANSB, nrA, nrB, ncA, &ALPHA, A, nrA, B, nrB, &BETA, C, nrA);

    return;    
  }
}


/* ******************************************************************************************* */
/* ALT_ALTminusB:  A = A - B,                                                                  */
/*                 where A is symmetric matrix stored in packed form as a lower triangle       */
/*                 and B is symmetric matrix stored unpacked                                   */
/*                                                                                             */
/* ******************************************************************************************* */
//
// A[LT(nrow)]:
// B[nrow x nrow]:
// nrow
//
void
ALT_ALTminusB(double *A,  const double *B,  const int &nrow)
{
  static int j, i;
  static double *AP;
  static const double *BP;
  
  AP = A;
  BP = B;
  for (j = 0; j < nrow; j++){
    for (i = 0; i < j; i++) BP++;
    for (i = j; i < nrow; i++){
      *AP -= *BP;
      AP++;
      BP++;
    }
  }

  return;
}


/* ******************************************************************************** */
/* add2diag:   Add a vector b to the diagonal of the symmetric matrix whose         */
/*             lower triangle is stored in A                                        */
/*                                                                                  */
/* ******************************************************************************** */
//
// A[LT(nrow)]:
// b[nrow]:
// nrow
//
void
ALT_addb2diag(double *A,  const double *b,  const int &nrow)
{
  static int j;
  static double *AP;  
  static const double *bP;

  AP = A;
  bP = b;
  for (j = nrow; j > 0; j--){
    *AP += *bP;
    AP += j;
    bP++;
  }

  return;
}


/* ******************************************************************************************* */
/* ALT_BLTremoveRowCol:  Remove a specific row and column from the symmetric matrix stored     */
/*                       in packed format as lower triangle                                    */
/*                       and return the value of that row/column + diagonal element            */
/*                                                                                             */
/* ALT_pp_BLTremoveRowCol: The same, except that resulted matrix is added to                   */
/*                         the input matrix A                                                  */
/* ******************************************************************************************* */
//
// A[LT(nrow-1)]:     Resulted symmetric matrix
// a[nrow]:           positions 1,...,nrow-1 = removed row (=column) without the diagonal element
//                    position 0 = removed diagonal element
// B[LT(nrow)]:       Input matrix
// nrow
// iremove:           Index of the row (and column) to be removed                   
//
void
ALT_BLTremoveRowCol(double *A,  double *a,  const double *B,  const int &nrow,  const int &iremove)
{
  static int i, j;
  static double *AP;
  static double *aP;
  static const double *BP;

  if (nrow <= 1 || iremove >= nrow){
    REprintf("nrow=%d,  iremove=%d", nrow, iremove);
    throw returnR("Error in AK_BLAS_LAPACK.cpp: ALT_BLTremoveRowCol.", 1);
  }

  AP = A;
  aP = a + 1;
  BP = B;

  /** j = 0, ..., iremove-1 **/
  j = 0;
  while (j < iremove){
    i = j;
    while (i < iremove){          /* i = j, ..., iremove-1 */
      *AP = *BP;
      AP++;
      BP++;
      i++;
    }
    *aP = *BP;                    /* i = iremove */
    aP++;
    BP++;
    i++;
    while (i < nrow){             /* i = iremove+1, ..., nrow-1 */
      *AP = *BP;
      AP++;
      BP++;
      i++;
    }   
    j++;
  }

  /** j = iremove **/
  i = j;
  *a = *BP;                        /* i = iremove, removed diagonal element */
  BP++;
  i++;
  while (i < nrow){                /* i = iremove+1, ..., nrow-1 */
    *aP = *BP;
    aP++;
    BP++;
    i++;
  }   
  j++;

  /** j = iremove+1, ..., nrow-1 **/
  while (j < nrow){
    i = j;
    while (i < nrow){              /* i = j (> iremove), ..., nrow-1 */
      *AP = *BP;
      AP++;
      BP++;
      i++;
    }   
    j++;
  }

  return;
}

void
ALT_pp_BLTremoveRowCol(double *A,  double *a,  const double *B,  const int &nrow,  const int &iremove)
{
  static int i, j;
  static double *AP;
  static double *aP;
  static const double *BP;

  if (nrow <= 1 || iremove >= nrow){
    REprintf("nrow=%d,  iremove=%d", nrow, iremove);
    throw returnR("Error in AK_BLAS_LAPACK.cpp: ALT_pp_BLTremoveRowCol.", 1);
  }

  AP = A;
  aP = a + 1;
  BP = B;

  /** j = 0, ..., iremove-1 **/
  j = 0;
  while (j < iremove){
    i = j;
    while (i < iremove){          /* i = j, ..., iremove-1 */
      *AP += *BP;
      AP++;
      BP++;
      i++;
    }
    *aP = *BP;                    /* i = iremove */
    aP++;
    BP++;
    i++;
    while (i < nrow){             /* i = iremove+1, ..., nrow-1 */
      *AP += *BP;
      AP++;
      BP++;
      i++;
    }   
    j++;
  }

  /** j = iremove **/
  i = j;
  *a = *BP;                        /* i = iremove, removed diagonal element */
  BP++;
  i++;
  while (i < nrow){                /* i = iremove+1, ..., nrow-1 */
    *aP = *BP;
    aP++;
    BP++;
    i++;
  }   
  j++;

  /** j = iremove+1, ..., nrow-1 **/
  while (j < nrow){
    i = j;
    while (i < nrow){              /* i = j (> iremove), ..., nrow-1 */
      *AP += *BP;
      AP++;
      BP++;
      i++;
    }   
    j++;
  }

  return;
}


/* ******************************************************************************************* */
/* ALT_BLT_min1b_minb1_plusb:  A = B[-(0,0)] - 1*b' - b*1' + b[0,0]*1*1'                       */
/*     where B is symmetric matrix and                                                         */
/*     B[-(0,0)] is submatrix of B obtained by removing a specific row and column              */
/*     b is that removed row/column without diagonal element                                   */
/*     b[0,0] is removed diagonal element                                                      */
/*                                                                                             */
/*                                                                                             */
/* ALT_pp_BLT_min1b_minb1_plusb: The same, except that resulted matrix is added to             */
/*                               the input matrix A                                            */
/* ******************************************************************************************* */
//
// A[LT(nrow-1)]:     Resulted symmetric matrix
// a[nrow]:           positions 1,...,nrow-1 = removed row (=column) without the diagonal element
//                    position 0 = removed diagonal element
// B[LT(nrow)]:       Input matrix
// nrow
// iremove:           Index of the row (and column) to be removed                   
//
extern "C"{
  void
  ALT_BLT_min1b_minb1_plusb(double *A,  double *a,  const double *B,  const int &nrow,  const int &iremove)
  {
    static int j, i;
    static double *AP;
    static const double *a1P, *a2P;    

    ALT_BLTremoveRowCol(A, a, B, nrow, iremove);

    AP = A;
    a1P = a + 1;

    /** j = 0, ..., iremove-1 **/
    j = 0;
    while (j < iremove){
      a2P = a1P;
      i = j;
      while (i < iremove){          /* i = j, ..., iremove-1 */
        *AP += *a - (*a1P) - (*a2P);
        AP++;
        a2P++;
        i++;
      }

      i++;                          /* skip the row that has been removed */

      while (i < nrow){             /* i = iremove+1, ..., nrow-1 */
        *AP += *a - (*a1P) - (*a2P);
        AP++;
        a2P++;
        i++;
      }   
      a1P++;
      j++;
    }

    j++;                           /* skip the column that has been removed */

    /** j = iremove+1, ..., nrow-1 **/
    while (j < nrow){
      a2P = a1P;
      i = j;
      while (i < nrow){            /* i = j (> iremove), ..., nrow-1 */
        *AP += *a - (*a1P) - (*a2P);
        AP++;
        a2P++;
        i++;
      }   
      a1P++;
     j++;
    } 

    return;  
  }

  void
  ALT_pp_BLT_min1b_minb1_plusb(double *A,  double *a,  const double *B,  const int &nrow,  const int &iremove)
  {
    static int j, i;
    static double *AP;
    static const double *a1P, *a2P;    

    ALT_pp_BLTremoveRowCol(A, a, B, nrow, iremove);

    AP = A;
    a1P = a + 1;

    /** j = 0, ..., iremove-1 **/
    j = 0;
    while (j < iremove){
      a2P = a1P;
      i = j;
      while (i < iremove){          /* i = j, ..., iremove-1 */
        *AP += *a - (*a1P) - (*a2P);
        AP++;
        a2P++;
        i++;
      }

      i++;                          /* skip the row that has been removed */

      while (i < nrow){             /* i = iremove+1, ..., nrow-1 */
        *AP += *a - (*a1P) - (*a2P);
        AP++;
        a2P++;
        i++;
      }   
      a1P++;
      j++;
    }

    j++;                           /* skip the column that has been removed */

    /** j = iremove+1, ..., nrow-1 **/
    while (j < nrow){
      a2P = a1P;
      i = j;
      while (i < nrow){            /* i = j (> iremove), ..., nrow-1 */
        *AP += *a - (*a1P) - (*a2P);
        AP++;
        a2P++;
        i++;
      }   
      a1P++;
     j++;
    } 

    return;  
  }
}


/*********************************************************************************************************************************/
/*********************************************************************************************************************************/
/*** LAPACK                                                                                                                    ***/
/*********************************************************************************************************************************/
/*********************************************************************************************************************************/

/* ******************************************************************************** */
/* chol_dpotrf:   Cholesky decomposition of SPD matrix Q fully stored               */
/*                in an array nrow*nrow                                             */
/*                wrapper to LAPACK 'dpptrf'                                        */
/*                                                                                  */
/*                                                                                  */
/* ******************************************************************************** */
//
// Q[nrow x nrow]:   INPUT:  matrix Q (in column major order)
//                   OUTPUT: decomposition L stored in the lower triangle (rest is unchanged), where Q = L*t(L)
// nrow:             dimension of Q
// info:             information on the success of the decomposition
//          = 0:  successful exit
//          < 0:  if INFO = -i, the i-th argument had an illegal value
//          > 0:  if INFO = i, the leading minor of order i is not
//                positive definite, and the factorization could not be
//                completed.
//                
//
extern "C"{
  void
  chol_dpotrf(double *Q,  const int *nrow,  int *info)
  {
    static const char *UPLO = "L";
    F77_CALL(dpotrf)(UPLO, nrow, Q, nrow, info);
    
    return;
  }
}


/* ******************************************************************************** */
/* chol_dpptrf:   Cholesky decomposition of SPD matrix Q stored in packed format    */
/*                as lower triangle                                                 */
/*                wrapper to LAPACK 'dpptrf'                                        */
/*                                                                                  */
/*                                                                                  */
/* ******************************************************************************** */
//
// Q[LT(nrow)]:   INPUT:  lower triangle (in column major order) or the matrix Q
//                OUTPUT: decomposition L, where Q = L*t(L)
// nrow:          dimension of Q
// info:          information on the success of the decomposition
//          = 0:  successful exit
//          < 0:  if INFO = -i, the i-th argument had an illegal value
//          > 0:  if INFO = i, the leading minor of order i is not
//                positive definite, and the factorization could not be
//                completed.
//                
//
extern "C"{
  void
  chol_dpptrf(double *Q,  const int *nrow,  int *info)
  {
    static const char *UPLO = "L";
    F77_CALL(dpptrf)(UPLO, nrow, Q, info);
    
    return;
  }
}


/* ********************************************************************************* */
/* chol_dpptrfPD:  Cholesky decomposition and attempts to make the matrix PD by      */
/*                 adding a small number to diagonal if needed                       */
/*                                                                                   */
/* ********************************************************************************* */
//
// Q[LT(nrow)]:      lower triangle of the symmetric matrix
// Qtemp[LT(nrow)]:  storage space to backup the original matrix
// nrow
// Attempt:          OUTPUT: number of attempts performed (0 if the original matrix is PD) 
// nAttempt:         maximal number of attempts to make the matrix PD               
// eps:              small number to be added to all diagonal elements if the matrix is not PD
// info:             info from the last 'chol_dpptrf' call
//
void
chol_dpptrfPD(double *Q,  double *Qtemp,  const int *nrow,  int *Attempt,  const int *nAttempt,  const double *eps,  int *info)
{
  double add[1] = {*eps};

  /** Backup the original matrix **/
  copyLT(Qtemp, Q, *nrow);

  /** Cholesky decomposition of the original matrix **/
  *Attempt = 0;
  chol_dpptrf(Q, nrow, info);

  /** If needed, inflate the diagonal **/
  while (*info && *Attempt < *nAttempt){
    copyLT_add2diag(Q, Qtemp, add, *nrow);
    (*Attempt)++;
    *add += *eps;
    chol_dpptrf(Q, nrow, info);
  }

  if (*info){  /** Restore the original matrix **/
    copyLT(Q, Qtemp, *nrow);
  }

  return;
}


/* ******************************************************************************** */
/* chol_dpptri:   Inverse of the SPD matrix from its Cholesky decomposition         */
/*                wrapper to LAPACK 'dpptri'                                        */
/*                                                                                  */
/*                                                                                  */
/* ******************************************************************************** */
//
// Q[LT(nrow)]:   INPUT:  lower triangle (in column major order) of L, where Q = L*t(L)
//                OUTPUT: lower triangle of Q^{-1}
// nrow:          dimension of Q
// info:          information on the success of the decomposition
//          = 0:  successful exit
//          < 0:  if INFO = -i, the i-th argument had an illegal value
//          > 0:  if INFO = i, the (i,i) element of the factor U or L is
//               zero, and the inverse could not be computed.                
//
extern "C"{
  void
  chol_dpptri(double *Q,  const int *nrow,  int *info)
  {
    static const char *UPLO = "L";
    F77_CALL(dpptri)(UPLO, nrow, Q, info);
    
    return;
  }
}


/* ********************************************************************************* */
/* chol_solve_system:  Solve Q*X = B                                                 */
/*                     from  Q = L*t(L)                                              */
/*                     for B being a matrix of neq right-hand sides                  */
/*                                                                                   */
/* ********************************************************************************* */
//
// x[nx*neq]:     INPUT:   right-hand sides of the equations (in columns), matrix nx x neq stored in column major order, that is matrix B
//                OUTPUT:  solutions (in columns), that is matrix Q^{-1}*B
// L[LT(nx)]:     Cholesky decomposition of Q (lower triangle only)
// nx:            number of rows of X and B
// neq:           number of columns of X and B
//
extern "C"{
  void
  chol_solve_system(double *x,  const double *L,  const int *nx,  const int *neq)
  {
    double *xP = x;

    /*** Equation 1 ***/
    chol_solve_forward(xP, L, nx);
    chol_solve_backward(xP, L, nx);

    /*** Equations 2,...,neq ***/
    for (int j = 1; j < *neq; j++){
      xP += (*nx);
      chol_solve_forward(xP, L, nx);
      chol_solve_backward(xP, L, nx);
    }

    return;
  }
}


/* ********************************************************************************* */
/* chol_solve_forward_system:  Solve L*X = B                                         */
/*                             for B being a matrix of neq right-hand sides          */
/*                                                                                   */
/* ********************************************************************************* */
//
// x[nx*neq]:     INPUT:   right-hand sides of the equations (in columns), matrix nx x neq stored in column major order, that is matrix B
//                OUTPUT:  solutions (in columns), that is matrix L^{-1}*B
// L[LT(nx)]:     lower triangular matrix (lower triangle only)
// nx:            number of rows of X and B
// neq:           number of columns of X and B
//
extern "C"{
  void
  chol_solve_forward_system(double *x,  const double *L,  const int *nx,  const int *neq)
  {
    double *xP = x;

    /*** Equation 1 ***/
    chol_solve_forward(xP, L, nx);

    /*** Equations 2,...,neq ***/
    for (int j = 1; j < *neq; j++){
      xP += (*nx);
      chol_solve_forward(xP, L, nx);
    }

    return;
  }
}


/* ********************************************************************************* */
/* chol_solve_backward_system:  Solve t(L)*X = B                                     */
/*                              for B being a matrix of neq right-hand sides         */
/*                                                                                   */
/* ********************************************************************************* */
//
// x[nx*neq]:     INPUT:   right-hand sides of the equations (in columns), matrix nx x neq stored in column major order, that is matrix B
//                OUTPUT:  solutions (in columns), that is matrix t(L)^{-1}*B
// L[LT(nx)]:     lower triangular matrix (lower triangle only)
// nx:            number of rows of X and B
// neq:           number of columns of X and B
//
extern "C"{
  void
  chol_solve_backward_system(double *x,  const double *L,  const int *nx,  const int *neq)
  {
    double *xP = x;

    /*** Equation 1 ***/
    chol_solve_backward(xP, L, nx);

    /*** Equations 2,...,neq ***/
    for (int j = 1; j < *neq; j++){
      xP += (*nx);
      chol_solve_backward(xP, L, nx);
    }

    return;
  }
}


/* ********************************************************************************* */
/* chol_solve_forward:  Solve L*x = b                                                */
/*                      for L lower triangular matrix                                */
/*                      (forward substitution)                                       */
/*                                                                                   */
/* ********************************************************************************* */
//
// x[nx]:        INPUT:   right-hand side b of the equation
//               OUTPUT:  solution, that is L^{-1}*b
// L[LT(nx)]:    lower triangle of L
// nx:           length of x
//
extern "C"{
  void
  chol_solve_forward(double *x,  const double *L,  const int *nx)
  {
    int j;
    double *xDoneP;
    const double *LP;

    double *xP = x;
    const double *LstartP = L;              /** L[0,0]  **/
    for (int i = 0; i < *nx; i++){
      LP = LstartP;                         /** L[i, 0] **/
      xDoneP = x;
      for (j = 0; j < i; j++){
        *xP -= (*LP) * (*xDoneP);
        LP += *nx - j - 1;
        xDoneP++;
      }
      *xP /= (*LP);
      LstartP++;
      xP++;
    }
 
    return;
  }
}


/* ********************************************************************************* */
/* chol_solve_backward:  Solve t(L)*x = b                                            */
/*                      for L lower triangular matrix                                */
/*                      (backward substitution)                                      */
/*                                                                                   */
/* ********************************************************************************* */
//
// x[nx]:        INPUT:   right-hand side b of the equation
//               OUTPUT:  solution, that is t(L)^{-1}*b
// L[LT(nx)]:    lower triangle of L
// nx:           length of x
//
extern "C"{
  void
  chol_solve_backward(double *x,  const double *L,  const int *nx)
  {
    int j;
    double *xDoneP;

    double *xEndP = x + (*nx) - 1;
    double *xP = xEndP;
    const double *LP = L + ((*nx)*(*nx+1))/2 - 1;       /** L[nx-1,nx-1] **/
    for (int i = *nx; i > 0; i--){
      xDoneP = xEndP;
      for (j = *nx; j > i; j--){
        *xP -= (*LP) * (*xDoneP);
        LP--;
        xDoneP--;
      }
      *xP /= (*LP);
      LP--;
      xP--;
    }

    return;
  }
}


/* **************************************************************************************************************** */
/* eigen_dspev:   Eigen values and vectors of symmetric matrix Q stored in packed format                            */
/*                wrapper to LAPACK 'dspev'                                                                         */
/*                                                                                                                  */
/*                                                                                                                  */
/* **************************************************************************************************************** */
//
// Q[LT(nrow)]:          INPUT:  lower triangle (in column major order) or the matrix Q
//                       OUTPUT: Q is overwritten by values generated during the reduction to tridiagonal form
// evalues[nrow]:        OUTPUT: computed eigenvalues (in ascending order)
// evectors[nrow*nrow]:  OUTPUT: orthonormal eigenvectors in columns
// nrow:          dimension of Q
// work[3*nrow]:  working array
// info:          information on the success of the decomposition
//          = 0:  successful exit
//          < 0:  if INFO = -i, the i-th argument had an illegal value
//          > 0:  if INFO = i, the leading minor of order i is not
//                positive definite, and the factorization could not be
//                completed.

//                
//
extern "C"{
  void
  eigen_dspev(double *Q,  double *evalues,  double *evectors,  const int *nrow,  double *work,  int *info)
  {
    static const char *UPLO = "L";
    static const char *JOBZ = "V";   /** both eigenvalues and eigenvectors **/
    F77_CALL(dspev)(JOBZ, UPLO, nrow, Q, evalues, evectors, nrow, work, info);
    
    return;
  }
}


/*********************************************************************************************************************************/
/*********************************************************************************************************************************/
/*** UTILITIES                                                                                                                 ***/
/*********************************************************************************************************************************/
/*********************************************************************************************************************************/

/* **************************************************************************************************************** */
/* printArray:                                                                                                      */
/*                                                                                                                  */
/* **************************************************************************************************************** */
void
printArray(const double *a,  const int &length)
{
  const double *aP = a;
  Rprintf("%5.5g", *aP);
  for (int j = 1; j < length; j++){
    aP++;
    Rprintf(",  %5.5g", *aP);
  }
  Rprintf("\n");

  return;
}


/* **************************************************************************************************************** */
/* printIArray:                                                                                                      */
/*                                                                                                                  */
/* **************************************************************************************************************** */
void
printIArray(const int *a,  const int &length)
{
  const int *aP = a;
  Rprintf("%d", *aP);
  for (int j = 1; j < length; j++){
    aP++;
    Rprintf(",  %d", *aP);
  }
  Rprintf("\n");

  return;
}


/* **************************************************************************************************************** */
/* printLT:                                                                                                         */
/*                                                                                                                  */
/* **************************************************************************************************************** */
void
printLT(const double *Q,  const int &nrow)
{
  int i, j, r, c, diagI;
  double x;
  for (i = 0; i < nrow; i++){
    for (j = 0; j < nrow; j++){

      if (i >= j){
        r = i;
        c = j;
      }
      else{
        r = j;
        c = i;
      }
      diagI = (c * (2*nrow - c + 1))/2;
      x = Q[diagI + r - c];
      Rprintf("%5g  ", (fabs(x) < _AK_ZERO ? 0 : x));      
    }
    Rprintf("\n");
  }

  return;
}


/* **************************************************************************************************************** */
/* printLT4R:                                                                                                       */
/*                                                                                                                  */
/* **************************************************************************************************************** */
void
printLT4R(const double *Q,  const int &nrow)
{
  int i, j, r, c, diagI;
  double x;
  Rprintf("matrix(c(");
  for (i = 0; i < nrow; i++){
    for (j = 0; j < nrow; j++){

      if (i >= j){
        r = i;
        c = j;
      }
      else{
        r = j;
        c = i;
      }
      diagI = (c * (2*nrow - c + 1))/2;
      x = Q[diagI + r - c];
      if (i > 0 || j > 0) Rprintf(", ");
      Rprintf("%5.5g", (fabs(x) < _AK_ZERO ? 0 : x));      
    }
    Rprintf("\n");
  }
  Rprintf("), nrow=%d, ncol=%d, byrow=TRUE)\n", nrow, nrow);

  return;
}


/* **************************************************************************************************************** */
/* printMatrix:                                                                                                     */
/*                                                                                                                  */
/* **************************************************************************************************************** */
void
printMatrix(const double *Q,  const int &nrow,  const int &ncol)
{
  int i, j;
  double x;
  const double *QP1, *QP2;

  QP1 = Q;
  for (i = 0; i < nrow; i++){
    QP2 = QP1;
    for (j = 0; j < ncol; j++){
      x = *QP2;
      Rprintf("%5g  ", (fabs(x) < _AK_ZERO ? 0 : x));
      QP2 += nrow;
    }
    Rprintf("\n");
    QP1++;
  }
  Rprintf("\n");

  return;
}

}  /*** end of namespace AK_BLAS_LAPACK ***/
