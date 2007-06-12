/*** Cholesky2.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  09/11/2006 (taken from MatrixLT.cpp)
//
//    PURPOSE: Routines for Cholesky decomposition of symmetric positive (semi-)definite matrices
//
/* ************************************************************************************************* */

#ifndef _CHOLESKY_TWEE_CPP_
#define _CHOLESKY_TWEE_CPP_

#include "Cholesky2.h"

/*** In all functions below A = LL' and:                                                       ***/
/***                                                                                           ***/
/***  _a[LT(_nrow)]:       Lower triangle of the input matrix stored in column major order     ***/
/***  _atemp[LT(_nrow)]:   Working array                                                       ***/
/***  _nrow:               Number of rows of the input/output matrix                           ***/
/***  _diagI[_nrow]:       Indeces of diagonal elements                                        ***/
/***                                                                                           ***/
/*** ========================================================================================  ***/
//
// cholesky0:      INPUT: A,     OUTPUT: L
// chinv0:         INPUT: L,     OUTPUT: Either A^{-1} or L^{-1} 
// chol_inv0:      INPUT: A,     OUTPUT: Either A^{-1} or L^{-1} or both            
// chol_solve0:    INPUT: A, c,  OUTPUT: L, L^{-1}c, A^{-1}c  
// choleskyPD0:
// chol_invPD0:
// chol_solvePD0:
//


/* ********************************************************************************* */
/* cholesky: Cholesky decomposition A = LL',                                         */
/*           where L is lower triangular with non-negative entries on a diagonal     */
/*           if A is positive semidefinite                                           */
/*           - matrix L is written to the array _a                                   */
/*           - function returns the rank of A                                        */
/*                                                                                   */
/* INPUT: backup: if 1, then the original matrix is not modified in the case         */
/*                it is not positive definite                                        */
/*                DEFAULT value is 0                                                 */
/*                                                                                   */
/*  * backup=0 and A is not positive semidefinite, we return A = LDL'                */
/*    where L is lower triangular with 1 on a diagonal and D is diagonal             */
/*                                                                                   */
/*  * If A is not positive semidefinite, returned rank is equal to -rank             */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
int
cholesky0(dataType *_a,  dataType *_atemp,  const int &_nrow,  const int *_diagI,  const int &backup)
{
    dataType temp;
    dataType eps, pivot;
    int  i, j, k;
    int rank, nonneg;
    dataType *a, *b;

    if (_nrow == 1){
      rank = 1*(*_a > _AK_TOL_CHOL) + (-1)*(*_a < -_AK_TOL_CHOL);
      if (rank) *_a = sqrt(*_a);
      return rank;
    }

    /*** Backup the original matrix for the case rank < _nrow ***/
    if (backup){
      a = _a;
      b = _atemp;
      for (j = 0; j < _nrow; j++){
        for (i = j; i < _nrow; i++){
          *b = *a;
          a++;
          b++;
        }
      }
    }

    /*** Find the rank ***/
    nonneg= 1;                                                   // At the beginning, I believe that A is positive definite
    eps = 0;
    a = _a; 
    for (j = 0; j < _nrow; j++) {                               // loop over a diagonal of the matrix
      if (fabs(*a) > eps)  eps = fabs(*a);                      // eps = max(A(i,i))
      a += _nrow - j;
    }
    eps *= _AK_TOL_CHOL;                                        // eps = toler * max(A(i,i))

    rank = 0;
    a = _a; 
    for (j = 0; j < _nrow; j++) {
      pivot = *a;
      if (pivot < -8*eps) nonneg= -1;
      if (fabs(pivot) < eps) {
        *a = 0;
        a += _nrow - j;                                        // skip to the next diagonal element
      }
      else  {
        a++;
        rank++;
	for (i = (j+1); i < _nrow; i++) {                     // loop over the jth column
	  temp = (*a)/pivot;
	  *a = temp;
	  _a[_diagI[i]] -= temp*temp*pivot;
	  for (k = (i+1); k < _nrow; k++) _a[_diagI[i] + k - i] -= temp*_a[_diagI[j] + k - j];
          a++;
	}
      }
    }
    rank *= nonneg;

    /*** Restore the original matrix for the case it is not positive definite and return ***/
    if (backup && rank < _nrow){
      a = _a;
      b = _atemp;
      for (j = 0; j < _nrow; j++){
        for (i = j; i < _nrow; i++){
          *a = *b;
          a++;
          b++;
        }
      }
      return rank;
    }

    /*** Compute the square root of the diagonal ***/
    a = _a; 
    if (nonneg == 1){
      for (j = 0; j < _nrow; j++){
        *a = sqrt(*a);
        a += _nrow - j;
      }        
    }

    /*** Multiply each column under the diagonal by the diagonal element ***/
    a = _a;
    for (j = 0; j < _nrow - 1; j++){
      temp = *a;                                   // diagonal element
      a++;
      for (i = j + 1; i < _nrow; i++ ){
        *a *= temp;
        a++;
      }
    }

    return rank;
}


/* ***************************************************************************************** */
/* chinv: Inversion of a symmetric, positive definite matrix from its Cholesky               */
/*        decomposition (previously computed)                                                */
/*                                                                                           */
/*        - INPUT: array _a contains L (Cholesky decomposition of the original matrix)       */
/*        - OUTPUT: array _a contains either A^{-1} or L^{-1}                                */ 
/*                                                                                           */
/* * onlyCholInv: if != 0, inversion of the L matrix is returned instead                     */
/*                of the inversion of the original matrix                                    */
/*                (default value is 0)                                                       */
/*                                                                                           */
/* ----------------------------------------------------------------------------------------- */
template <typename dataType>
void
chinv0(dataType *_a,  dataType *_atemp,  const int &_nrow,  const int *_diagI,  const int &onlyCholInv)
{
  dataType temp;
  int i, j, k;
  dataType *a, *b;

  if (_nrow == 1){
    *_a = (onlyCholInv ? 1/(*_a) : 1/((*_a)*(*_a)));
    return;
  }

  /*** Divide each column under the diagonal by the diagonal element. ***/
  a = _a;
  for (j = 0; j < _nrow - 1; j++){
    temp = *a;                                   // diagonal element
    if (*a != 0){
      a++;
      for (i = j + 1; i < _nrow; i++ ){
        *a /= temp;
        a++;
      }
    }
    else{
      a += _nrow - j;                            // skip to the next diagonal element
    }
  }

  /*** (Almost) invert the cholesky in the lower triangle.                                                           ***/
  /*** Do not multiply the rows in the lower strict triangle by the inverse of the diagonal element on that row.     ***/
  /*** This corresponds to L^{-1} where L had ones on the diagonal.                                                  ***/
  a = _a;
  for (j = 0; j < _nrow; j++){                                  // loop over columns
    if (*a > 0) {                                               // if A[j,j] > 0
      *a = 1/(*a);                                              // this line inverts a diagonal
      a++;
      for (i= (j+1); i < _nrow; i++){                           // loop over the rows of the j-th column
        *a = -(*a);
  	for (k = 0; k < j; k++){                                // sweep operator
          _a[_diagI[k] + i - k] += (*a) * _a[_diagI[k] + j - k];     
        }
        a++;
      }
    }
    else{
      a += _nrow - j;                                           // skip to the next diagonal element
    }
  }

  if (onlyCholInv){
    a = _a;
    for (i = 0; i < _nrow; i++){                 // loop over rows (each of them must be multiplied by the corresponding diagonal element)
      b = _a + i;                                // pointer to A(i, 0)
      if (*a == 0){
        for (j = 0 ; j < i; j++){
          *b = 0;
          b += _nrow - j - 1;                    // skip to A(i, j+1)
        }
      }
      else{
        for (j = 0 ; j < i; j++){
          *b *= (*a);
          b += _nrow - j - 1;                    // skip to A(i, j+1)
        }
      }
      a += _nrow - i;                            // skip to the next diagonal element
    }
    return;
  }


  /*** Lower triangle now contains inverse of cholesky.                                        ***/                   
  /*** Calculate (L^{-1})'DL^{-1} (inverse of cholesky decomp. process),                       ***/
  /***   where D is a diagonal matrix with 1/(diag(L)*diag(L)) and L had ones on a diagonal.   ***/
  /*** This will be the inverse of the original matrix.                                        ***/

  /*** Store L^{-1}D matrix in a temporary array.   ***/
  a = _a;
  b = _atemp;
  for (j = 0; j < _nrow; j++){
    *b = (*a) * (*a);
    a++;
    b++;
    for (i = j + 1; i < _nrow; i++){
      *b = *a;
      a++;
      b++;
    }
  }

  /*** Compute L^{-1}DL^{-1}.   ***/
  a = _a;
  b = _atemp;
  for (j = 0; j < _nrow; j++) {                              // loop over columns
    if (*a == 0){                                            // singular column -> fill it by zeros (diagonal included)
      for (i = j; i < _nrow; i++){
        *a = 0;
        a++;
      }
      b += _nrow - j;
    }
    else {
      *a = *b;                                              // initialize the jth diagonal element
      if (_nrow - 1 > j){
        i = _nrow - 1;
        _a[_diagI[j] + i - j] = _atemp[_diagI[i]] * _atemp[_diagI[j] + i - j];        // last row of the jth column
        for (k = i - 1; k > j; k--){                                                  // initialize rows j + 1, ..., n - 2 of the jth column
          _a[_diagI[j] + k - j] = _a[_diagI[j] + i - j] * _atemp[_diagI[k] + i - k];
        }
        _a[_diagI[j]] += _a[_diagI[j] + i - j] * _atemp[_diagI[k] + i - j];           // update the jth diagonal element       
      }
      for (i = _nrow - 2; i > j; i--){
        temp = _atemp[_diagI[i]] * _atemp[_diagI[j] + i - j];
        _a[_diagI[j] + i - j] += temp;                                   // update the (i,j)th element
        for (k = i - 1; k >= j; k--){                                    // update rows j , ..., i - 1 of the jth column
          _a[_diagI[j] + k - j] += temp * _atemp[_diagI[k] + i - k];
        }
      }
      a += _nrow - j;
      b += _nrow - j;
    }
  }

  return;
}


/* ********************************************************************************* */
/* chol_inv: Inversion of a symmetric, positive definite matrix                      */
/*           via Cholesky decomposition A = LL'                                      */
/*                                                                                   */
/*  what: 0 = compute only A^{-1} (DEFAULT)                                          */
/*        1 = compute only L^{-1}                                                    */
/*        2 = compute both L^{-1} and A^{-1}                                         */
/*  Li:   NULL (DEFAULT) or array of length _nrow*(_nrow+1)/2                        */
/*                                                                                   */
/*   - INPUT: _a = original matrix A                                                 */
/*            Li = NULL or array of length _nrow*(_nrow+1)/2                         */
/*                                                                                   */
/*   - OUTPUT: if what = 0: _a = A^{-1}                                              */
/*                          Li not modified                                          */
/*             if what = 1: _a = L^{-1}                                              */
/*                          Li not modified                                          */
/*             if what = 2: _a = A^{-1}                                              */
/*                          Li = L^{-1}                                              */
/*                                                                                   */
/* * Function returns the rank of the matrix as described in 'cholesky' function.    */
/* * If rank < _nrow, neither array _a nor array Li are modified                     */
/*                                                                                   */
/* * This is in fact combination of 'cholesky' and 'chinv' functions,                */
/*   where checks for positive diagonal in L are removed since the                   */
/*   function inverts only positive definite matrices                                */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
int
chol_inv0(dataType *_a,  dataType *_atemp,  const int &_nrow,  const int *_diagI,  const int &what,  dataType *Li)
{
  dataType temp;
  dataType eps, pivot;
  int  i, j, k;
  int rank, nonneg;
  dataType *a, *b;

  if (_nrow == 1){
    rank = 1*(*_a > _AK_TOL_CHOL) + (-1)*(*_a < -_AK_TOL_CHOL);
    if (rank){
      switch (what){
      case 0:
        *_a = 1/(*_a);
        break;
      case 1:
        *_a = 1/sqrt(*_a);
        break;
      case 2:
        *_a = 1/(*_a);
        *Li = sqrt(*_a);
        break;
      default:
        throw returnR("Cholesky2.cpp: chol_inv0(..., what, Li) error. Unknown what argument", 1);
      }
    }
    return rank;
  }

  /*** Backup the original matrix for the case rank < _nrow ***/
  a = _a;
  b = _atemp;
  for (j = 0; j < _nrow; j++){
    for (i = j; i < _nrow; i++){
      *b = *a;
      a++;
      b++;
    }
  }

  /*** Find the rank ***/
  nonneg= 1;                                                   // At the beginning, I believe that A is positive definite
  eps = 0;
  a = _a; 
  for (j = 0; j < _nrow; j++) {                               // loop over a diagonal of the matrix
    if (fabs(*a) > eps)  eps = fabs(*a);                      // eps = max(A(i,i))
    a += _nrow - j;
  }
  eps *= _AK_TOL_CHOL;                                        // eps = toler * max(A(i,i))

  rank = 0;
  a = _a; 
  for (j = 0; j < _nrow; j++) {
    pivot = *a;
    if (pivot < -8*eps) nonneg= -1;
    if (fabs(pivot) < eps) {
      *a = 0;
      a += _nrow - j;                                        // skip to the next diagonal element
    }
    else  {
      a++;
      rank++;
      for (i = (j+1); i < _nrow; i++) {                     // loop over the jth column
        temp = (*a)/pivot;
        *a = temp;
        _a[_diagI[i]] -= temp*temp*pivot;
        for (k = (i+1); k < _nrow; k++) _a[_diagI[i] + k - i] -= temp*_a[_diagI[j] + k - j];
        a++;
      }
    }
  }
  rank *= nonneg;

  /*** Restore the original matrix for the case it is not positive definite and return ***/
  if (rank < _nrow){
    a = _a;
    b = _atemp;
    for (j = 0; j < _nrow; j++){
      for (i = j; i < _nrow; i++){
        *a = *b;
        a++;
        b++;
      }
    }
    return rank;
  }

  /*** Compute the square root of the diagonal ***/
  a = _a; 
  for (j = 0; j < _nrow; j++){
    *a = sqrt(*a);
    a += _nrow - j;
  }        

  /*** (Almost) invert the cholesky in the lower triangle.                                                           ***/
  /*** Do not multiply the rows in the lower strict triangle by the inverse of the diagonal element on that row.     ***/
  /*** This corresponds to L^{-1} where L had ones on the diagonal.                                                  ***/
  a = _a;
  for (j = 0; j < _nrow; j++){                                  // loop over columns
    *a = 1/(*a);                                                // this line inverts a diagonal
    a++;
    for (i= (j+1); i < _nrow; i++){                             // loop over the rows of the j-th column
      *a = -(*a);
      for (k = 0; k < j; k++){                                  // sweep operator
        _a[_diagI[k] + i - k] += (*a) * _a[_diagI[k] + j - k];     
      }
      a++;
    }
  }

  switch (what){
  case 0:
    break;
  case 1:             /*** Invert cholesky (store it in _a) ***/
    a = _a;
    for (i = 0; i < _nrow; i++){                 // loop over rows (each of them must be multiplied by the corresponding diagonal element)
      b = _a + i;                                // pointer to A(i, 0)
      for (j = 0 ; j < i; j++){                  // loop over columns of the i-th row
        *b *= (*a);
        b += _nrow - j - 1;                      // skip to A(i, j+1)
      }
      a += _nrow - i;                            // skip to the next diagonal element
    }
    return rank;
  case 2:             /*** Invert cholesky (store it in Li) ***/
    a = _a;
    b = Li;
    for (j = 0; j < _nrow; j++){
      for (i = j; i < _nrow; i++){
        *b = *a;
        a++;
        b++;
      }
    }
    a = Li;
    for (i = 0; i < _nrow; i++){                 // loop over rows (each of them must be multiplied by the corresponding diagonal element)
      b = Li + i;                                // pointer to A(i, 0)
      for (j = 0 ; j < i; j++){                  // loop over columns of the i-th row
        *b *= (*a);
        b += _nrow - j - 1;                      // skip to A(i, j+1)
      }
      a += _nrow - i;                            // skip to the next diagonal element
    }
    break;
  default:
    throw returnR("Cholesky2.cpp: chol_inv0(..., what, Li) error. Unknown what argument", 1);
  }
  

  /*** Lower triangle now contains inverse of cholesky.                                        ***/                   
  /*** Calculate (L^{-1})'DL^{-1} (inverse of cholesky decomp. process),                       ***/
  /***   where D is a diagonal matrix with 1/(diag(L)*diag(L)) and L had ones on a diagonal.   ***/
  /*** This will be the inverse of the original matrix.                                        ***/

  /*** Store L^{-1}D matrix in a temporary array.   ***/
  a = _a;
  b = _atemp;
  for (j = 0; j < _nrow; j++){
    *b = (*a) * (*a);
    a++;
    b++;
    for (i = j + 1; i < _nrow; i++){
      *b = *a;
      a++;
      b++;
    }
  }

  /*** Compute L^{-1}DL^{-1}.   ***/
  a = _a;
  b = _atemp;
  for (j = 0; j < _nrow; j++) {                              // loop over columns
    *a = *b;                                                 // initialize the jth diagonal element
    if (_nrow - 1 > j){
      i = _nrow - 1;
      _a[_diagI[j] + i - j] = _atemp[_diagI[i]] * _atemp[_diagI[j] + i - j];        // last row of the jth column
      for (k = i - 1; k > j; k--){                                                  // initialize rows j + 1, ..., n - 2 of the jth column
        _a[_diagI[j] + k - j] = _a[_diagI[j] + i - j] * _atemp[_diagI[k] + i - k];
      }
      _a[_diagI[j]] += _a[_diagI[j] + i - j] * _atemp[_diagI[k] + i - j];           // update the jth diagonal element       
    }
    for (i = _nrow - 2; i > j; i--){
      temp = _atemp[_diagI[i]] * _atemp[_diagI[j] + i - j];
      _a[_diagI[j] + i - j] += temp;                                   // update the (i,j)th element
      for (k = i - 1; k >= j; k--){                                    // update rows j , ..., i - 1 of the jth column
        _a[_diagI[j] + k - j] += temp * _atemp[_diagI[k] + i - k];
      }
    }
    a += _nrow - j;
    b += _nrow - j;
  }

  return rank;
}


/* ********************************************************************************* */
/* chol_solve: Solve Ax = c, where A is a symmetric, positive definite matrix        */
/*             via Cholesky decomposition A = LL'                                    */
/*                                                                                   */
/*   - INPUT: _a = original matrix A (left-hand side of the equation)                */
/*             c = right-hand side of the equation (array of length _nrow)           */
/*   - OUTPUT:      _a = L                                                           */
/*              c(OUT) = L^{-1}*c(IN)                                                */
/*                   x = t(L^{-1})*c(OUT)                                            */
/*                     = A^{-1}*c(IN)     = solution (array of length_nrow)          */
/*                                                                                   */
/* * Function returns the rank of the matrix A as described in 'cholesky' function.  */
/* * If rank < _nrow, arrays _a, c, x are not modified                               */
/*                                                                                   */
/* Generalization added on 13/09/2006:                                               */
/* -----------------------------------                                               */
/*     add_u: DEFAULT value is 0                                                     */
/*   solve_u: DEFAULT value is 0                                                     */
/*         u: DEFAULT value is NULL                                                  */
/*                                                                                   */
/*   with add_u = 1:                                                                 */
/*     * allows sampling from N(m, A^{-1}), where m = A^{-1}*c                       */
/*       that is where m is the solution to Am = c                                   */
/*                                                                                   */
/*     * on OUTPUT:     _a = L                                                       */
/*                  c(OUT) = L^{-1}*c(IN) + u                                        */
/*                  x = t(L^{-1})*c(OUT)                                             */
/*                    = A^{-1}*c(IN) + t(L^{-1})*u                                   */
/*                    = m + t(L^{-1})*u                                              */
/*                                                                                   */
/*                    that is if u is sampled from N(0, eye)                         */
/*                    then x is sampled from N(m, A^{-1})                            */
/*                                                                                   */
/*  with solve_u = 1 && add_u = 1:                                                   */
/*     * additionally solves t(L)*x2 = u                                             */
/*       and returns it in u                                                         */
/*     * that is on OUTPUT:     u(OUT) = t(L^{-1})*u(IN)                             */
/*                                     = x(OUT) - m                                  */
/*     * this will be useful for computation of the density evaluated in x           */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
int
chol_solve0(dataType *_a,      dataType *_atemp,    const int &_nrow,  const int *_diagI,
            dataType *c,       dataType *x, 
            const int &add_u,  const int &solve_u,  dataType *u)
{
  dataType temp;
  dataType eps, pivot;
  int  i, j, k;
  int rank, nonneg;
  dataType *a, *b, *cP, *xP, *xDoneP;
  const int *dI;

  if (_nrow == 1){
    rank = 1*(*_a > _AK_TOL_CHOL) + (-1)*(*_a < -_AK_TOL_CHOL);
    if (rank){
      *_a = sqrt(*_a);
      *c /= *_a;
      if (add_u){
        *c += *u;
        *x = (*c)/(*_a);
        if (solve_u){
          *u /= *_a;
        }
      }
      else{  
        *x = (*c)/(*_a);
      }      
    }
    return rank;
  }

  /*** Backup the original matrix for the case rank < _nrow ***/
  a = _a;
  b = _atemp;
  for (j = 0; j < _nrow; j++){
    for (i = j; i < _nrow; i++){
      *b = *a;
      a++;
      b++;
    }
  }

  /*** Find the rank ***/
  nonneg= 1;                                                   // At the beginning, I believe that A is positive definite
  eps = 0;
  a = _a; 
  for (j = 0; j < _nrow; j++) {                               // loop over a diagonal of the matrix
    if (fabs(*a) > eps)  eps = fabs(*a);                      // eps = max(A(i,i))
    a += _nrow - j;
  }
  eps *= _AK_TOL_CHOL;                                        // eps = toler * max(A(i,i))

  rank = 0;
  a = _a; 
  for (j = 0; j < _nrow; j++) {
    pivot = *a;
    if (pivot < -8*eps) nonneg= -1;
    if (fabs(pivot) < eps) {
      *a = 0;
      a += _nrow - j;                                        // skip to the next diagonal element
    }
    else  {
      a++;
      rank++;
      for (i = (j+1); i < _nrow; i++) {                     // loop over the jth column
        temp = (*a)/pivot;
        *a = temp;
        _a[_diagI[i]] -= temp*temp*pivot;
        for (k = (i+1); k < _nrow; k++) _a[_diagI[i] + k - i] -= temp*_a[_diagI[j] + k - j];
        a++;
      }
    }
  }
  rank *= nonneg;

  /*** Restore the original matrix for the case it is not positive definite and return ***/
  if (rank < _nrow){
    a = _a;
    b = _atemp;
    for (j = 0; j < _nrow; j++){
      for (i = j; i < _nrow; i++){
        *a = *b;
        a++;
        b++;
      }
    }
    return rank;
  }

  /*** Compute the square root of the diagonal ***/
  a = _a; 
  for (j = 0; j < _nrow; j++){
    *a = sqrt(*a);
    a += _nrow - j;
  }        

  /*** Multiply each column under the diagonal by the diagonal element ***/
  a = _a;
  for (j = 0; j < _nrow - 1; j++){
    temp = *a;                                   // diagonal element
    a++;
    for (i = j + 1; i < _nrow; i++ ){
      *a *= temp;
      a++;
    }
  }

  /*** Solve Ly = c (forward substitution) ***/
  /*** write solution to c                 ***/
  cP = c;
  for (i = 0; i < _nrow; i++){
    a = _a + i;
    xDoneP = c;
    for (j = 0; j < i; j++){
      *cP -= (*a) * (*xDoneP);
      a += _nrow - j - 1;
      xDoneP++;
    }
    *cP /= (*a);
    cP++;
  }

  /*** Add u if required ***/
  if (add_u){
    cP = c;
    xP = u;
    for (i = 0; i < _nrow; i++){
      *cP += *xP;
      cP++;
      xP++;
    }
  }

  /*** Solve L'x = y or L'x = y + u (backward substitution) ***/
  /*** write solution to x                                  ***/
  dI = _diagI + _nrow;
  cP = c + _nrow;
  xP = x + _nrow;
  for (j = _nrow-1; j >= 0; j--){
    dI--;
    cP--;
    xP--;

    a = _a + (*dI);
    temp = *a;
    *xP = *cP;
    for (i = j+1; i < _nrow; i++){
      a++;
      *xP -= (*xDoneP) * (*a);
      xDoneP++;
    }
    *xP /= temp;
    xDoneP = xP;
  }

  /*** Solve L'x2 = u (backward substitution)               ***/
  /*** write solution to u                                  ***/
  if (add_u && solve_u){
    dI = _diagI + _nrow;
    cP = u + _nrow;
    xP = u + _nrow;
    for (j = _nrow-1; j >= 0; j--){
      dI--;
      cP--;
      xP--;

      a = _a + (*dI);
      temp = *a;
      *xP = *cP;
      for (i = j+1; i < _nrow; i++){
        a++;
        *xP -= (*xDoneP) * (*a);
        xDoneP++;
      }
      *xP /= temp;
      xDoneP = xP;
    }
  }

  return rank;
}


/* ********************************************************************************* */
/* choleskyPD: Call to 'cholesky'. If the matrix is not positive definite,           */
/*             this function tries to make it PD by adding a small multiple          */
/*             of the minimal positive diagonal element to all diagonal elements.    */
/*                                                                                   */
/* INPUT: nAttempt:  maximal number of attempts to make the matrix PD                */
/*             eps:  multiple of the minimal abs(diagonal element) to be added       */
/*                   to all diagonal elements                                        */
/*                                                                                   */
/* OUTPUT: Attempt:  number of attempts performed (0 if the original matrix is PD)   */
/*                                                                                   */
/* RETURN: Rank of the last matrix with inflated diagonal                            */
/*         rank = _nrow indicates that the last attempt to make the matrix PD        */
/*                was successfull                                                    */
/*                                                                                   */
/*                                                                                   */
/* ********************************************************************************* */
template <typename dataType>
int
choleskyPD0(dataType *_a,   dataType *_atemp,     const int &_nrow,    const int *_diagI,
            int &Attempt,   const int &nAttempt,  const dataType &eps)
{
  int i, rank;
  dataType *a, minDiag, add;
  const dataType *aC;

  Attempt = 0;
  rank = cholesky0(_a, _atemp, _nrow, _diagI, 1);
  if (rank < _nrow){

    /*** Find minimum of the positive entries on the diagonal ***/
    minDiag = 0;
    aC = _a;
    for (i = 0; i < _nrow; i++){
      if (*aC > _AK_ZERO && *aC < minDiag) minDiag = *aC;
      aC += _nrow - i;
    }
    add = (minDiag < _AK_ZERO ? eps : eps*minDiag);

    /*** Try to make the matrix positive definite ***/
    while (rank < _nrow && Attempt < nAttempt){
      Attempt ++;
      a = _a;
      for (i = 0; i < _nrow; i++){
        *a += add;
        a += _nrow - i;
      }
      rank = cholesky0(_a, _atemp, _nrow, _diagI, 1);
    }

  }

  return rank;
}


/* ********************************************************************************* */
/* chol_invPD: Call to 'chol_inv'. If the matrix is not positive definite,           */
/*             this function tries to make it PD by adding a small multiple          */
/*             of the minimal positive diagonal element to all diagonal elements.    */
/*                                                                                   */
/* INPUT:     what:  see 'chol_inv'                                                  */
/*              Li:  see 'chol_inv'                                                  */
/*        nAttempt:  maximal number of attempts to make the matrix PD                */
/*             eps:  multiple of the minimal abs(diagonal element) to be added       */
/*                   to all diagonal elements                                        */
/*                                                                                   */
/* OUTPUT: Attempt:  number of attempts performed (0 if the original matrix is PD)   */
/*             ETC:  see 'chol_inv'                                                  */
/*                                                                                   */
/* RETURN: Rank of the last matrix with inflated diagonal                            */
/*         rank = _nrow indicates that the last attempt to make the matrix PD        */
/*                was successfull                                                    */
/*                                                                                   */
/*                                                                                   */
/* ********************************************************************************* */
template <typename dataType>
int
chol_invPD0(dataType *_a,         dataType *_atemp,     const int &_nrow,  const int *_diagI,
            int &Attempt,         const int &what,      dataType *Li, 
            const int &nAttempt,  const dataType &eps)
{
  int i, rank;
  dataType *a, minDiag, add;
  const dataType *aC;

  Attempt = 0;
  rank = chol_inv0(_a, _atemp, _nrow, _diagI, what, Li);
  if (rank < _nrow){

    /*** Find minimum of the diagonal ***/
    minDiag = 0;
    aC = _a;
    for (i = 0; i < _nrow; i++){
      if (*aC > _AK_ZERO && *aC < minDiag) minDiag = *aC;
      aC += _nrow - i;
    }
    add = (minDiag < _AK_ZERO ? eps : eps*minDiag);

    /*** Try to make the matrix positive definite ***/
    while (rank < _nrow && Attempt < nAttempt){
      Attempt ++;
      a = _a;
      for (i = 0; i < _nrow; i++){
        *a += add;
        a += _nrow - i;
      }
      rank = chol_inv0(_a, _atemp, _nrow, _diagI, what, Li);      
    }

  }

  return rank;
}


/* ********************************************************************************* */
/* chol_solvePD: Call to 'chol_solve'. If the matrix is not positive definite,       */
/*             this function tries to make it PD by adding a small multiple          */
/*             of the minimal positive diagonal element to all diagonal elements.    */
/*                                                                                   */
/* INPUT:        c:  see 'chol_solve'                                                */
/*        nAttempt:  maximal number of attempts to make the matrix PD                */
/*             eps:  multiple of the minimal abs(diagonal element) to be added       */
/*                   to all diagonal elements                                        */
/*           add_u:  see 'chol_solve'                                                */
/*         solve_u:  see 'chol_solve'                                                */
/*               u:  see 'chol_solve'                                                */                          
/*                                                                                   */
/* OUTPUT: Attempt:  number of attempts performed (0 if the original matrix is PD)   */
/*             ETC:  see 'chol_solve'                                                */
/*                                                                                   */
/* RETURN: Rank of the last matrix with inflated diagonal                            */
/*         rank = _nrow indicates that the last attempt to make the matrix PD        */
/*                was successfull                                                    */
/*                                                                                   */
/*                                                                                   */
/* ********************************************************************************* */
template <typename dataType>
int
chol_solvePD0(dataType *_a,         dataType *_atemp,     const int &_nrow,   const int *_diagI,
              int &Attempt,         dataType *c,          dataType *x, 
              const int &nAttempt,  const dataType &eps,
              const int &add_u,     const int &solve_u,   dataType *u)
{
  int i, rank;
  dataType *a, minDiag, add;
  const dataType *aC;

  Attempt = 0;
  rank = chol_solve0(_a, _atemp, _nrow, _diagI, c, x, add_u, solve_u, u);
  if (rank < _nrow){

    /*** Find minimum of the positive entries on the diagonal ***/
    minDiag = 0;
    aC = _a;
    for (i = 0; i < _nrow; i++){
      if (*aC > _AK_ZERO && *aC < minDiag) minDiag = *aC;
      aC += _nrow - i;
    }
    add = (minDiag < _AK_ZERO ? eps : eps*minDiag);

    /*** Try to make the matrix positive definite ***/
    while (rank < _nrow && Attempt < nAttempt){
      Attempt ++;
      a = _a;
      for (i = 0; i < _nrow; i++){
        *a += add;
        a += _nrow - i;
      }
      rank = chol_solve0(_a, _atemp, _nrow, _diagI, c, x, add_u, solve_u, u);
    }

  }

  return rank;
}

#endif
