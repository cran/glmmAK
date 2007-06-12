/*** MatrixLT.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  07/08/2006
//
//                             BASICS:  07/08/2006
//   print, printI, fillBy, mat2array:  08/08/2006
//                           cholesky:  10/08/2006
//        chinv, chol_inv, chol_solve:  11/08/2006
//           chol_invPD, chol_solvePD:  13/08/2006
//                         add_b2diag:  13/09/2006
//                         choleskyPD:  14/09/2006
//                             add_A2:  18/09/2006
//                             tLxAxL:  21/09/2006
//                          array2mat:  27/09/2006
//                            print4R:  12/10/2006
//                           sqrtDiag:  23/10/2006
//                               LxtL:  24/10/2006
//
// PURPOSE: Functions to manipulate symmetric matrices stored as a lower triangle
//          (in column major order)
//
/* ********************************************************************************* */

#ifndef _MATRIX_LOWER_TRIANGLE_CPP_
#define _MATRIX_LOWER_TRIANGLE_CPP_

#include "MatrixLT.h"

/* ********************************************************************************* */
/* Constructors and destructors                                                      */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */

/***** Nonparametric constructor *****/
template <typename dataType>
MatrixLT<dataType>::MatrixLT()
  : _nrow(0), _length(0), _a(NULL), _atemp(NULL), _diagI(NULL)
{
}


/***** Parametric constructor 1  *****/
// 
// type:  type of the data supplied
//          0: a = lower triangle stored in column major order
//          1: a = full symmetric matrix stored in column major order
//
template <typename dataType>
MatrixLT<dataType>::MatrixLT(const int& nrow, const dataType* a, const int& type)
{
  int i, j;
  dataType *_aP; 

  if (nrow < 0) throw returnR("MatrixLT.cpp: MatrixLT::MatrixLT(nrow, a, type) error", 1);
  _nrow = nrow;
  _length = (_nrow * (_nrow+1))/2;

  if (_length){
    _diagI = (int*) calloc(_nrow, sizeof(int));
    if (!_diagI) throw returnR("Out of memory in MatrixLT.cpp: MatrixLT::MatrixLT(nrow, a, type)", 99);
    for (i = 0; i < _nrow; i++){
      _diagI[i] = (i * (2*_nrow - i + 1))/2;
    }

    _a = (dataType*) calloc(_length, sizeof(dataType));
    _atemp = (dataType*) calloc(_length, sizeof(dataType));
    if (!_a || !_atemp) throw returnR("Out of memory in MatrixLT.cpp: MatrixLT::MatrixLT(nrow, a, type)", 99);
    switch (type){
    case 0:
      for (i = 0; i < _length; i++){
        _a[i] = a[i];
      }
      break;
    case 1:
      _aP = _a;
      for (j = 0; j < _nrow; j++){
        for (i = j; i < _nrow; i++){
          *_aP = a[j*_nrow+i];
          _aP++;
        }
      }
      break;
    default:
      throw returnR("MatrixLT.cpp: MatrixLT::MatrixLT(nrow, a, type) error. Unknown type argument", 1);   
    }
  }
  else{
    _a = NULL;
    _atemp = NULL;
    _diagI = NULL;
  }
}  /** end of the parametric constructor 1 **/


/***** Parametric constructor 2  *****/
//
// Creates a matrix filled by zeros 
//
template <typename dataType>
MatrixLT<dataType>::MatrixLT(const int& nrow)
{
  int i;
 
  if (nrow < 0) throw returnR("MatrixLT.cpp: MatrixLT::MatrixLT(nrow) error", 1);
  _nrow = nrow;
  _length = (_nrow * (_nrow+1))/2;

  if (_length){
    _diagI = (int*) calloc(_nrow, sizeof(int));
    if (!_diagI) throw returnR("Out of memory in MatrixLT.cpp: MatrixLT::MatrixLT(nrow)", 99);
    for (i = 0; i < _nrow; i++){
      _diagI[i] = (i * (2*_nrow - i + 1))/2;
    }

    _a = (dataType*) calloc(_length, sizeof(dataType));
    _atemp = (dataType*) calloc(_length, sizeof(dataType));
    if (!_a || !_atemp) throw returnR("Out of memory in MatrixLT.cpp: MatrixLT::MatrixLT(nrow)", 99);
    for (i = 0; i < _length; i++){
      _a[i] = 0;
    }
  }
  else{
    _a = NULL;
    _atemp = NULL;
    _diagI = NULL;
  }

}  /** end of the parametric constructor 2 **/


/***** Parametric constructor 3  *****/
//
// Create a matrix from a vector as b*b^T
//
template <typename dataType>
MatrixLT<dataType>::MatrixLT(const int& nrow, const dataType* b)
{
  int i, j;
  const dataType *b1, *b2;
  dataType *aP;

  if (nrow < 0) throw returnR("MatrixLT.cpp: MatrixLT::MatrixLT(nrow, b) error", 1);
  _nrow = nrow;
  _length = (_nrow * (_nrow+1))/2;

  if (_length){
    _diagI = (int*) calloc(_nrow, sizeof(int));
    if (!_diagI) throw returnR("Out of memory in MatrixLT.cpp: MatrixLT::MatrixLT(nrow, b)", 99);
    for (i = 0; i < _nrow; i++){
      _diagI[i] = (i * (2*_nrow - i + 1))/2;
    }

    _a = (dataType*) calloc(_length, sizeof(dataType));
    _atemp = (dataType*) calloc(_length, sizeof(dataType));
    if (!_a || !_atemp) throw returnR("Out of memory in MatrixLT.cpp: MatrixLT::MatrixLT(nrow, b)", 99);
    b2 = b;
    aP = _a;
    for (j = 0; j < _nrow; j++){
      b1 = b + j;
      for (i = j; i < _nrow; i++){
        *aP = (*b1) * (*b2);
        aP++;
        b1++;
      }
      b2++;
    }
  }
  else{
    _a = NULL;
    _atemp = NULL;
    _diagI = NULL;
  }

}  /** end of the parametric constructor 3 **/


/***** Copy constructor *****/
template <typename dataType>
MatrixLT<dataType>::MatrixLT(const MatrixLT<dataType>& A)
{
  int i;

  _nrow = A.nrow();
  _length = A.length();

  if (_length){
    _diagI = (int*) calloc(_nrow, sizeof(int));
    if (!_diagI) throw returnR("Out of memory in MatrixLT.cpp: MatrixLT::MatrixLT(A)", 99);
    for (i = 0; i < _nrow; i++){
      _diagI[i] = A.diagI(i);
    }

    _a = (dataType*) calloc(_length, sizeof(dataType));
    _atemp = (dataType*) calloc(_length, sizeof(dataType));
    if (!_a || !_atemp) throw returnR("Out of memory in MatrixLT.cpp: MatrixLT::MatrixLT(A)", 99);
    for (i = 0; i < _length; i++){
      _a[i] = A.a(i);
      _atemp[i] = A._atemp[i];
    }
  }
  else{
    _a = NULL;
    _atemp = NULL;
    _diagI = NULL;
  }
}


/***** Assignment operator *****/
template <typename dataType>
MatrixLT<dataType>&
MatrixLT<dataType>::operator=(const MatrixLT<dataType>& A)
{
  int i;

  if (_length){
    free(_a);
    free(_atemp);
    free(_diagI);
  }

  _nrow = A.nrow();
  _length = A.length();

  if (_length){
    _diagI = (int*) calloc(_nrow, sizeof(int));
    if (!_diagI) throw returnR("Out of memory in MatrixLT.cpp: MatrixLT::operator=", 99);
    for (i = 0; i < _nrow; i++){
      _diagI[i] = A.diagI(i);
    }

    _a = (dataType*) calloc(_length, sizeof(dataType));
    _atemp = (dataType*) calloc(_length, sizeof(dataType));
    if (!_a || !_atemp) throw returnR("Out of memory in MatrixLT.cpp: MatrixLT::operator=", 99);
    for (i = 0; i < _length; i++){
      _a[i] = A.a(i);
      _atemp[i] = A._atemp[i];
    }
  }
  else{
    _a = NULL;
    _atemp = NULL;
    _diagI = NULL;
  }

  return *this;
}


/***** Destructor *****/
template <typename dataType>
MatrixLT<dataType>::~MatrixLT()
{
  if (_length){
    free(_a);
    free(_atemp);
    free(_diagI);
  }
}


/* ********************************************************************************* */
/* Utilities                                                                         */
/*                                                                                   */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */

/***** Check for NaN (works well only for real matrices!!!)                          */
/*       false for Inf, -Inf, NA, NaN                                                */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
bool
MatrixLT<dataType>::is_finite() const
{
  const dataType *a;
  for (int i = 0; i < _length; i++){
    if (!R_finite(*a)) return false;
    a++;
  }

  return true;
}

/***** sqrtDiag: Compute square roots of diagonal elements ****************************************/
/*                                                                                                */
/*                                                                                                */
/*                                                                                                */
/* ---------------------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixLT<dataType>::sqrtDiag(dataType *result) const
{
  static int i;
  static const dataType *aP;
  static dataType *resP;

  aP = _a;
  resP = result;
  for (i = _nrow; i > 0; i--){
    *resP = sqrt(*aP);
    aP += i;
    resP++;
  }
  
  return;
}


/* ********************************************************************************* */
/* fillBy: Fill the matrix by a specific value                                       */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixLT<dataType>::fillBy(const dataType& val)
{
  dataType *a;
  a = _a;
  for (int i = 0; i < _length; i++){
    *a = val;
    a++;
  }

  return;
} 


/* ********************************************************************************* */
/* mat2array: copy matrix to the array (in column major order)                       */
/*                                                                                   */
/* type:  type of the array required                                                 */
/*          0: put only a lower triangle to the array                                */
/*          1: put the full matrix to the array                                      */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixLT<dataType>::mat2array(dataType* a, const int& type) const
{
  int i, j;
  const dataType* _aP = _a;

  switch (type){
  case 0:
    for (i = 0; i < _length; i++){
      a[i] = _a[i];
    }
    break;
  case 1:
    for (j = 0; j < _nrow; j++){
      a[j*_nrow + j] = *_aP;
      _aP++;
      for (i = j+1; i < _nrow; i++){
        a[j*_nrow + i] = *_aP;
        a[i*_nrow + j] = *_aP;
        _aP++;
      }
    }
    break;
  default:
    throw returnR("MatrixLT.cpp: MatrixLT::mat2array(a, type) error. Unknown type argument", 1);   
  }

  return;
}


/* ********************************************************************************* */
/* array2mat: copy array to *this, without checking dimensions!!!                    */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixLT<dataType>::array2mat(const dataType* a)
{
  const dataType *aP = a;
  dataType *tP = _a;
  for (int i = 0; i < _length; i++){
    *tP = *aP;
    tP++;
    aP++;
  }
    
  return;
}


/* ********************************************************************************* */
/* Linear algebra                                                                    */
/*                                                                                   */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
/*                                                                                   */
/*   !!! In most cases, no checks are performed concerning dimensionality            */
/*       of input and output matrices !!!                                            */
/*   !!! In most cases, it is assumed that both input and output matrices            */
/*       are not of _length=0 !!!                                                    */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */

/***** add_b2diag: Add vector to the diagonal of (*this)                         *****/
/*                                                                                   */
/* * (*this) is a matrix                                                             */
/* * b must be a vector of the length = this->nrow()                                 */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixLT<dataType>::add_b2diag(const dataType* b)
{
  dataType *a = _a;
  const dataType *bP = b;
  for (int i = 0; i < _nrow; i++){
    *a += (*bP);
    a += _nrow - i;
    bP++;
  }

  return;
}


/***** add_A2: Add matrix to (*this)                                             *****/
/*                                                                                   */
/* * (*this) is a matrix                                                             */
/* * A must be an array of the same length as this->a()                              */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixLT<dataType>::add_A2(const dataType* A)
{
  dataType *a = _a;
  const dataType *AP = A;
  for (int i = 0; i < _length; i++){
    *a += (*AP);
    a++;
    AP++;
  }

  return;
}


/* ********************************************************************************* */
/* Matrix decompositions                                                             */
/*                                                                                   */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */

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
MatrixLT<dataType>::cholesky(const int &backup)
{
  int rank = cholesky0(_a, _atemp, _nrow, _diagI, backup);
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
MatrixLT<dataType>::chinv(const int& onlyCholInv)
{
  chinv0(_a, _atemp, _nrow, _diagI, onlyCholInv);
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
MatrixLT<dataType>::chol_inv(const int &what, dataType *Li)
{
  int rank = chol_inv0(_a, _atemp, _nrow, _diagI, what, Li);
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
MatrixLT<dataType>::chol_solve(dataType *c, dataType *x, const int &add_u, const int &solve_u, dataType *u)
{
  int rank = chol_solve0(_a, _atemp, _nrow, _diagI, c, x, add_u, solve_u, u);
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
MatrixLT<dataType>::choleskyPD(int &Attempt, const int &nAttempt, const dataType &eps)
{
  int rank = choleskyPD0(_a, _atemp, _nrow, _diagI, Attempt, nAttempt, eps);
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
MatrixLT<dataType>::chol_invPD(int &Attempt,        const int &what,      dataType *Li, 
                               const int &nAttempt, const dataType &eps)
{
  int rank = chol_invPD0(_a, _atemp, _nrow, _diagI, Attempt, what, Li, nAttempt, eps);
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
MatrixLT<dataType>::chol_solvePD(int &Attempt,        dataType *c,          dataType *x, 
                                 const int &nAttempt, const dataType &eps,
                                 const int &add_u,    const int &solve_u,   dataType *u)
{
  int rank = chol_solvePD0(_a, _atemp, _nrow, _diagI, Attempt, c, x, nAttempt, eps, add_u, solve_u, u);
  return rank;
}


/* ********************************************************************************* */
/* print: print symmetric matrix stored as a lower triangle                          */
/*             in column major order                                                 */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixLT<dataType>::print(const int &info) const
{
  int i, j;
  dataType x;
  if (info) Rprintf("Symmetric matrix (nrow=%d):\n", _nrow);
  if (_nrow > 0){
    for (i = 0; i < _nrow; i++){
      for (j = 0; j < _nrow; j++){
        x = a(i, j);
        Rprintf("%5g  ", (fabs(x) < _AK_ZERO ? 0 : x));      
      }
      Rprintf("\n");
    }
  }

  return;
}

template <typename dataType>
void
MatrixLT<dataType>::print4R() const
{
  int i, j;
  dataType x;
  if (_nrow > 0){
    Rprintf("matrix(c(");
    Rprintf("%5g", (fabs(*_a) < _AK_ZERO ? 0 : *_a));      
    for (i = 1; i < _nrow; i++){
      x = a(i, 0);
      Rprintf(", %5g", (fabs(x) < _AK_ZERO ? 0 : x));
    }  

    for (j = 1; j < _nrow; j++){
      for (i = 0; i < _nrow; i++){
        x = a(i, j);
        Rprintf(", %5g", (fabs(x) < _AK_ZERO ? 0 : x));      
      }
    }

    Rprintf("), nrow=%d, ncol=%d)\n", _nrow, _nrow);
  }

  return;
}


template <typename dataType>
void
MatrixLT<dataType>::printI(const int &info) const
{
  int i, j;
  if (info) Rprintf("Symmetric matrix (nrow=%d):\n", _nrow);
  if (_nrow > 0){
    for (i = 0; i < _nrow; i++){
      for (j = 0; j < _nrow; j++){
        Rprintf("%5d  ", a(i, j));
      }
      Rprintf("\n");
    }
  }

  return;
}


/* ********************************************************************************* */
/* Linear algebra, not part of MatrixLT class                                        */
/*                                                                                   */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */


/***** tLxAxL: Compute t(L)*A*L, where L is lower triangular matrix and A a symmetric matrix ******/
/*                                                                                                */
/*  tLAL: RESULT, array of length LT(p)                                                           */
/*   tLA: working array of length p*p                                                             */
/*        on OUTPUT it contains t(L)*A stored in COLUMN major order                               */
/*     L: lower triangular matrix in an array of length LT(p)                                     */
/*     A: symmetric matrix whose lower triangle is stored in an array of length LT(p)             */
/* diagI: indeces of diagonal elements, diagI[i] = (i*(2*p-i+1))/2                                */
/*     p: number of rows and columns of all matrices                                              */
/*                                                                                                */
/* ---------------------------------------------------------------------------------------------- */
template <typename dataType>
void
tLxAxL(dataType *tLAL, dataType *tLA, const dataType *L,  const dataType *A,  const int *diagI, const int *p)
{
  static int i, j, k;
  static const dataType *LP, *AP;
  static dataType *tLAP, *tLAP2, *tLALP;

  /*** Compute t(L)*A ***/
  tLAP = tLA;  
  for (j = 0; j < *p; j++){
    LP = L;
    for (i = 0; i <= j; i++){
      *tLAP = (*LP) * A[diagI[i] + j - i];
      LP++;
      for (k = i + 1; k <= j; k++){
        *tLAP += (*LP) * A[diagI[k] + j - k];
        LP++;
      }
      for (k = j + 1; k < *p; k++){
        *tLAP += (*LP) * A[diagI[j] + k - j];
        LP++;
      }
      tLAP++;
    }

    for (i = j + 1; i < *p; i++){
      *tLAP = (*LP) * A[diagI[j] + i - j];
      LP++;
      for (k = i + 1; k < *p; k++){
        *tLAP += (*LP) * A[diagI[j] + k - j];
        LP++;
      }
      tLAP++;
    }
  }


  /*** Compute t(L)*A*L ***/
  tLALP = tLAL;
  LP = L - (*p) - 1;
  tLAP = tLA;
  for (j = 0; j < *p; j++){
    LP += (*p) - j + 1;          /* go to L[j,j]                    */   
    tLAP += j;                   /* go to a[j,j] of tLA             */
    for (i = j; i < *p; i++){
      tLAP2 = tLAP;
      *tLALP = (*tLAP) * (*LP);
      LP++;
      tLAP2 += *p;
      for (k = j + 1; k < *p; k++){
        *tLALP += (*tLAP2) * (*LP);
        LP++;
        tLAP2 += *p;
      }
      LP -= (*p) - j;           /* go back to L[j,j]          */
      tLAP++;                   /* go to a[i+1,j] of tLA      */
      tLALP++;
    }
  }

  return;
}


/***** LxtL: Compute L*t(L), where L is lower triangular matrix ***********************************/
/*                                                                                                */
/*   LtL: RESULT, array of length LT(p)                                                           */
/*     L: lower triangular matrix in an array of length LT(p)                                     */
/*     p: number of rows and columns of all matrices                                              */
/*                                                                                                */
/* Idea of the algorithm, see p. 46 of red notes.                                                 */
/*                                                                                                */
/* ---------------------------------------------------------------------------------------------- */
template <typename dataType>
void
LxtL(dataType *LtL, const dataType *L, const int *p)
{
  static int i, j, i2;
  static dataType *LtLP, *startLtLP;
  static const dataType *LP, *LP2;

  /*** Column 0 of LtL and initialization of the remaining columns of LtL ***/
  startLtLP = LtL;
  LtLP = startLtLP;
  LP = L; 
  for (i = 0; i < *p; i++){
    LP2 = LP;
    for (i2 = i; i2 < *p; i2++){
      *LtLP = (*LP) * (*LP2);
      LP2++;
      LtLP++;
    }
    LP++;
  }
  startLtLP += *p;

  /*** Column 1,...,p-1 of LtL ***/
  for (j = 1; j < *p; j++){
    LtLP = startLtLP;
    for (i = j; i < *p; i++){
      LP2 = LP;
      for (i2 = i; i2 < *p; i2++){
        *LtLP += (*LP) * (*LP2);
        LP2++;
        LtLP++;
      }
      LP++;
    }
    startLtLP += *p - j;
  }
  
  return;
}


#endif
