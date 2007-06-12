/*** MatrixRect.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  08/08/2006
//
//                                                  BASICS:  08/08/2006
//                        print, printI, fillBy, mat2array:  08/08/2006
//                                          Ab, bA, bA2row:  09/08/2006
//                                        vec2array, add_A:  13/08/2006
//                                             copy, rNorm:  13/09/2006
//                                         xTx, subtract_A:  14/09/2006
//                               sum, anyNonNeg, BAcolProd:  17/09/2006
//               bA2, BAcolProd2, Ab2, add_A2, subtract_A2:  18/09/2006
//                                  max, min, b_plus_rowsA:  19/09/2006
//                                                    mean:  26/09/2006
//                  square (2 overloaded versions), divide:  21/03/2007
//                                                 
//
// PURPOSE: Functions to manipulate general matrices stored in an array
//          (in column major order)
//
/* ********************************************************************************* */

#ifndef _MATRIX_RECTANGLE_CPP_
#define _MATRIX_RECTANGLE_CPP_

#include "MatrixRect.h"

/* ********************************************************************************* */
/* Constructors and destructors                                                      */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */

/***** Nonparametric constructor *****/
template <typename dataType>
MatrixRect<dataType>::MatrixRect()
  : _nrow(0), _ncol(0), _length(0), _a(NULL)
{
}


/***** Parametric constructor 1  *****/
// 
// type:  type of the data supplied
//          0: a = matrix in column major order
//          1: a = matrix in row major order
//
template <typename dataType>
MatrixRect<dataType>::MatrixRect(const int& nrow, const int& ncol, const dataType* a, const int& type)
{
  int i, j;
  dataType *_aP;
 
  if (nrow < 0 || ncol < 0) throw returnR("MatrixRect.cpp: MatrixRect::MatrixRect(nrow, ncol, a, type) error", 1);
  _nrow = nrow;
  _ncol = ncol;
  _length = _nrow * _ncol;
  if (!_length){
    _nrow = 0;
    _ncol = 0;
  }

  if (_length){
    _a = (dataType*) calloc(_length, sizeof(dataType));
    if (!_a) throw returnR("Out of memory in MatrixRect.cpp: MatrixRect::MatrixRect(nrow, ncol, a, type)", 99);
    switch (type){
    case 0:
      for (i = 0; i < _length; i++){
        _a[i] = a[i];
      }
      break;
    case 1:
      _aP = _a;
      for (j = 0; j < _ncol; j++){
        for (i = 0; i < _nrow; i++){
          *_aP = a[i*_ncol+j];
          _aP++;
        }
      }
      break;
    default:
      throw returnR("MatrixRect.cpp: MatrixRect::MatrixRect(nrow, ncol, a, type) error. Unknown type argument", 1);   
    }
  }
  else{
    _a = NULL;
  }
}  /** end of the parametric constructor 1 **/


/***** Parametric constructor 2  *****/
//
// Creates a matrix filled by zeros 
//
template <typename dataType>
MatrixRect<dataType>::MatrixRect(const int& nrow, const int& ncol)
{
  int i;
 
  if (nrow < 0 || ncol < 0) throw returnR("MatrixRect.cpp: MatrixRect::MatrixRect(nrow, ncol) error", 1);
  _nrow = nrow;
  _ncol = ncol;
  _length = _nrow * _ncol;
  if (!_length){
    _nrow = 0;
    _ncol = 0;
  }

  if (_length){
    _a = (dataType*) calloc(_length, sizeof(dataType));
    if (!_a) throw returnR("Out of memory in MatrixRect.cpp: MatrixRect::MatrixRect(nrow, ncol)", 99);
    for (i = 0; i < _length; i++){
      _a[i] = 0;
    }
  }
  else{
    _a = NULL;
  }

}  /** end of the parametric constructor 2 **/


/***** Parametric constructor 3  *****/
//
// Create a matrix from 2 vectors as 
//    either 1) type=0: b * c^T (DEFAULT CHOICE)
//                        INPUT: b[nrow]
//                               c[ncol]
//    or     2) type=1: b * c (elementwise product)
//                        INPUT: b[nrow*ncol]
//                               c[nrow*ncol]
//
template <typename dataType>
MatrixRect<dataType>::MatrixRect(const int& nrow, const int& ncol, const dataType* b, const dataType* c, const int &type)
{
  int i, j;
  const dataType *b1, *c1;
  dataType *aP;

  if (nrow < 0 || ncol < 0) throw returnR("MatrixRect.cpp: MatrixRect::MatrixRect(nrow, ncol, b, c, type) error", 1);
  _nrow = nrow;
  _ncol = ncol;
  _length = _nrow * _ncol;
  if (!_length){
    _nrow = 0;
    _ncol = 0;
  }

  if (_length){
    _a = (dataType*) calloc(_length, sizeof(dataType));
    if (!_a) throw returnR("Out of memory in MatrixRect.cpp: MatrixRect::MatrixRect(nrow, ncol, b, c, type)", 99);

    switch (type){
    case 0:
      c1 = c;
      aP = _a;
      for (j = 0; j < _ncol; j++){
        b1 = b;
        for (i = 0; i < _nrow; i++){
          *aP = (*b1) * (*c1);
          aP++;
          b1++;
        }
        c1++;
      }
      break;

    case 1:
      b1 = b;
      c1 = c;
      aP = _a;
      for (j = 0; j < _ncol; j++){
        for (i = 0; i < _nrow; i++){
          *aP = (*b1) * (*c1);
          b1++;
          c1++;
          aP++;
        }
      }
      break;

    default:
      throw returnR("MatrixRect.cpp: MatrixRect::MatrixRect(nrow, ncol, b, c, type) error", 1);
    }
  }
  else{
    _a = NULL;
  }

}  /** end of the parametric constructor 3 **/


/***** Copy constructor *****/
template <typename dataType>
MatrixRect<dataType>::MatrixRect(const MatrixRect<dataType>& A)
{
  int i;

  _nrow = A.nrow();
  _ncol = A.ncol();
  _length = A.length();

  if (_length){
    _a = (dataType*) calloc(_length, sizeof(dataType));
    if (!_a) throw returnR("Out of memory in MatrixRect.cpp: MatrixRect::MatrixRect(A)", 99);
    for (i = 0; i < _length; i++){
      _a[i] = A.a(i);
    }
  }
  else{
    _a = NULL;
  }
}


/***** Assignment operator *****/
template <typename dataType>
MatrixRect<dataType>&
MatrixRect<dataType>::operator=(const MatrixRect<dataType>& A)
{
  int i;

  if (_length){
    free(_a);
  }

  _nrow = A.nrow();
  _ncol = A.ncol();
  _length = A.length();

  if (_length){
    _a = (dataType*) calloc(_length, sizeof(dataType));
    if (!_a) throw returnR("Out of memory in MatrixRect.cpp: MatrixRect::operator=", 99);
    for (i = 0; i < _length; i++){
      _a[i] = A.a(i);
    }
  }
  else{
    _a = NULL;
  }

  return *this;
}


/***** Destructor *****/
template <typename dataType>
MatrixRect<dataType>::~MatrixRect()
{
  if (_length){
    free(_a);
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
MatrixRect<dataType>::is_finite() const
{
  const dataType *a;
  for (int i = 0; i < _length; i++){
    if (!R_finite(*a)) return false;
    a++;
  }

  return true;
}

/***** Copy (part) of 1 matrix to the (part) of *this                                */
/*       * operates directly on the arrays                                           */
/*       * no checks concerning dimensionality performed!!!                          */
/*                                                                                   */
/*       length:  how many elements should be taken from b                           */
/*        bfrom:  where to start at b (DEFAULT = 0)                                  */
/*     thisfrom:  where to start at (*this) (DEFAULT = 0)                            */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::copy(const MatrixRect<dataType> *b,  const int &length,  const int &bfrom,  const int &thisfrom)
{
  int i;

  const dataType *bP;
  dataType *aP;
  bP = b->aconst() + bfrom;
  aP = _a + thisfrom;
  for (i = 0; i < length; i++){
    *aP = *bP;
    bP++;
    aP++;
  }

  return;
}

/***** Fill in the vector by random numbers generated from N(mean, sd) *****/
/*                                                                         */
/* ----------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::rNorm(const double& mean,  const double& sd)
{
  int i;
  dataType *aP = _a;
  for (i = 0; i < _length; i++){
    *aP = rnorm(mean, sd);
    aP++;
  }

  return;
}


/* ********************************************************************************* */
/* Matrix operators                                                                  */
/*                                                                                   */
/*   !!! In most cases, no checks are performed concerning dimensionality            */
/*       of input and output matrices !!!                                            */
/*   !!! In most cases, it is assumed that both input and output matrices            */
/*       are not of _length=0 !!!                                                    */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */


/***** xTx: scalar product of two vectors (t(x)*x) *****/
/*                                                                                    */
/* * (*this) can be both row or column vector                                         */
/*                                                                                    */
/* ---------------------------------------------------------------------------------- */
template <typename dataType>
dataType
MatrixRect<dataType>::xTx() const
{
  static int i;
  const dataType *xP = _a;
  double RESULT;

  RESULT = (*xP)*(*xP);
  for (i = 1; i < _length; i++){
    xP++;
    RESULT += (*xP)*(*xP);
  }
  
  return RESULT;
}


/***** Ab: Product of a general (Rect) matrix and a vector *****/
/*                                                                                   */
/* * (*this) can be both row or column vector                                        */
/* * length of (*this) must be equal to A->nrow()                                    */
/* * b can be both row or column vector                                              */
/* * from b we use b[bstart],...,b[bstart+A->ncol()-1]                               */
/* * default for bstart is 0                                                         */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::Ab(const MatrixRect<dataType>* A, const MatrixRect<dataType>* b, const int& bstart)
{
  static int i, j;
  const dataType *AP, *bP;
  dataType *aP;

  aP = _a;
  for (i = 0; i < A->nrow(); i++){    
    AP = A->aconst() + i;
    bP = b->aconst()+bstart;
    *aP = 0;
    for (j = 0; j < A->ncol(); j++){
      *aP += (*AP) * (*bP);
      AP += A->nrow();
      bP++;
    }
    aP++;
  }

  return;
}


/***** Ab2: Product of a general (Rect) matrix and a vector, vector given as an array *****/
/*                                                                                        */
/*   b: array of length A->ncol()                                                         */
/*                                                                                        */
/* * (*this) can be both row or column vector                                             */
/* * length of (*this) must be equal to A->nrow()                                         */
/* * b can be both row or column vector                                                   */
/*                                                                                        */
/* -------------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::Ab2(const MatrixRect<dataType>* A, const dataType* b)
{
  static int i, j;
  const dataType *AP, *bP;
  dataType *aP;

  aP = _a;
  for (i = 0; i < A->nrow(); i++){    
    AP = A->aconst() + i;
    bP = b;
    *aP = 0;
    for (j = 0; j < A->ncol(); j++){
      *aP += (*AP) * (*bP);
      AP += A->nrow();
      bP++;
    }
    aP++;
  }

  return;
}


/***** bA: Product of a vector and a general (Rect) matrix                               *****/
/*                                                                                           */
/*     TYPICAL USAGE: compute linear predictor                                               */
/*             b: regression coefficients                                                    */ 
/*             A: design matrix, each column = covariates for a specific observation         */
/*                                                                                           */
/*                    (*this): vector of length n                                            */
/*   b[bstart,...,bstart+p-1]: vector of length p                                            */
/*                          A: matrix p x n                                                  */
/*                                                                                           */
/*     DEFAULTS: bstart = 0                                                                  */
/*                                                                                           */
/*     RESULT:                                                                               */
/*         (*this)[i] = t(b) * A[,i]                                                         */
/*                                                                                           */
/* ----------------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::bA(const MatrixRect<dataType>* b, const MatrixRect<dataType>* A, const int& bstart)
{
  static int i, j;
  const dataType *AP, *bP;
  dataType *aP;

  aP = _a;
  AP = A->aconst();
  for (i = 0; i < A->ncol(); i++){    
    bP = b->aconst()+bstart;
    *aP = (*AP) * (*bP);
    AP++;
    bP++;
    for (j = 1; j < A->nrow(); j++){
      *aP += (*AP) * (*bP);
      AP++;
      bP++;
    }
    aP++;
  }

  return;
}


/***** bA2: Product of a vector b and a general (Rect) matrix A                          *****/
/*          each column of A is used C-times                                                 */
/*                                                                                           */
/*     TYPICAL USAGE: compute linear predictor, part not proportional w.r.t. odds            */
/*                    in the cumulative logit model                                          */
/*             b: regression coefficients                                                    */ 
/*             A: design matrix, each column = covariates for a specific observation         */
/*                                                                                           */
/*                    (*this): matrix C x n                                                  */
/* b[bstart,...,bstart+C*q-1]: vector of length C*q                                          */
/*                          A: matrix q x n                                                  */
/*                                                                                           */
/*     DEFAULTS: bstart = 0                                                                  */
/*                                                                                           */
/*     RESULT:                                                                               */
/*         (*this)[c,i] = t(b[c*q,...,(c+1)*q-1]) * A[,i]                                    */
/*              that is, each column of A is used C-times                                    */
/*                                                                                           */
/* ----------------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::bA2(const MatrixRect<dataType>* b, const MatrixRect<dataType>* A, const int& bstart)
{
  static int i, j, c;
  const dataType *AP, *bP;
  dataType *aP;

  /*  const int C = _nrow;  */

  aP = _a;
  AP = A->aconst();
  for (i = 0; i < A->ncol(); i++){    
    bP = b->aconst() + bstart;

    /* c = 0 */
    *aP = (*AP) * (*bP);
    AP++;
    bP++;
    for (j = 1; j < A->nrow(); j++){
      *aP += (*AP) * (*bP);
      AP++;
      bP++;
    }
    aP++;
    
    /* c = 1, ... */
    for (c = 1; c < _nrow; c++){
      AP -= A->nrow();
      *aP = (*AP) * (*bP);
      AP++;
      bP++;
      for (j = 1; j < A->nrow(); j++){
        *aP += (*AP) * (*bP);
        AP++;
        bP++;
      }
      aP++;      
    }
  }

  return;
}


/***** bA2row: Product of a vector and a general (Rect) matrix                   *****/
/*                                                                                   */
/*     TYPICAL USAGE: originally used in a loop instead of bA2 function              */ 
/* * (*this) is a matrix                                                             */
/* * result is filled in the 'row'th row of (*this)                                  */
/* * this->ncol() must be equal to A->ncol()                                         */
/* * b can be both row or column vector                                              */
/* * from b we use b[bstart],...,b[bstart+A->nrow()-1]                               */
/* * default for bstart is 0                                                         */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
//template <typename dataType>
//void
//MatrixRect<dataType>::bA2row(const int& row, const MatrixRect<dataType>* b, const MatrixRect<dataType>* A, const int& bstart)
//{
//  static int i, j;
//  const dataType *AP, *bP;
//  dataType *aP;
//
//  aP = _a + row;
//  AP = A->aconst();
//  for (j = 0; j < A->ncol(); j++){    
//    bP = b->aconst() + bstart;
//    *aP = 0;
//    for (i = 0; i < A->nrow(); i++){
//      *aP += (*AP) * (*bP);
//      AP++;
//      bP++;
//    }
//    aP += _nrow;
//  }
//
//  return;
//}


/***** BAcolProd: Scalar product of (recycled) columns of B and A                        *****/
/*                                                                                           */
/*     TYPICAL USAGE: compute linear predictor based on cluster-specific random effects      */
/*             B: random effects, each column = r.eff. for 1 cluster                         */ 
/*            ni: numbers of observations per each cluster                                   */
/*             A: design matrix for random effects, each column = covariates                 */
/*                for a specific observation                                                 */
/*                                                                                           */
/*                     (*this): vector of length n                                           */
/*  B[bstart,...,bstart+p-1, ]: matrix p x N                                                 */
/*                          ni: vector of length N                                           */
/*                              sum(ni) = n                                                  */
/*                           A: matrix p x n                                                 */
/*                                                                                           */
/* * from each column of B we take the elements bstart, ..., bstart + A->nrow() - 1          */
/* * the i-th column of B is used ni[i]-times                                                */
/*                                                                                           */
/* DEFAULTS: bstart = 0                                                                      */
/*                                                                                           */
/* RESULT: this[i] = scalar product of B[,istar] and A[,i], i=0,...,n-1                      */
/*                 = t(B[,istar]) * A[,i]                                                    */
/*                   where istar is determined from ni                                       */
/*                                                                                           */
/* ----------------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::BAcolProd(const MatrixRect<dataType> *B, const MatrixRect<int> *ni, const MatrixRect<dataType> *A, const int &bstart)
{
  static int i, k, l;
  const dataType *AP, *BP;
  const int *niP;
  dataType *aP;

  const int skip = B->nrow() - A->nrow();

  AP = A->aconst();
  BP = B->aconst() + bstart;
  niP = ni->aconst();
  aP = _a;
  for (k = 0; k < ni->length(); k++){

    /* l = 0 */
    *aP = (*AP)*(*BP);
    AP++;
    BP++;
    for (i = 1; i < A->nrow(); i++){
      *aP += (*AP)*(*BP);
      AP++;
      BP++;
    }
    aP++;

    /* l = 1,... */
    for (l = 1; l < *niP; l++){
      BP -= A->nrow();
      *aP = (*AP)*(*BP);
      AP++;
      BP++;
      for (i = 1; i < A->nrow(); i++){
        *aP += (*AP)*(*BP);
        AP++;
        BP++;
      }
      aP++;
    }

    BP += skip;
    niP++;
  }

  return;
}


/***** BAcolProd2: Scalar product of columns of B and A                                  *****/
/*                                                                                           */
/*     TYPICAL USAGE: compute linear predictor based on cluster-specific random effects,     */
/*                    part not proportional w.r.t. odds in the cumulative logit model        */
/*                                                                                           */
/*             B: random effects, each column = r.eff. for 1 cluster                         */ 
/*            ni: numbers of observations per each cluster                                   */
/*             A: design matrix for random effects, each column = covariates                 */
/*                for a specific observation                                                 */
/*                                                                                           */
/*                       (*this): matrix C x n                                               */
/*  B[bstart,...,bstart+C*q-1, ]: matrix C*q x N                                             */
/*                            ni: vector of length N                                         */
/*                                sum(ni) = n                                                */
/*                             A: matrix q x n                                               */
/*                                                                                           */
/* * from each column of B we take the elements bstart, ..., bstart + A->nrow() - 1          */
/* * the i-th column of B is used ni[i]-times                                                */
/*                                                                                           */
/* DEFAULTS: bstart = 0                                                                      */
/*                                                                                           */
/* RESULT: this[c,i] = scalar product of B[c*q,...,(c+1)*q-1,istar] and A[,i], i=0,...,n-1   */
/*                   = t(B[,istar]) * A[,i]                                                  */
/*                   where istar is determined from ni                                       */
/*                                                                                           */
/* ----------------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::BAcolProd2(const MatrixRect<dataType> *B, const MatrixRect<int> *ni, const MatrixRect<dataType> *A, const int &bstart)
{
  static int i, k, l, c;
  const dataType *AP, *BP;
  const int *niP;
  dataType *aP;

  /*  const int C = _nrow;  */
  const int back = _nrow * A->nrow();
  const int skip = B->nrow() - back;

  AP = A->aconst();
  BP = B->aconst() + bstart;
  niP = ni->aconst();
  aP = _a;
  for (k = 0; k < ni->length(); k++){

    /* l = 0 */
      /* c = 0 */
    *aP = (*AP)*(*BP);
    AP++;
    BP++;
    for (i = 1; i < A->nrow(); i++){
      *aP += (*AP)*(*BP);
      AP++;
      BP++;
    }
    aP++;

    /* c = 1,... */
    for (c = 1; c < _nrow; c++){
      AP -= A->nrow();
      *aP = (*AP)*(*BP);
      AP++;
      BP++;
      for (i = 1; i < A->nrow(); i++){
        *aP += (*AP)*(*BP);
        AP++;
        BP++;
      }
      aP++;
    }

    /* l = 1,... */    
    for (l = 1; l < *niP; l++){      
      BP -= back;
      
      /* c = 0 */
        *aP = (*AP)*(*BP);
        AP++;
        BP++;
        for (i = 1; i < A->nrow(); i++){
          *aP += (*AP)*(*BP);
          AP++;
          BP++;
        }
        aP++;

      /* c = 1, ... */
      for (c = 1; c < _nrow; c++){
        AP -= A->nrow();
        *aP = (*AP)*(*BP);
        AP++;
        BP++;
        for (i = 1; i < A->nrow(); i++){
          *aP += (*AP)*(*BP);
          AP++;
          BP++;
        }
        aP++;
      }
    }

    BP += skip;
    niP++;
  }

  return;
}


/***** add_A: Add matrix to (*this)                                              *****/
/*                                                                                   */
/* * (*this) is a matrix                                                             */
/* * A must be a matrix of the same dimension as (*this)                             */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::add_A(const MatrixRect<dataType>* A)
{
  static int i;
  dataType *a = _a;
  const dataType *a2 = A->aconst();
  for (i = 0; i < _length; i++){
    *a += (*a2);
    a++;
    a2++;
  }

  return;
}


/***** add_A2: Add matrix stored in an array to (*this)                          *****/
/*                                                                                   */
/* * (*this) is a matrix                                                             */
/* * A must be an array of the same length as (*this)                                */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::add_A2(const dataType* A)
{
  static int i;
  dataType *a = _a;
  const dataType *a2 = A;
  for (i = 0; i < _length; i++){
    *a += (*a2);
    a++;
    a2++;
  }

  return;
}


/***** subtract_A: Subtract matrix from (*this)                                  *****/
/*                                                                                   */
/* * (*this) is a matrix                                                             */
/* * A must be a matrix of the same dimension as (*this)                             */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::subtract_A(const MatrixRect<dataType>* A)
{
  static int i;
  dataType *a = _a;
  const dataType *a2 = A->aconst();
  for (i = 0; i < _length; i++){
    *a -= (*a2);
    a++;
    a2++;
  }

  return;
}


/***** subtract_A2: Subtract matrix given as an array from (*this)               *****/
/*                                                                                   */
/* * (*this) is a matrix                                                             */
/* * A must be an array of the same length as (*this)                                */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::subtract_A2(const dataType* A)
{
  static int i;
  dataType *a = _a;
  const dataType *a2 = A;
  for (i = 0; i < _length; i++){
    *a -= (*a2);
    a++;
    a2++;
  }

  return;
}


/***** add_const: Add constant to (the part of) (*this)                          *****/
/*                                                                                   */
/* * (*this) is a matrix                                                             */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::add_const(const dataType &c, const int &start, const int &stop)
{
  static int i;
  dataType *a = _a + start;
  for (i = 0; i <= stop - start; i++){
    *a += c;
    a++;
  }

  return;
}

template <typename dataType>
void
MatrixRect<dataType>::add_const(const dataType &c)
{
  static int i;
  dataType *a = _a;
  for (i = 0; i < _length; i++){
    *a += c;
    a++;
  }

  return;
}


/***** b_plus_rowsA: Add vector to the rows of A and store it in (*this)         *****/
/*                                                                                   */
/* (*this):  matrix p x n                                                            */
/* b:        vector of length n                                                      */
/* A:        matrix p x n                                                            */
/*                                                                                   */
/*   OUTPUT:  (*this)[i,] = b + A[i,]                                                */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::b_plus_rowsA(const dataType *b,  const MatrixRect<dataType> *A)
{
  static int i, j;
  dataType *aP = _a;
  const dataType *bP = b;
  const dataType *AP = A->aconst();
  for (j = 0; j < A->ncol(); j++){
    for (i = 0; i < A->nrow(); i++){
      *aP = *bP + (*AP);
      aP++;
      AP++;
    }
    bP++;
  }
  
  return;
}


/***** square: Square all elements of (*this)                                    *****/
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::square()
{
  static int i;
  dataType *aP = _a;
  for (i = 0; i < _length; i++){
    *aP = (*aP) * (*aP);
    aP++;
  }
  
  return;
}


/***** square: Assign to (*this) squares of b                                    *****/
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::square(const dataType *b)
{
  static int i;
  dataType *aP       = _a;
  const dataType *bP = b;
  for (i = 0; i < _length; i++){
    *aP = (*bP) * (*bP);
    aP++;
    bP++;
  }
  
  return;
}


/***** square: Divide all elements of (*this) by b                               *****/
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::divide(const dataType &b)
{
  static int i;
  dataType *aP = _a;
  for (i = 0; i < _length; i++){
    *aP /= b; 
    aP++;
  }
  
  return;
}


/***** sum: Compute the sum of elements of (*this)                               *****/
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
dataType
MatrixRect<dataType>::sum() const
{
  static int i;
  const dataType *a = _a;
  dataType RET = 0;
  for (i = 0; i < _length; i++){
    RET += *a;
    a++;
  }

  return RET;
}


/***** max: Find the maximum of elements of (*this)                              *****/
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
dataType
MatrixRect<dataType>::max() const
{
  static int i;
  const dataType *a = _a;
  dataType RET = *a;
  for (i = 1; i < _length; i++){
    a++;
    if (*a > RET) RET = *a;
  }

  return RET;
}


/***** min: Find the minimum of elements of (*this)                              *****/
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
dataType
MatrixRect<dataType>::min() const
{
  static int i;
  const dataType *a = _a;
  dataType RET = *a;
  for (i = 1; i < _length; i++){
    a++;
    if (*a < RET) RET = *a;
  }

  return RET;
}


/***** mean: Compute the mean over rows or columns of the matrix (*this)         *****/
/*                                                                                   */
/*  margin:  1 ... means for each row of (*this) are computed                        */
/*           Mean must have the length this->ncol()                                  */
/*                                                                                   */
/*           2 ... means for each column of (*this) are computed                     */
/*           Mean must have the length this->nrow()                                  */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::mean(dataType *Mean,  const int &margin) const
{
  static int i, j;
  dataType *MP = Mean;
  const dataType *aP;

  switch (margin){
  case 1:              /*** mean over rows ***/
    for (i = 0; i < _nrow; i++){
      *MP = 0;
      aP = _a + i;
      for (j = 0; j < _ncol; j++){
        *MP += *aP;
        aP += _nrow;
      }
      *MP /= _ncol;
      MP++;
    }
    break;

  case 2:              /*** mean over columns ***/
    aP = _a;
    for (j = 0; j < _ncol; j++){
      *MP = 0;
      for (i = 0; i < _nrow; i++){
        *MP += *aP;
        aP++;
      }
      *MP /= _nrow;
      MP++;
    }
    break;

  default:
    throw returnR("MatrixRect.cpp: MatrixRect::mean(Mean, margin) error. Incorrect margin argument", 1);    
  }

  return;
}


/***** anyNonNeg: Check whether any element of (*this) is non-negative           *****/
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
bool
MatrixRect<dataType>::anyNonNeg() const
{
  static int i;
  const dataType *a = _a;
  for (i = 0; i < _length; i++){
    if (*a <= 0) return true;
    a++;
  }

  return false;
}


/* ********************************************************************************* */
/* fillBy: Fill the matrix by a specific value                                       */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::fillBy(const dataType& val)
{
  static int i;
  dataType *a;
  a = _a;
  for (i = 0; i < _length; i++){
    *a = val;
    a++;
  }

  return;
} 


/* ********************************************************************************* */
/* mat2array: copy matrix to the array                                               */
/*                                                                                   */
/* type:  major order of the filling                                                 */
/*          0: column major order (default)                                          */
/*          1: row major order                                                       */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::mat2array(dataType* a, const int& type) const
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
    _aP = _a;
    for (j = 0; j < _ncol; j++){
      for (i = 0; i < _nrow; i++){
         a[i*_ncol+j] = *_aP;
        _aP++;
      }
    }
    break;
  default:
    throw returnR("MatrixRect.cpp: MatrixRect::mat2array(a, type) error. Unknown type argument", 1);   
  }

  return;
}

/* ********************************************************************************* */
/* vec2array: copy (part of the) vector to the array                                 */
/*                                                                                   */
/*       a: array of the length stop - start + 1                                     */
/*   start: first element of this->_a to be copied to a                              */
/*    stop: last element of this->_a to be copied to a                               */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::vec2array(dataType *a, const int &start, const int &stop) const
{
  int i;

  if (start < 0 || start >= _length || stop < 0 || stop >= _length || stop < start)
    throw returnR("MatrixRect.cpp: MatrixRect::vec2array(a, start, stop) error. Incorrect start or stop argument", 1);

  dataType *aa = a;
  const dataType *_aP = _a + start;
  for (i = start; i <= stop; i++){
    *aa = *_aP;
    aa++;
    _aP++;
  }

  return;
}


/* ********************************************************************************* */
/* print: print a matrix stored in an array                                          */
/*        in column major order                                                      */
/*                                                                                   */
/*    info: 0/1, default is 1                                                        */
/*  colmax: if negative, colmax is set to _ncol, default is -1                       */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
MatrixRect<dataType>::print(const int &info,  const int &colmax) const
{
  int i, j, cmax;
  dataType x;

  if (info) Rprintf("Matrix (nrow=%d, ncol=%d):\n", _nrow, _ncol);
  if (_length > 0){
    if (_nrow==1 || _ncol==1){
      cmax = (colmax < 0 ? _length : (colmax <= _length ? colmax : _length));
      for (i = 0; i < cmax; i++){
        x = _a[i];
	Rprintf("%5g  ", (fabs(x) < _AK_ZERO ? 0 : x)); 
      }
      Rprintf("\n");
    }
    else{
      cmax = (colmax < 0 ? _ncol : (colmax <= _ncol ? colmax : _ncol));
      for (i = 0; i < _nrow; i++){
        for (j = 0; j < cmax; j++){
          x = a(i, j);
          Rprintf("%5g  ", (fabs(x) < _AK_ZERO ? 0 : x));      
        }
        Rprintf("\n");
      }
    }
  }

  return;
}


template <typename dataType>
void
MatrixRect<dataType>::printI(const int &info,  const int &colmax) const
{
  int i, j, cmax;
  if (info) Rprintf("Matrix (nrow=%d, ncol=%d):\n", _nrow, _ncol);
  if (_length > 0){
    if (_nrow==1 || _ncol==1){
      cmax = (colmax < 0 ? _length : (colmax <= _length ? colmax : _length));
      for (i = 0; i < cmax; i++){
	Rprintf("%5d  ", _a[i]); 
      }
      Rprintf("\n");
    }
    else{
      cmax = (colmax < 0 ? _ncol : (colmax <= _ncol ? colmax : _ncol));
      for (i = 0; i < _nrow; i++){
        for (j = 0; j < cmax; j++){
          Rprintf("%5d  ", a(i, j));      
        }
        Rprintf("\n");
      }
    }
  }

  return;
}


#endif
