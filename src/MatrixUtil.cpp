/*** MatrixUtil.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  13/09/2006
//
//                             Ab:  13/09/2006
//                            Ab2:  18/09/2006
//                     NSampleVar:  20/09/2006
//
// PURPOSE: Functions to manipulate matrices
//
/* ********************************************************************************* */

#ifndef _MATRIX_UTILITIES_CPP_
#define _MATRIX_UTILITIES_CPP_

#include "MatrixUtil.h"


/***** Ab: Product of a symmetric (LT) matrix and a vector *****/
/*                                                                                   */
/* * (*c) can be both row or column vector                                           */
/* * length of (*c) must be equal to A->nrow()                                       */
/* * b can be both row or column vector                                              */
/* * from b we use b[bstart],...,b[bstart+A->ncol()-1]                               */
/* * default for bstart is 0                                                         */
/*                                                                                   */
/* OUTPUT: c                                                                         */
/* INPUT:  A, b, bstart                                                              */
/*                                                                                   */
/* --------------------------------------------------------------------------------- */
template <typename dataType>
void
Ab(MatrixRect<dataType> *c, const MatrixLT<dataType> *A, const MatrixRect<dataType> *b, const int &bstart)
{
  int i, j;
  const dataType *AP, *bP0, *bP1;
  dataType *cP0, *cP1;

  cP0 = c->a();
  for (i = 0; i < c->length(); i++){
    *cP0 = 0;
    cP0++;
  }

  bP0 = b->aconst() + bstart;
  cP0 = c->a();
  AP = A->aconst();
  for (j = 0; j < A->ncol(); j++){
    *cP0 += (*AP) * (*bP0);
    AP++;
    bP1 = bP0 + 1;
    cP1 = cP0 + 1;
    for (i = j+1; i < A->nrow(); i++){
      *cP0 += (*AP) * (*bP1);
      *cP1 += (*AP) * (*bP0);
      AP++;
      bP1++;
      cP1++;
    }
    bP0++;
    cP0++;
  }

  return;
}


/***** Ab2: Product of a symmetric (LT) matrix and a vector which is given as an array *****/
/*                                                                                         */
/* * (*c) can be both row or column vector                                                 */
/* * length of (*c) must be equal to A->nrow()                                             */
/* * b can be both row or column vector                                                    */
/* * from b we use b[bstart],...,b[bstart+A->ncol()-1]                                     */
/* * default for bstart is 0                                                               */
/*                                                                                         */
/*   c: array of length A->nrow()                                                          */
/*   b: array of length A->ncol()                                                          */
/*                                                                                         */
/* OUTPUT: c                                                                               */
/* INPUT:  A, b                                                                            */
/*                                                                                         */
/* --------------------------------------------------------------------------------------- */
template <typename dataType>
void
Ab2(dataType *c, const MatrixLT<dataType> *A, const dataType *b)
{
  int i, j;
  const dataType *AP, *bP0, *bP1;
  dataType *cP0, *cP1;

  cP0 = c;
  for (i = 0; i < A->nrow(); i++){
    *cP0 = 0;
    cP0++;
  }

  bP0 = b;
  cP0 = c;
  AP = A->aconst();
  for (j = 0; j < A->ncol(); j++){
    *cP0 += (*AP) * (*bP0);
    AP++;
    bP1 = bP0 + 1;
    cP1 = cP0 + 1;
    for (i = j+1; i < A->nrow(); i++){
      *cP0 += (*AP) * (*bP1);
      *cP1 += (*AP) * (*bP0);
      AP++;
      bP1++;
      cP1++;
    }
    bP0++;
    cP0++;
  }

  return;
}


/***** NSampleVar:  Compute sum((b[i] - alpha)*(b[i] - alpha)'), i=1,...,N             *****/
/*                                                                                         */
/*       Var:  result, symmetric matrix q x q                                              */
/*   b_alpha:  working array of length q*N                                                 */
/*         b:  matrix  q x N, stored in an array                                           */
/*     alpha:  vector of length q                                                          */
/*                                                                                         */
/*                                                                                         */
/* --------------------------------------------------------------------------------------- */
template <typename dataType>
void
NSampleVar(MatrixLT<dataType> *Var, dataType *b_alpha, const dataType *b, const dataType *alpha, const int &N)
{
  static int i, j, k;
  static const dataType *bP, *bP2, *alphaP;
  static dataType *fillP;

  /*** b[i] - alpha ***/
  bP = b;
  fillP = b_alpha;
  for (j = 0; j < N; j++){
    alphaP = alpha;
    for (i = 0; i < Var->nrow(); i++){
      *fillP = *bP - (*alphaP);
      bP++;
      fillP++;
      alphaP++;
    }
  }

  /*** sum((b[i] - alpha)*(b[i] - alpha)') ***/
  Var->fillBy(0);
  bP = b_alpha;
  for (j = 0; j < N; j++){
    fillP = Var->a();
    for (k = 0; k < Var->ncol(); k++){
      bP2 = bP;
      for (i = k; i < Var->nrow(); i++){
        *fillP += *bP * (*bP2);
        fillP++;
        bP2++;
      }
      bP++;
    }
  }

  return;
}


#endif
