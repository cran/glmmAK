/*** MatrixUtil.h ***/

#ifndef _MATRIX_UTILITIES_H_
#define _MATRIX_UTILITIES_H_

#include <cmath>

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

#include "MatrixRect.h"
#include "MatrixLT.h"

template <typename dataType>
void
Ab(MatrixRect<dataType> *c, const MatrixLT<dataType> *A, const MatrixRect<dataType> *b, const int &bstart=0);

template <typename dataType>
void
Ab2(dataType *c, const MatrixLT<dataType> *A, const dataType *b);

template <typename dataType>
void
NSampleVar(MatrixLT<dataType> *Var, dataType *b_alpha, const dataType *b, const dataType *alpha, const int &N);

  /*** Necessary for template instantiation with some compilers. ***/
#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || defined (EXPLICIT_TEMPLATE_INSTANTIATION)
#include "MatrixUtil.cpp"
#endif                     /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif
