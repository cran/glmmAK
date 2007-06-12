/*** Cholesky2.h ***/

#ifndef _CHOLESKY_TWEE_H_
#define _CHOLESKY_TWEE_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"


template <typename dataType>
int
cholesky0(dataType *_a,  dataType *_atemp,  const int &_nrow,  const int *_diagI,  const int &backup=0);

template <typename dataType>
void
chinv0(dataType *_a,  dataType *_atemp,  const int &_nrow,  const int *_diagI,  const int &onlyCholInv=0);

template <typename dataType>
int
chol_inv0(dataType *_a,  dataType *_atemp,  const int &_nrow,  const int *_diagI,  const int &what=0,  dataType *Li=NULL);

template <typename dataType>
int
chol_solve0(dataType *_a,        dataType *_atemp,      const int &_nrow,  const int *_diagI,
            dataType *c,         dataType *x, 
            const int &add_u=0,  const int &solve_u=0,  dataType *u=NULL);

template <typename dataType>
int
choleskyPD0(dataType *_a,   dataType *_atemp,        const int &_nrow,         const int *_diagI,  
            int &Attempt,   const int &nAttempt=10,  const dataType &eps=0.01);

template <typename dataType>
int
chol_invPD0(dataType *_a,            dataType *_atemp,     const int &_nrow,  const int *_diagI,
            int &Attempt,            const int &what,      dataType *Li=NULL, 
            const int &nAttempt=10,  const dataType &eps=0.01);

template <typename dataType>
int
chol_solvePD0(dataType *_a,            dataType *_atemp,          const int &_nrow,  const int *_diagI, 
              int &Attempt,            dataType *c,               dataType *x, 
              const int &nAttempt=10,  const dataType &eps=0.01,
              const int &add_u=0,      const int &solve_u=0,      dataType *u=NULL);

  /*** Necessary for template instantiation with some compilers. ***/
#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || defined (EXPLICIT_TEMPLATE_INSTANTIATION)
#include "Cholesky2.cpp"
#endif                     /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif
