/*** MatrixLT.h ***/

#ifndef _MATRIX_LOWER_TRIANGLE_H_
#define _MATRIX_LOWER_TRIANGLE_H_

#include <cmath>

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "Cholesky2.h"

template <typename dataType=double>
class MatrixLT
{
  private:
    int _nrow;          // number of rows     
    int _length;        // length of the array
    dataType* _a;       // array where the elements of the lower triangle are stored
    dataType* _atemp;   // array of the same length as _a, used to store temporary results (e.g., during inversion of the matrix)
    int* _diagI;        // indeces of diagonal elements of the matrix in the array _a, i.e., _a[_diagI[i]] = a[i,i]

  public:

  /***** Constructors and destructors *****/
    MatrixLT();

    MatrixLT(const int& nrow, const dataType* a, const int& type);
    MatrixLT(const int& nrow);
    MatrixLT(const int& nrow, const dataType* b);

    MatrixLT(const MatrixLT<dataType>& A);
    MatrixLT& operator=(const MatrixLT<dataType>& A);
    ~MatrixLT();

  /*****  Access functions to various components of the matrix *****/
    inline int
    nrow() const { return _nrow;}

    inline int
    ncol() const { return _nrow;}

    inline int
    length() const { return _length;}

    /* Direct access to the array allowing its modification */
    inline dataType*
    a() { return _a;}

    inline const dataType*
    aconst() const { return _a;}

    /* Index of the i-th diagonal element in the array _a */
    inline const int*
    diagI() const { return _diagI; }
 
    inline int
    diagI(const int& i) const
    {
      if (i < 0 || i >= _nrow){
        REprintf("i=%d, _nrow=%d\n", i, _nrow);
        throw returnR("MatrixLT.h: diagI error", 1);
      }
      return (_diagI[i]);
      // return (i * (2*_nrow - i + 1))/2;
    }
  
    inline dataType
    a(const int& i) const { 
      if (i < 0 || i >= _length) throw returnR("MatrixLT.h: Function 'a(i)' error", 1);
      return _a[i];
    }

    inline dataType
    a(const int& i,  const int& j) const {
      if (i < 0 || i >= _nrow || j < 0 || j >= _nrow) throw returnR("MatrixLT.h: Function 'a(i, j)' error", 1);
 
      int r, c;
      if (i >= j){
        r = i;
        c = j;
      }
      else{
        r = j;
        c = i;
      }
      return(_a[_diagI[c] + r - c]);
    }  

  /***** Utilities *****/
    bool
    is_finite() const;

    void
    sqrtDiag(dataType *result) const;

  /*****  Input-output functions *****/
    void
    fillBy(const dataType& val);

    void
    mat2array(dataType* a, const int& type) const;

    void
    array2mat(const dataType* a);

    void
    print(const int &info=1) const;

    void
    print4R() const;

    void
    printI(const int &info=1) const;


  /*****  Linear algebra *****/
    void
    add_b2diag(const dataType *b);

    void
    add_A2(const dataType* A);


  /*****  Decompositions *****/
    int
    cholesky(const int &backup=0);

    void
    chinv(const int &onlyCholInv=0);

    int
    chol_inv(const int &what=0, dataType *Li=NULL);
   
    int
    chol_solve(dataType *c, dataType *x, const int &add_u=0, const int &solve_u=0, dataType *u=NULL);

    int
    choleskyPD(int &Attempt, const int &nAttempt=10, const dataType &eps=0.01);

    int
    chol_invPD(int &Attempt,           const int &what=0,           dataType *Li=NULL, 
               const int &nAttempt=10, const dataType &eps=0.01);

    int
    chol_solvePD(int &Attempt,           dataType *c,               dataType *x, 
                 const int &nAttempt=10, const dataType &eps=0.01,
                 const int &add_u=0,     const int &solve_u=0,      dataType *u=NULL);

};  /* end of the class MatrixLT */

template <typename dataType>
void
tLxAxL(dataType *tLAL, dataType *tLA, const dataType *L,  const dataType *A,  const int *diagI, const int *p);

template <typename dataType>
void
LxtL(dataType *LtL, const dataType *L, const int *p);

  /*** Necessary for template instantiation with some compilers. ***/
#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || defined (EXPLICIT_TEMPLATE_INSTANTIATION)
#include "MatrixLT.cpp"
#endif                     /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif

