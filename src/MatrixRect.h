/*** MatrixRect.h ***/

#ifndef _MATRIX_RECTANGLE_H_
#define _MATRIX_RECTANGLE_H_

#include <cmath>

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

template <typename dataType=double>
class MatrixRect
{
  private:
    int _nrow;        // number of rows     
    int _ncol;        // number of columns
    int _length;      // length of the array
    dataType* _a;     // array where the elements of the matrix are stored in column major order

  public:

  /***** Constructors and destructors *****/
    MatrixRect();

    MatrixRect(const int& nrow, const int& ncol, const dataType* a, const int& type=0);
    MatrixRect(const int& nrow, const int& ncol);
    MatrixRect(const int& nrow, const int& ncol, const dataType* b, const dataType* c, const int &type=0);

    MatrixRect(const MatrixRect<dataType>& A);
    MatrixRect<dataType>& operator=(const MatrixRect<dataType>& A);
    ~MatrixRect();

  /*****  Access functions to various components of the matrix *****/
    inline int
    nrow() const { return _nrow;}

    inline int
    ncol() const { return _ncol;}

    inline int
    length() const { return _length;}

    /* Direct access to the array allowing its modification */
    inline dataType*
    a() { return _a;}

    inline const dataType*
    aconst() const { return _a;}
  
    inline dataType
    a(const int& i) const { 
      if (i < 0 || i >= _length) throw returnR("MatrixRect.h: Function 'a(i)' error", 1);
      return _a[i];
    }

    inline dataType
    a(const int& i,  const int& j) const {
      if (i < 0 || i >= _nrow || j < 0 || j >= _ncol) throw returnR("MatrixRect.h: Function 'a(i, j)' error", 1);
      return(_a[j*_nrow + i]);
    }  

  /***** Utilities *****/
    bool
    is_finite() const;

    void
    copy(const MatrixRect<dataType> *b,  const int &length,  const int &bfrom=0,  const int &thisfrom=0);
  
    void
    rNorm(const double& mean=0,  const double& sd=1);

  /*****  Input-output functions *****/
    void
    fillBy(const dataType& val);

    void
    mat2array(dataType* a, const int& type=0) const;

    void
    vec2array(dataType *a, const int &start, const int &stop) const;

    void
    print(const int &info=1, const int &colmax=-1) const;

    void
    printI(const int &info=1, const int &colmax=-1) const;


  /*****  Linear algebra and other operations *****/
    dataType
    xTx() const;

    void
    Ab(const MatrixRect<dataType>* A, const MatrixRect<dataType>* b, const int& bstart=0);

    void
    Ab2(const MatrixRect<dataType>* A, const dataType* b);

    void
    bA(const MatrixRect<dataType>* b, const MatrixRect<dataType>* A, const int& bstart=0);

    void
    bA2(const MatrixRect<dataType>* b, const MatrixRect<dataType>* A, const int& bstart=0);

    // void
    // bA2row(const int& row, const MatrixRect<dataType>* b, const MatrixRect<dataType>* A, const int& bstart=0);

    void
    BAcolProd(const MatrixRect<dataType> *B, const MatrixRect<int> *ni, const MatrixRect<dataType> *A, const int &bstart=0);

    void
    BAcolProd2(const MatrixRect<dataType> *B, const MatrixRect<int> *ni, const MatrixRect<dataType> *A, const int &bstart=0);

    void
    add_A(const MatrixRect<dataType>* A);

    void
    add_A2(const dataType* A);

    void
    subtract_A(const MatrixRect<dataType>* A);

    void
    subtract_A2(const dataType* A);

    void
    add_const(const dataType &c);

    void
    add_const(const dataType &c, const int &start, const int &stop);

    void
    b_plus_rowsA(const dataType *b,  const MatrixRect<dataType> *A);

    void
    square();

    void
    square(const dataType *b);

    void
    divide(const dataType &b);

    dataType
    sum() const;

    dataType
    max() const;

    dataType
    min() const;

    void
    mean(dataType *Mean,  const int &margin) const;

    bool
    anyNonNeg() const;

};  /* end of the class MatrixRect */

  /*** Necessary for template instantiation with some compilers. ***/
#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || defined (EXPLICIT_TEMPLATE_INSTANTIATION)
#include "MatrixRect.cpp"
#endif                     /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif

