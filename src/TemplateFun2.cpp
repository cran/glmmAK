/*** TemplateFun2.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  14/09/2006
//
//         writeToFile_1:  14/09/2006
//        writeToFile2_1:  27/09/2006
//   writeToFile2_1_diag:  02/11/2006
//
//
// PURPOSE: Some template functions.
//          * partially taken from templatefun.cpp and templatefun_GS.cpp of the 'bayesSurv' package
//
/* ********************************************************************************* */

#ifndef _TEMPLATE_FUN_TWEE_CPP_
#define _TEMPLATE_FUN_TWEE_CPP_

#include "TemplateFun2.h"

// =============================================================================================
// ***** writeToFile_1: Function to write a numeric array as one row to a file             *****
// =============================================================================================
template <typename dd>
void
writeToFile_1(const dd* array,  const int& nC,     std::ofstream& ofile,
              const int& prec,  const int& width)
{
  const dd *arrP = array;
  
  for (int j = 0; j < nC; j++){
    if (*arrP >= FLT_MAX){                                                 
      ofile << std::setw(width) << "1e50";
      ofile << "   ";
    }
    else{
      if (*arrP < 1 && *arrP > -1 && *arrP != 0){
        ofile << std::scientific << std::setw(width) << std::setprecision(prec) << *arrP;
        ofile << "   ";
      }
      else{
        ofile << std::fixed << std::setw(width) << std::setprecision(prec) << *arrP;
        ofile << "   ";
      }
    }
    arrP++;
  }
  ofile << std::endl;

  return;
}


// =============================================================================================
// ***** writeToFile2_1: Function to write 2 numeric arrays as one row to a file           *****
// =============================================================================================
template <typename dd>
void
writeToFile2_1(const dd* array1,  const int& nC1,
               const dd* array2,  const int& nC2,    std::ofstream& ofile,
               const int& prec,   const int& width)
{
  int j;
  const dd *arrP = array1;  
  for (j = 0; j < nC1; j++){
    if (*arrP >= FLT_MAX){                                                 
      ofile << std::setw(width) << "1e50";
      ofile << "   ";
    }
    else{
      if (*arrP < 1 && *arrP > -1 && *arrP != 0){
        ofile << std::scientific << std::setw(width) << std::setprecision(prec) << *arrP;
        ofile << "   ";
      }
      else{
        ofile << std::fixed << std::setw(width) << std::setprecision(prec) << *arrP;
        ofile << "   ";
      }
    }
    arrP++;
  }

  arrP = array2;  
  for (j = 0; j < nC2; j++){
    if (*arrP >= FLT_MAX){                                                 
      ofile << std::setw(width) << "1e50";
      ofile << "   ";
    }
    else{
      if (*arrP < 1 && *arrP > -1 && *arrP != 0){
        ofile << std::scientific << std::setw(width) << std::setprecision(prec) << *arrP;
        ofile << "   ";
      }
      else{
        ofile << std::fixed << std::setw(width) << std::setprecision(prec) << *arrP;
        ofile << "   ";
      }
    }
    arrP++;
  }

  ofile << std::endl;

  return;
}


// =============================================================================================
// ***** writeToFile2_1_diag: Function to write 2 numeric arrays as one row to a file      *****
//       * arrays to be written are assumed to form a diagonal in the symmetric matrix    
//         whose lower triangles are stored in array1 and array2
// =============================================================================================
//
// array1[LT(nC1)]:    lower triangle of the symmetric matrix stored in COLUMN major order
// nC1:                number of rows of the matrix stored in array1
// array2[LT(nC2)]:    lower triangle of the symmetric matrix stored in COLUMN major order
// nC2:                number of rows of the matrix stored in array1
//
//
template <typename dd>
void
writeToFile2_1_diag(const dd* array1,  const int& nC1,
                    const dd* array2,  const int& nC2,    std::ofstream& ofile,
                    const int& prec,   const int& width)
{
  int j;
  const dd *arrP = array1;  
  for (j = nC1; j > 0; j--){
    if (*arrP >= FLT_MAX){                                                 
      ofile << std::setw(width) << "1e50";
      ofile << "   ";
    }
    else{
      if (*arrP < 1 && *arrP > -1 && *arrP != 0){
        ofile << std::scientific << std::setw(width) << std::setprecision(prec) << *arrP;
        ofile << "   ";
      }
      else{
        ofile << std::fixed << std::setw(width) << std::setprecision(prec) << *arrP;
        ofile << "   ";
      }
    }
    arrP += j;
  }

  arrP = array2;  
  for (j = nC2; j > 0; j--){
    if (*arrP >= FLT_MAX){                                                 
      ofile << std::setw(width) << "1e50";
      ofile << "   ";
    }
    else{
      if (*arrP < 1 && *arrP > -1 && *arrP != 0){
        ofile << std::scientific << std::setw(width) << std::setprecision(prec) << *arrP;
        ofile << "   ";
      }
      else{
        ofile << std::fixed << std::setw(width) << std::setprecision(prec) << *arrP;
        ofile << "   ";
      }
    }
    arrP += j;
  }

  ofile << std::endl;

  return;
}


#endif

