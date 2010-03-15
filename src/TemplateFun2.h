/*** TemplateFun2.h ***/

#ifndef _TEMPLATE_FUN_TWEE_H_
#define _TEMPLATE_FUN_TWEE_H_

#include <R.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "AK_Error.h"


template <typename dd>
void
writeToFile_1(const dd* array,      const int& nC,        std::ofstream& ofile,
              const int& prec = 6,  const int& width = 1);

template <typename dd>
void
writeToFile2_1(const dd* array1,   const int& nC1,
               const dd* array2,   const int& nC2,         std::ofstream& ofile,
               const int& prec=6,  const int& width=1);

template <typename dd>
void
writeToFile2_1_diag(const dd* array1,  const int& nC1,
                    const dd* array2,  const int& nC2,    std::ofstream& ofile,
                    const int& prec,   const int& width);


  /*** Necessary for template instantiation with some compilers. ***/
#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || defined (EXPLICIT_TEMPLATE_INSTANTIATION)
#include "TemplateFun2.cpp"
#endif                     /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif

