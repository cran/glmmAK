/*** summary_BiGspline.h ***/

#ifndef _SUMMARY_BI_GSPLINE_H_
#define _SUMMARY_BI_GSPLINE_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_BLAS_LAPACK.h"
#include "In_Output2.h"
#include "Quantile2.h"

#include "summary_Gspline1.h"

namespace summary_BiGsplineA{

extern "C"{

  void
  eval_BiGspline(double *average,           double *average0,         double *average1,
                 double *value,             double *value0,           double *value1,
                 double *weight,            double *weight0,          double *weight1,
                 double *knots_tau0,        double *knots_tau1,
                 double *sigma_tau0,        double *sigma_tau1,
                 const int *ind_w_effect,
                 const double *grid0,       const double *grid1,      const int *ngrid,         
                 const int *standard,
                 const double *knots0,      const double *knots1,       
                 const double *sigma0,      const double *sigma1,       
                 const int *nknots,
                 const double *intcpt,      const double *tau,        const int *logw);

}  /*** end of extern "C" ***/

void
readWeightsFromFilesBi(int *k_effect,           int *ind_w_effect,           double *weight,  
                       const int &skip,         const int &row,              const int *total_length,
                       std::ifstream &indfile,  const std::string& sindpath,
                       std::ifstream &wfile,    const std::string& swpath);

}  /*** end of the namespace summary_BiGsplineA ***/

#endif
