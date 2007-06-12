/*** summary_Gspline1.h ***/

#ifndef _SUMMARY_GSPLINE_EEN_H_
#define _SUMMARY_GSPLINE_EEN_H_

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

namespace summary_Gspline1A{

extern "C"{

  void
  summary_Gspline1(double *value,           double *summary,         
                   const double *prob,      const int *nprob,        const int *retValue,
                   const double *grid,      const int *ngrid,        const int *standard,
                   const double *knots,     const double *sigma,     const int *nknots,
                   const double *intcpt,    const double *tau,       const int *niter,        
                   char **wpath,            const int *skip,         const int *by,
                   const int *logw,         const int *nwrite,       int *err);

  void
  eval_Gspline1(double *average,           double *value,            
                double *weight,            double *knots_tau,        double *sigma_tau,
                const double *grid,        const int *ngrid,         const int *standard,
                const double *knots,       const double *sigma,      const int *nknots,
                const double *intcpt,      const double *tau,        const int *logw);
}

void
readWeightsFromFiles1(double *weight,  
                      const int &skip,       const int &row,            const int *nknots,
                      std::ifstream &wfile,  const std::string& swpath);

}  /*** end of the namespace summary_Gspline1A ***/

#endif
