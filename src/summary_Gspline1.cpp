/*** summary_Gspline1.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  18/01/2007
//
//         summary_Gspline1:  18/01/2007
//            eval_Gspline1:  18/01/2007
//    readWeightsFromFiles1:  19/01/2007
//
// PURPOSE: Summary for the sample from the UNIVARIATE G-spline
//          (pointwise posterior mean and/or quantiles)
//
/* ********************************************************************************* */

#include "summary_Gspline1.h"

namespace summary_Gspline1A{

extern "C"{

// value:                    If (retValue || nprob) value[ngrid, niter] = values of all evaluated G-splines
//                           Else                   value[ngrid]        = values of the last evaluated G-spline
//                           
// summary[ngrid, 1+nprob]:  Pointwise averages
//                           + quantiles (if nprob>0)
// prob[nprob]:              Probabilities of required quantiles
// nprob:                    Number of quantiles required
// retValue:                 Do we want to return the values of G-spline at each iteration?
//
// grid[ngrid]:              Grid of values where the G-spline is to be evaluated
// ngrid:                    Number of grid values
//
// standard:                 If <> 0 then values of the standardized (zero-mean, unit-variance) G-splines are computed
//
// knots[nknots]:            Knots
// sigma[nknots]:            Basis standard deviations (may differ from component to component)
// nknots:                   Number of mixture components
//
// intcpt[niter]:          Sampled values of the G-spline intercepts
// tau[niter]:             Sampled values of the G-spline scale parameters
//
// niter:                    Number of sampled G-splines to be evaluated
// wpath:                    Path to the file with (log-)weights
// skip:                     How many rows are to be skipped in the file with (log-)weights
//                           If header is included in the file, then skip should be >= 1
// by:                       Possible thinning of the file with (log-)weights
// logw:                     0 if weights are given directly, <>0 if log-weights are supplied
// nwrite:                   Frequency of informing the user about the progress

void
summary_Gspline1(double *value,           double *summary,         
                 const double *prob,      const int *nprob,        const int *retValue,
                 const double *grid,      const int *ngrid,        const int *standard,
                 const double *knots,     const double *sigma,     const int *nknots,
                 const double *intcpt,    const double *tau,       const int *niter,        
                 char **wpath,            const int *skip,         const int *by,
                 const int *logw,         const int *nwrite,       int *err)
{
  try{
    *err = 0;

    const double _ZERO_ = 0.0;

    int i, ix, nwrite2;
    double *valueP, *sumP;
    const double *intcptP, *tauP;

    /*** Some input checks ***/
    if (*niter <= 0){
      REprintf("niter=%d\n", *niter);
      throw returnR("Error in summary_Gspline1.cpp: summary_Gspline1(). niter must be positive", 1);
    }
    if (*skip < 0){
      REprintf("skip=%d\n", *skip);
      throw returnR("Error in summary_Gspline1.cpp: summary_Gspline1(). Incorrect value of skip", 1);
    }
    if (*by <= 0 || *by > *niter){
      REprintf("niter=%d,  by=%d\n", *niter, *by);
      throw returnR("Error in summary_Gspline1.cpp: summary_Gspline1(). Incorrect value of by", 1);
    }
    if (*nknots <= 0){
      REprintf("nknots=%d\n", *nknots);
      throw returnR("Error in summary_Gspline1.cpp: summary_Gspline1(). nknots must be positive", 1);
    }
    nwrite2 = *nwrite > *niter ? *niter : *nwrite;
    for (i = 0; i < *nknots; i++){
      if (sigma[i] <= 0){
        REprintf("sigma[%d]=%g\n", i, sigma[i]);
        throw returnR("Error in summary_Gspline1.cpp: summary_Gspline1(). sigma must be positive", 1);
      }
    }

    /*** Allocate the space for weights, knots_tau, sigma_tau ***/
    double *weight    = (double*) calloc(*nknots, sizeof(double));
    double *knots_tau = (double*) calloc(*nknots, sizeof(double));
    double *sigma_tau = (double*) calloc(*nknots, sizeof(double));
    if (!weight || !knots_tau || !sigma_tau) throw returnR("Not enough memory available in summary_Gspline1.cpp: summary_Gspline1().", 99);

    /*** Open file with simulated G-spline (log-)weights and skip rows at the beginning of the file that are to be skipped ***/
    std::string swpath = *wpath;
    std::ifstream wfile;
    In_Output::open_File_toRead(wfile, swpath, *skip);    /** Possible skipping of the header should be included in the value of skip **/

    /*** Reset averages ***/
    AK_BLAS_LAPACK::fillArray(summary, &_ZERO_, *ngrid);

    /*** Read the first set of G-spline weights ***/
    summary_Gspline1A::readWeightsFromFiles1(weight, 0, *skip, nknots, wfile, swpath);
    valueP  = value;
    intcptP = intcpt;
    tauP    = tau;
    //Rprintf("*** Weights read:\n");
    //AK_BLAS_LAPACK::printArray(weight, *nknots);    

    /*** Compute the value of the first G-spline ***/
    summary_Gspline1A::eval_Gspline1(summary, valueP, weight, knots_tau, sigma_tau, grid, ngrid, standard,
                                      knots, sigma, nknots, intcptP, tauP, logw);

    /*** Compute the values of the remaining G-splines ***/
    int backs      = 0;
    int by_1       = *by - 1;
    int jump_value = (*retValue || *nprob) ? *ngrid : 0;

    Rprintf("Values evaluated: ");
    for (int iter = 1; iter < *niter; iter++){    
      summary_Gspline1A::readWeightsFromFiles1(weight, by_1, iter, nknots, wfile, swpath);
      valueP += jump_value;
      intcptP++;
      tauP++;
      summary_Gspline1A::eval_Gspline1(summary, valueP, weight, knots_tau, sigma_tau, grid, ngrid, standard,
                                        knots, sigma, nknots, intcptP, tauP, logw);

      if (!((iter+1) % nwrite2) || (iter+1) == *niter){
        for (i = 0; i < backs; i++) Rprintf("\b");
        Rprintf("%d", iter+1);
        backs = int(log10(double(iter+1))) + 1;
      }
    }      
    Rprintf("\n");

    /*** Close files with simulated G-splines (log-)weights ***/
    wfile.close();

    /*** Cleaning ***/
    free(weight);
    free(knots_tau);
    free(sigma_tau);

    /*** MCMC averages ***/
    sumP = summary;
    for (ix = 0; ix < *ngrid; ix++){
      *sumP /= *niter;
      sumP++;
    }

    /*** Quantiles (if required) ***/
    if (*nprob){
      Rprintf("Computing quantiles\n");
      Quantile::Quantile(sumP, value, ngrid, niter, prob, nprob);
    }
 
    return;
  }
  catch(returnR rr){
    *err = rr.errflag();
    return;
  }
}


/***** eval_Gspline1:  Evaluate a UNIVARIATE Gspline in a grid of values     *****/
/***                                                                           ***/
/*** ========================================================================= ***/
//  
// g(y) = sum[j](w[j]*(sigma[j]*tau)^(-1)*phi((y-intcpt-tau*mu[j])/(sigma[j]*tau)))
//  
// average[ngrid]:     OUTPUT:  INPUT plus value of the G-spline evaluated in a grid
// value[ngrid]:       OUTPUT:  Value of the G-spline evaluated in a grid
// weight[nknots]:     INPUT:   Weights or log-weights
//                     OUTPUT:  Weights (still multiplied by sumexpa)
// knots_tau[nknots]:  OUTPUT (if !standard):  Values of knots*tau 
// sigma_tau[nknots]:  OUTPUT (if !standard):  Values of sigma[j]*tau
//
// standard:           If <> 0 then standardized version (zero-mean, unit-variance) of the G-spline
//                     is computed
//
// sigma[nknots]:      INPUT:   Basis standard deviations for each component
//
// ======================================================================================
//
// REMARK:  EY    = mean_g = sum(w[j]*mu[j])
//         var(Y) = sum{w[j]*(sigma[j]^2 + mu[j]^2)} - mean_g^2
//                = sum{w[j]*(sigma[j]^2 + mu[j]^2 - mean_g^2)}
//                = sum{w[j]*(sigma[j]^2 + (mu[j] - mean_g)^2)}
//
void
eval_Gspline1(double *average,           double *value,            
              double *weight,            double *knots_tau,        double *sigma_tau,
              const double *grid,        const int *ngrid,         const int *standard,
              const double *knots,       const double *sigma,      const int *nknots,
              const double *intcpt,      const double *tau,        const int *logw)
{
  static int ix, k;
  static double *wP, *knots_tauP, *sigma_tauP, *averP, *valP;  
  static const double *knotsP, *sigmaP, *gridP;
  static double val_gk, mean_g, sd_g, knot_mean;
  static double sumexpa;

  /*** Compute exp(log-weights) if necessary ***/
  /*** (do not divide them by sumexpa)       ***/ 
  sumexpa = 0.0;
  if (*logw){
    wP = weight;
    for (k = 0; k < *nknots; k++){
      *wP = exp_AK(*wP);
      sumexpa += *wP;
      wP++;      
    }
  }
  else{              /** only sum-up the weights **/
    wP = weight;
    for (k = 0; k < *nknots; k++){
      sumexpa += *wP;
      wP++;      
    }
  }

  //Rprintf("*** Weights (multiplied by %g):\n", sumexpa);
  //AK_BLAS_LAPACK::printArray(weight, *nknots);
  
  if (*standard){
    /*** Compute the mean of the G-spline (without added intercept) ***/
    mean_g = 0.0;
    wP     = weight;
    knotsP = knots;
    for (k = 0; k < *nknots; k++){
      mean_g += (*wP)*(*knotsP);                                         // += w[k]*mu[k]
      wP++;
      knotsP++;
    }
    mean_g /= sumexpa;

    /*** Compute the variance of the G-spline (not multiplied by the overall variance) ***/
    sd_g   = 0.0;
    wP     = weight;
    knotsP = knots;
    sigmaP = sigma;
    for (k = 0; k < *nknots; k++){
      knot_mean = *knotsP - mean_g;
      sd_g     += *wP * (*sigmaP * *sigmaP + knot_mean * knot_mean);
      wP++;
      knotsP++;
      sigmaP++;
    }
    sd_g /= sumexpa;

    /*** Compute the mean and the variance of the G-spline (without added intercept) ***/
    /*** NUMERICALLY LESS STABLE POSSIBILITY                                         ***/
    //mean_g = 0.0;
    //sd_g   = 0.0;
    //wP     = weight;
    //knotsP = knots;
    //sigmaP = sigma;
    //for (k = 0; k < *nknots; k++){
    //  mean_g += (*wP)*(*knotsP);                                         // += w[k]*mu[k]
    //  sd_g   += (*wP)*((*knotsP)*(*knotsP) + (*sigmaP)*(*sigmaP));       // += w[k]*(mu[k]^2 + sigma[k]^2)
    //  wP++;
    //  knotsP++;
    //  sigmaP++;
    //}
    //mean_g /= sumexpa;
    //sd_g   /= sumexpa;    
    //sd_g -= mean_g*mean_g;

    if (sd_g <= 0){
      throw returnR("Error in summary_Gspline1.cpp: eval_Gspline1. Variance of the G-spline is not positive", 1);
    }
    sd_g = sqrt(sd_g);

    //sd_g *= (*tau);                // Taking into account the scale tau
    //mean_g *= (*tau);              // Taking into account the scale tau
    //mean_g += (*intcpt);           // Taking into account the intercept

    /*** Compute the value of the G-spline in the grid ***/
    valP  = value;
    averP = average;
    gridP = grid;
    for (ix = 0; ix < *ngrid; ix++){
      *valP = 0.0;
      knotsP = knots;
      sigmaP = sigma;
      wP     = weight;
      for (k = 0; k < *nknots; k++){   
        val_gk = (sd_g*(*gridP) + mean_g - (*knotsP))/(*sigmaP);
        val_gk = -0.5 * val_gk * val_gk;
        val_gk = val_gk < _AK_EMIN ? 0.0 : (*wP)*exp(val_gk)/(*sigmaP);
        *valP += val_gk;
        knotsP++;
        sigmaP++;
        wP++;                
      }
      *valP /= sumexpa;
      *valP *= sd_g;
      *valP *= _AK_invSQRT_TWO_PI;
      *averP += *valP;
      valP++;
      averP++;
      gridP++;
    }
  }
  else{
    /*** Compute knots_tau and sigma_tau ***/  
    knots_tauP = knots_tau;
    sigma_tauP = sigma_tau;
    knotsP  = knots;
    sigmaP  = sigma;
    for (k = 0; k < *nknots; k++){
      *knots_tauP = *knotsP * (*tau);
      *sigma_tauP = *sigmaP * (*tau);
      knots_tauP++;
      sigma_tauP++;
      knotsP++;
      sigmaP++;
    }   

    /*** Compute the value of the G-spline in the grid ***/
    valP  = value;
    averP = average;
    gridP = grid;
    for (ix = 0; ix < *ngrid; ix++){
      *valP = 0.0;
      knots_tauP = knots_tau;
      sigma_tauP = sigma_tau;
      wP         = weight;
      for (k = 0; k < *nknots; k++){   
        val_gk = ((*gridP) - (*intcpt) - (*knots_tauP))/(*sigma_tauP);
        val_gk = -0.5 * val_gk * val_gk;
        val_gk = val_gk < _AK_EMIN ? 0.0 : (*wP)*exp(val_gk)/(*sigma_tauP);
        *valP += val_gk;
        knots_tauP++;
        sigma_tauP++;
        wP++;                
      }
      *valP /= sumexpa;
      *valP *= _AK_invSQRT_TWO_PI;
      *averP += *valP;
      valP++;
      averP++;
      gridP++;
    }
  }
  
  return;
}


}  /*** end of extern "C" ***/

  
/*** readWeightsFromFiles1:  Read (log-)weights from the file (G-spline)                                       ***/
/***                                                                                                           ***/
/***     It is assumed that each row of the file contains (log-)weights for all G-spline components            ***/
/***                                                                                                           ***/
/*** ========================================================================================================= ***/
//
// weight[nknots]:  OUTPUT: Read (log-)weights
//                          There should always be nknots values on each row of the file wfile
//
// skip:   Number of rows that should be skipped before reading real data
// row:    How many rows of the data have already been read before the function call?
//         (only used to generate possible error messages)
//
void
readWeightsFromFiles1(double *weight,  
                      const int &skip,       const int &row,            const int *nknots,
                      std::ifstream &wfile,  const std::string& swpath)
{
  try{
    static int j, ihelp;
    static char ch;
    static std::string errmes;
    static double *wP;

    //Rprintf("*** skip=%d,  nknots=%d\n", skip, *nknots);

    /**  Skip rows that are to be skipped  **/
    for (j = 0; j < skip; j++){
      wfile.get(ch);        
      while (ch != '\n') wfile.get(ch);
    }

    /** Read G-spline weights **/
    if (wfile.eof()){
      ihelp = row + 1;
      errmes = std::string("Error in summary_Gspline1.cpp: readWeightsFromFiles1(): Reached end of file ") + swpath + " before "
               + char(ihelp) + std::string(" sets of G-spline weights were read.");
      throw returnR(errmes, 99);
    }

    wP = weight;
    for (j = 0; j < *nknots; j++){
      wfile >> *wP;
      //Rprintf("w[%d]=%g,  ", j, *wP);
      wP++;
    }
    wfile.get(ch);                 
    while (ch != '\n') wfile.get(ch);

    return;
  }
  catch(returnR){
    throw;
  }  
}

}  /*** end of the namespace summary_Gspline1A ***/
