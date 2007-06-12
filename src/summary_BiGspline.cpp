/*** summary_BiGspline.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  03/04/2007
//
//         summary_BiGspline:  04/04/2007
//            eval_BiGspline:  04/04/2007
//    readWeightsFromFilesBi:  04/04/2007
//
// PURPOSE: Summary for the sample from the BIVARIATE G-spline
//          (pointwise posterior mean and/or quantiles)
//          Also summary for marginal G-splines is computed
//
/* ********************************************************************************* */

#include "summary_BiGspline.h"

namespace summary_BiGsplineA{

extern "C"{

// value:                    If (retValue || nprob) value[ngrid[0]*ngrid[1], niter] = values of all evaluated G-splines
//                           Else                   value[ngrid[0]*ngrid[1]]        = values of the last evaluated G-spline
// value0:                   If (retValue || nprob) value[ngrid[0], niter] = values of all evaluated marginal G-splines (margin 0)
//                           Else                   value[ngrid[0]]        = values of the last evaluated marginal G-spline (margin 0)
// value1:                   If (retValue || nprob) value[ngrid[1], niter] = values of all evaluated marginal G-splines (margin 1)
//                           Else                   value[ngrid[1]]        = values of the last evaluated marginal G-spline (margin 1)
//                           
// summary[ngrid[0]*ngrid[1], 1+nprob]:  Pointwise averages
//                                       + quantiles (if nprob>0)
// summary0[ngrid[0], 1+nprob]:  Pointwise averages for margin 0
//                               + quantiles (if nprob>0)
// summary1[ngrid[1], 1+nprob]:  Pointwise averages for margin 1
//                               + quantiles (if nprob>0)
//
// prob[nprob]:              Probabilities of required quantiles
// nprob:                    Number of quantiles required
// retValue:                 Do we want to return the values of G-spline at each iteration?
//
// grid0[ngrid[0]]:          Grid of values for margin 0 where the G-spline is to be evaluated
// grid0[ngrid[1]]:          Grid of values for margin 1 where the G-spline is to be evaluated
// ngrid[2]:                 Numbers of grid values
//
// standard:                 If <> 0 then values of the standardized (zero-mean, unit-variance) G-splines are computed
//
// knots0[nknots[0]]:        Knots for margin 0
// knots1[nknots[1]]:        Knots for margin 1
//
// sigma0[nknots[0]]:        Basis standard deviations for margin 0 (may differ from component to component)
// sigma1[nknots[1]]:        Basis standard deviations for margin 1 (may differ from component to component)
//
// nknots[2]:                Numbers of mixture components
//
// intcpt[2*niter]:          Sampled values of the G-spline intercepts
// tau[2*niter]:             Sampled values of the G-spline scale parameters
//
// niter:                    Number of sampled G-splines to be evaluated
// wpath:                    Path to the file with (log-)weights
// indpath:                  Path to the file with indeces of "non-zero" components
//
// skip:                     How many rows are to be skipped in the files wpath, indpath
//                           If header is included in the file, then skip should be >= 1
// by:                       Possible thinning of the file with (log-)weights
//
// logw:                     0 if weights are given directly, <>0 if log-weights are supplied in wpath
// is_indfile:               If <> 0 then file indpath is used
//                           Else         file indpath is ignored and it is assumed that wpath contains on each row
//                                        the weights of all (zero as well) components    
//
// nwrite:                   Frequency of informing the user about the progress
//
void
summary_BiGspline(double *value,           double *value0,          double *value1,
                  double *summary,         double *summary0,        double *summary1,
                  const double *prob,      const int *nprob,        const int *retValue,
                  const double *grid0,     const double *grid1,     const int *ngrid,        const int *standard,
                  const double *knots0,    const double *knots1,    
                  const double *sigma0,    const double *sigma1,    const int *nknots,
                  const double *intcpt,    const double *tau,       const int *niter,        
                  char **wpath,            char **indpath,          const int *skip,         const int *by,
                  const int *logw,         const int *is_indfile,   const int *nwrite,       int *err)
{
  try{
    *err = 0;

    const double _ZERO_ = 0.0;
    const int _ONE_ = 1;

    int i, ix, nwrite2, total_length, total_ngrid;
    double *valueP, *valueP0, *valueP1, *sumP;
    const double *intcptP, *tauP;

    /*** Some input checks ***/
    if (*niter <= 0){
      REprintf("niter=%d\n", *niter);
      throw returnR("Error in summary_BiGspline.cpp: summary_BiGspline(). niter must be positive", 1);
    }
    if (*skip < 0){
      REprintf("skip=%d\n", *skip);
      throw returnR("Error in summary_BiGspline.cpp: summary_BiGspline(). Incorrect value of skip", 1);
    }
    if (*by <= 0 || *by > *niter){
      REprintf("niter=%d,  by=%d\n", *niter, *by);
      throw returnR("Error in summary_BiGspline.cpp: summary_BiGspline(). Incorrect value of by", 1);
    }
    if (nknots[0] <= 0 || nknots[1] <= 0){
      REprintf("nknots[0]=%d,  nknots[1]=%d\n", nknots[0], nknots[1]);
      throw returnR("Error in summary_BiGspline.cpp: summary_BiGspline(). nknots must be positive", 1);
    }
    total_length = nknots[0]*nknots[1];
    total_ngrid  = ngrid[0]*ngrid[1];
    nwrite2 = *nwrite > *niter ? *niter : *nwrite;
    for (i = 0; i < nknots[0]; i++){
      if (sigma0[i] <= 0){
        REprintf("sigma0[%d]=%g\n", i, sigma0[i]);
        throw returnR("Error in summary_BiGspline.cpp: summary_BiGspline(). sigma must be positive", 1);
      }
    }
    for (i = 0; i < nknots[1]; i++){
      if (sigma1[i] <= 0){
        REprintf("sigma1[%d]=%g\n", i, sigma1[i]);
        throw returnR("Error in summary_BiGspline.cpp: summary_BiGspline(). sigma must be positive", 1);
      }
    }

    /*** Allocate the space for weights, knots_tau, sigma_tau ***/
    double *weight     = (double*) calloc(total_length, sizeof(double));
    double *weight0    = (double*) calloc(nknots[0], sizeof(double));
    double *weight1    = (double*) calloc(nknots[1], sizeof(double));
    double *knots_tau0 = (double*) calloc(nknots[0], sizeof(double));
    double *knots_tau1 = (double*) calloc(nknots[1], sizeof(double));
    double *sigma_tau0 = (double*) calloc(nknots[0], sizeof(double));
    double *sigma_tau1 = (double*) calloc(nknots[1], sizeof(double));
    if (!weight || !weight0 || !weight1 || !knots_tau0 || !knots_tau1 || !sigma_tau0 || !sigma_tau1) throw returnR("Not enough memory available in summary_BiGspline.cpp: summary_BiGspline().", 99);

    /*** Allocate the space for ind_w_effect and fill it by 1 ***/
    int *ind_w_effect = (int*) calloc(total_length, sizeof(int));
    if (!ind_w_effect) throw returnR("Not enough memory available in summary_BiGspline.cpp: summary_BiGspline().", 99);
    AK_BLAS_LAPACK::fillArrayI(ind_w_effect, &_ONE_, total_length);

    /*** Open file with simulated G-spline (log-)weights and skip rows at the beginning of the file that are to be skipped ***/
    std::string swpath = *wpath;
    std::ifstream wfile;
    In_Output::open_File_toRead(wfile, swpath, *skip);    /** Possible skipping of the header should be included in the value of skip **/

    /*** File with indeces of non-zero weights and skip rows at the beginning of the file that are to be skipped ***/
    std::string sindpath = *indpath;
    std::ifstream indfile;

    /*** Reset averages ***/
    AK_BLAS_LAPACK::fillArray(summary, &_ZERO_, total_ngrid);
    AK_BLAS_LAPACK::fillArray(summary0, &_ZERO_, ngrid[0]);
    AK_BLAS_LAPACK::fillArray(summary1, &_ZERO_, ngrid[1]);


    int backs       = 0;
    int by_1        = *by - 1;
    int jump_value  = (*retValue || *nprob) ? total_ngrid : 0;
    int jump_value0 = (*retValue || *nprob) ? ngrid[0] : 0;
    int jump_value1 = (*retValue || *nprob) ? ngrid[1] : 0;
    int iter;
    int k_effect[1];

    if (*is_indfile){
      In_Output::open_File_toRead(indfile, sindpath, *skip);    /** Possible skipping of the header should be included in the value of skip **/

      /*** Read the first set of G-spline weights ***/
      summary_BiGsplineA::readWeightsFromFilesBi(k_effect, ind_w_effect, weight, 0, *skip, &total_length, indfile, sindpath, wfile, swpath);
      valueP  = value;
      valueP0 = value0;
      valueP1 = value1;
      intcptP = intcpt;
      tauP    = tau;
      //Rprintf("*** Weights read:\n");
      //AK_BLAS_LAPACK::printArray(weight, total_length);
      //Rprintf("ind_w_effect:\n");
      //for (int ii = 0; ii < total_length; ii++) Rprintf("%d(%d), ", ii, ind_w_effect[ii]);
      //Rprintf("\n");

      /*** Compute the value of the first G-spline ***/
      summary_BiGsplineA::eval_BiGspline(summary, summary0, summary1, valueP, valueP0, valueP1,
                                         weight, weight0, weight1, knots_tau0, knots_tau1, sigma_tau0, sigma_tau1, ind_w_effect,
                                         grid0, grid1, ngrid, standard,
                                         knots0, knots1, sigma0, sigma1, nknots, intcptP, tauP, logw);

      /*** Compute the values of the remaining G-splines ***/
      Rprintf("Values evaluated: ");
      for (iter = 1; iter < *niter; iter++){    
        summary_BiGsplineA::readWeightsFromFilesBi(k_effect, ind_w_effect, weight, by_1, iter, &total_length, indfile, sindpath, wfile, swpath);
        valueP  += jump_value;
        valueP0 += jump_value0;
        valueP1 += jump_value1;
        intcptP += 2;
        tauP    += 2;
        summary_BiGsplineA::eval_BiGspline(summary, summary0, summary1, valueP, valueP0, valueP1,
                                           weight, weight0, weight1, knots_tau0, knots_tau1, sigma_tau0, sigma_tau1, ind_w_effect,
                                           grid0, grid1, ngrid, standard,
                                           knots0, knots1, sigma0, sigma1, nknots, intcptP, tauP, logw);

        if (!((iter+1) % nwrite2) || (iter+1) == *niter){
          for (i = 0; i < backs; i++) Rprintf("\b");
          Rprintf("%d", iter+1);
          backs = int(log10(double(iter+1))) + 1;
        }
      }      
      Rprintf("\n");

      indfile.close();
    }
    else{           /** !(*is_indfile) **/
      /*** Read the first set of G-spline weights ***/
      summary_Gspline1A::readWeightsFromFiles1(weight, 0, *skip, &total_length, wfile, swpath);
      valueP  = value;
      valueP0 = value0;
      valueP1 = value1;
      intcptP = intcpt;
      tauP    = tau;
      //Rprintf("*** Weights read:\n");
      //AK_BLAS_LAPACK::printArray(weight, total_length);    

      /*** Compute the value of the first G-spline ***/
      summary_BiGsplineA::eval_BiGspline(summary, summary0, summary1, valueP, valueP0, valueP1,
                                         weight, weight0, weight1, knots_tau0, knots_tau1, sigma_tau0, sigma_tau1, ind_w_effect,
                                         grid0, grid1, ngrid, standard,
                                         knots0, knots1, sigma0, sigma1, nknots, intcptP, tauP, logw);

      /*** Compute the values of the remaining G-splines ***/
      Rprintf("Values evaluated: ");
      for (iter = 1; iter < *niter; iter++){    
        summary_Gspline1A::readWeightsFromFiles1(weight, by_1, iter, &total_length, wfile, swpath);
        valueP  += jump_value;
        valueP0 += jump_value0;
        valueP1 += jump_value1;
        intcptP += 2;
        tauP    += 2;
        summary_BiGsplineA::eval_BiGspline(summary, summary0, summary1, valueP, valueP0, valueP1,
                                           weight, weight0, weight1, knots_tau0, knots_tau1, sigma_tau0, sigma_tau1, ind_w_effect,
                                           grid0, grid1, ngrid, standard,
                                           knots0, knots1, sigma0, sigma1, nknots, intcptP, tauP, logw);

        if (!((iter+1) % nwrite2) || (iter+1) == *niter){
          for (i = 0; i < backs; i++) Rprintf("\b");
          Rprintf("%d", iter+1);
          backs = int(log10(double(iter+1))) + 1;
        }
      }      
      Rprintf("\n");
    }

    /*** Close file with simulated G-splines (log-)weights ***/
    wfile.close();

    /*** Cleaning ***/
    free(weight);               //Rprintf("weight CLEANED\n");
    free(weight0);              //Rprintf("weight0 CLEANED\n");
    free(weight1);              //Rprintf("weight1 CLEANED\n");
    free(knots_tau0);           //Rprintf("knots_tau0 CLEANED\n");
    free(knots_tau1);           //Rprintf("knots_tau1 CLEANED\n");
    free(sigma_tau0);           //Rprintf("sigma_tau0 CLEANED\n");
    free(sigma_tau1);           //Rprintf("sigma_tau1 CLEANED\n");
    free(ind_w_effect);         //Rprintf("ind_w_effect CLEANED\n");

    /*** MCMC averages and quantiles (if required) ***/
    sumP = summary;
    for (ix = 0; ix < total_ngrid; ix++){
      *sumP /= *niter;
      sumP++;
    }
    if (*nprob){
      Rprintf("Computing quantiles of the joint G-spline\n");
      Quantile::Quantile(sumP, value, &total_ngrid, niter, prob, nprob);
    }

    sumP = summary0;
    for (ix = 0; ix < ngrid[0]; ix++){
      *sumP /= *niter;
      sumP++;
    }
    if (*nprob){
      Rprintf("Computing quantiles of the marginal G-spline (margin 0)\n");
      Quantile::Quantile(sumP, value0, ngrid+0, niter, prob, nprob);
    }

    sumP = summary1;
    for (ix = 0; ix < ngrid[1]; ix++){
      *sumP /= *niter;
      sumP++;
    }
    if (*nprob){
      Rprintf("Computing quantiles of the marginal G-spline (margin 1)\n");
      Quantile::Quantile(sumP, value1, ngrid+1, niter, prob, nprob);
    }
 
    return;
  }
  catch(returnR rr){
    *err = rr.errflag();
    return;
  }
}


/***** eval_BiGspline:  Evaluate a BIVARIATE Gspline in a grid of values     *****/
/***                                                                           ***/
/*** ========================================================================= ***/
//  
// average[ngrid[0]*ngrid[1]]:  OUTPUT:  INPUT plus value of the G-spline evaluated in a grid
//                                       (array from a matrix in column major order)
// average0[ngrid[0]]:          OUTPUT:  INPUT plus value of the marginal G-spline (margin 0) evaluated in a grid
// average1[ngrid[1]]:          OUTPUT:  INPUT plus value of the marginal G-spline (margin 1) evaluated in a grid
//
// value[ngrid[0]*ngrid[1]]:    OUTPUT:  Value of the G-spline evaluated in a grid
//                                       (array from a matrix in column major order)
// value0[ngrid[0]]:            OUTPUT:  Value of the marginal G-spline (margin 0) evaluated in a grid
// value1[ngrid[1]]:            OUTPUT:  Value of the marginal G-spline (margin 1) evaluated in a grid
//
// weight[nknots[0]*nknots[1]]:     INPUT:   Weights or log-weights
//                                           Only the weights for which ind_w_effect is equal to 1 are considered,
//                                           remaining weights are assumed to be equal to zero 
//                                  OUTPUT:  Weights (still multiplied by sumexpa)
//                                           (array from a matrix in column major order)
//
// weight0[ngrid[0]]:  OUTPUT:  Marginal weights (still multiplied by sumexpa) for the first margin
// weight1[ngrid[1]]:  OUTPUT:  Marginal weights (still multiplied by sumexpa) for the second margin
//
// knots_tau0[nknots[0]]:  OUTPUT (if !standard):  Values of knots[0]*tau[0]
// knots_tau1[nknots[1]]:  OUTPUT (if !standard):  Values of knots[1]*tau[1] 
//
// sigma_tau0[nknots[0]]:  OUTPUT (if !standard):  Values of sigma0[j]*tau[0]
// sigma_tau1[nknots[1]]:  OUTPUT (if !standard):  Values of sigma1[j]*tau[1]
//
// ind_w_effect[nknots[0]*nknots[1]]:  0/1 indicators of "non-zero" components
//
// grid0[ngrid[0]]:    Grid values in the first margin
// grid1[ngrid[1]]:    Grid values in the second margin
// ngrid[2]:           Lengths of the grids
//
// standard:           If <> 0 then standardized version (zero-mean, unit-variance) of the G-spline
//                     is computed
//
// sigma0[nknots[0]]:      INPUT:   Basis standard deviations for each component in the first margin
// sigma1[nknots[1]]:      INPUT:   Basis standard deviations for each component in the second margin
//
// nknots[2]:        Numbers of knots in each margin
//
// intcpt[2]:        Overall intercept for each margin
// tau[2]:           Overall standard deviation for each margin
//
// ======================================================================================
//
// REMARK:  EY    = mean_g = sum(w[j]*mu[j])
//         var(Y) = sum{w[j]*(sigma[j]^2 + mu[j]^2)} - mean_g^2
//                = sum{w[j]*(sigma[j]^2 + mu[j]^2 - mean_g^2)}
//                = sum{w[j]*(sigma[j]^2 + (mu[j] - mean_g)^2)}
//
// REMARK: univariately, not standardized version:
//         g(y) = sum[j](w[j]*(sigma[j]*tau)^(-1)*phi((y-intcpt-tau*mu[j])/(sigma[j]*tau)))
//         and similarly bivariately
//
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
               const double *intcpt,      const double *tau,        const int *logw)
{
  static int ix0, ix1, k0, k1;
  static double *wP, *wP0, *wP1;
  static double *knots_tauP0, *knots_tauP1;
  static double *sigma_tauP0, *sigma_tauP1;
  static double *averP, *averP0, *averP1;
  static double *valP, *valP0, *valP1;  
  static const int *indP;
  static const double *knotsP0, *knotsP1;
  static const double *sigmaP0, *sigmaP1;
  static const double *gridP0, *gridP1;
  static double val_gk, val_gk0, val_gk1, knot_mean;
  static double mean_g0, mean_g1;
  static double sd_g0, sd_g1;
  static double sumexpa;

  /*** Compute exp(log-weights) if necessary           ***/
  /*** (do not divide them by sumexpa)                 ***/ 
  /*** and compute weights for marginal distributions  ***/
  /*** =============================================== ***/
  wP0 = weight0;
  for (k0 = 0; k0 < nknots[0]; k0++){
    *wP0 = 0.0;
    wP0++;
  }
  sumexpa = 0.0;

  indP = ind_w_effect;
  if (*logw){
    wP  = weight;
    wP1 = weight1; 
    for (k1 = 0; k1 < nknots[1]; k1++){
      *wP1 = 0.0;
      wP0 = weight0;
      for (k0 = 0; k0 < nknots[0]; k0++){
        *wP = *indP ? exp_AK(*wP) : 0.0;
        
        sumexpa += *wP;     
        *wP0    += *wP;
        *wP1    += *wP;

        indP++;
        wP++;
        wP0++;
      }
      wP1++;
    }
  }
  else{              /** only sum-up the weights **/
    wP  = weight;
    wP1 = weight1; 
    for (k1 = 0; k1 < nknots[1]; k1++){
      *wP1 = 0.0;
      wP0 = weight0;
      for (k0 = 0; k0 < nknots[0]; k0++){
        if (!(*indP)) *wP = 0.0;

        sumexpa += *wP;     
        *wP0    += *wP;
        *wP1    += *wP;

        indP++;
        wP++;
        wP0++;
      }
      wP1++;
    }
  }

  //Rprintf("*** Weights in margin0 (multiplied by %g):\n", sumexpa);
  //AK_BLAS_LAPACK::printArray(weight0, nknots[0]);
  //Rprintf("*** Weights in margin1 (multiplied by %g):\n", sumexpa);
  //AK_BLAS_LAPACK::printArray(weight1, nknots[1]);
  
  if (*standard){
    /*** Compute the mean of the G-spline in margin 0 (without added intercept) ***/
    /*** ====================================================================== ***/
    mean_g0 = 0.0;
    wP0     = weight0;
    knotsP0 = knots0;
    for (k0 = 0; k0 < nknots[0]; k0++){
      mean_g0 += (*wP0)*(*knotsP0);                                         // += w[k]*mu[k]
      wP0++;
      knotsP0++;
    }
    mean_g0 /= sumexpa;

    /*** Compute the mean of the G-spline in margin 1 (without added intercept) ***/
    /*** ====================================================================== ***/
    mean_g1 = 0.0;
    wP1     = weight1;
    knotsP1 = knots1;
    for (k1 = 0; k1 < nknots[1]; k1++){
      mean_g1 += (*wP1)*(*knotsP1);                                         // += w[k]*mu[k]
      wP1++;
      knotsP1++;
    }
    mean_g1 /= sumexpa;

    /*** Compute the variance of the G-spline in margin 0 (not multiplied by the overall variance) ***/
    /*** ========================================================================================= ***/
    sd_g0   = 0.0;
    wP0     = weight0;
    knotsP0 = knots0;
    sigmaP0 = sigma0;
    for (k0 = 0; k0 < nknots[0]; k0++){
      knot_mean = *knotsP0 - mean_g0;
      sd_g0     += *wP0 * (*sigmaP0 * *sigmaP0 + knot_mean * knot_mean);
      wP0++;
      knotsP0++;
      sigmaP0++;
    }
    sd_g0 /= sumexpa;

    /*** Compute the variance of the G-spline in margin 1 (not multiplied by the overall variance) ***/
    /*** ========================================================================================= ***/
    sd_g1   = 0.0;
    wP1     = weight1;
    knotsP1 = knots1;
    sigmaP1 = sigma1;
    for (k1 = 0; k1 < nknots[1]; k1++){
      knot_mean = *knotsP1 - mean_g1;
      sd_g1     += *wP1 * (*sigmaP1 * *sigmaP1 + knot_mean * knot_mean);
      wP1++;
      knotsP1++;
      sigmaP1++;
    }
    sd_g1 /= sumexpa;

    if (sd_g0 <= 0 || sd_g1 <= 0){
      throw returnR("Error in summary_BiGspline.cpp: eval_BiGspline. Variance of the G-spline is not positive", 1);
    }
    sd_g0 = sqrt(sd_g0);
    sd_g1 = sqrt(sd_g1);

    //sd_g0 *= tau[0];                // Taking into account the scale tau
    //sd_g1 *= tau[1];                // Taking into account the scale tau
    //mean_g0 *= tau[0];              // Taking into account the scale tau
    //mean_g0 += intcpt[0];           // Taking into account the intercept
    //mean_g1 *= tau[1];              // Taking into account the scale tau
    //mean_g1 += intcpt[1];           // Taking into account the intercept

    /*** Compute the value of the G-spline in the grid ***/
    /*** ============================================= ***/
    valP   = value;
    averP  = average;

    gridP1 = grid1;
    for (ix1 = 0; ix1 < ngrid[1]; ix1++){

      gridP0 = grid0;
      for (ix0 = 0; ix0 < ngrid[0]; ix0++){
        *valP = 0.0;
        wP    = weight;
        indP  = ind_w_effect;

        knotsP1 = knots1;
        sigmaP1 = sigma1;
        for (k1 = 0; k1 < nknots[1]; k1++){

          knotsP0 = knots0;
          sigmaP0 = sigma0;
          for (k0 = 0; k0 < nknots[0]; k0++){
            if (*indP){
              val_gk0 = (sd_g0*(*gridP0) + mean_g0 - (*knotsP0))/(*sigmaP0);
              val_gk0 = -0.5 * val_gk0 * val_gk0;

              val_gk1 = (sd_g1*(*gridP1) + mean_g1 - (*knotsP1))/(*sigmaP1);
              val_gk1 = -0.5 * val_gk1 * val_gk1;

              val_gk  = val_gk0 + val_gk1;
              val_gk  = val_gk < _AK_EMIN ? 0.0 : (*wP)*exp(val_gk)/(*sigmaP0 * *sigmaP1);
              *valP  += val_gk;
	    }
 
            knotsP0++;
            sigmaP0++;
            wP++;
            indP++;
          }
          
          knotsP1++;
          sigmaP1++;
        }

        *valP /= sumexpa;
        *valP *= sd_g0 * sd_g1;        
        *valP *= _AK_invTWO_PI;
        *averP += *valP;

        valP++;
        averP++;
        gridP0++;
      }
      gridP1++;
    }

    /*** Compute the value of the marginal G-spline in margin 0 ***/
    /*** ====================================================== ***/
    valP0  = value0;
    averP0 = average0;
    gridP0 = grid0;
    for (ix0 = 0; ix0 < ngrid[0]; ix0++){
      *valP0  = 0.0;
      wP0     = weight0;
      knotsP0 = knots0;
      sigmaP0 = sigma0;
      for (k0 = 0; k0 < nknots[0]; k0++){
        val_gk0 = (sd_g0*(*gridP0) + mean_g0 - (*knotsP0))/(*sigmaP0);
        val_gk0 = -0.5 * val_gk0 * val_gk0;

        val_gk0  = val_gk0 < _AK_EMIN ? 0.0 : (*wP0)*exp(val_gk0)/(*sigmaP0);
        *valP0  += val_gk0;

        knotsP0++;
        sigmaP0++;
        wP0++;
      }

      *valP0 /= sumexpa;
      *valP0 *= sd_g0;        
      *valP0 *= _AK_invSQRT_TWO_PI;
      *averP0 += *valP0;

      valP0++;
      averP0++;
      gridP0++;
    }

    /*** Compute the value of the marginal G-spline in margin 1 ***/
    /*** ====================================================== ***/
    valP1  = value1;
    averP1 = average1;
    gridP1 = grid1;
    for (ix1 = 0; ix1 < ngrid[1]; ix1++){
      *valP1  = 0.0;
      wP1     = weight1;
      knotsP1 = knots1;
      sigmaP1 = sigma1;
      for (k1 = 0; k1 < nknots[1]; k1++){
        val_gk1 = (sd_g1*(*gridP1) + mean_g1 - (*knotsP1))/(*sigmaP1);
        val_gk1 = -0.5 * val_gk1 * val_gk1;

        val_gk1  = val_gk1 < _AK_EMIN ? 0.0 : (*wP1)*exp(val_gk1)/(*sigmaP1);
        *valP1  += val_gk1;

        knotsP1++;
        sigmaP1++;
        wP1++;
      }

      *valP1 /= sumexpa;
      *valP1 *= sd_g1;        
      *valP1 *= _AK_invSQRT_TWO_PI;
      *averP1 += *valP1;

      valP1++;
      averP1++;
      gridP1++;
    }
  }
  else{    /***   !(*standard)   ***/
    /*** Compute knots_tau and sigma_tau in margin 0 ***/  
    /*** =========================================== ***/
    knots_tauP0 = knots_tau0;
    sigma_tauP0 = sigma_tau0;
    knotsP0     = knots0;
    sigmaP0     = sigma0;
    for (k0 = 0; k0 < nknots[0]; k0++){
      *knots_tauP0 = *knotsP0 * tau[0];
      *sigma_tauP0 = *sigmaP0 * tau[0];
      knots_tauP0++;
      sigma_tauP0++;
      knotsP0++;
      sigmaP0++;
    }   

    /*** Compute knots_tau and sigma_tau in margin 1 ***/  
    /*** =========================================== ***/
    knots_tauP1 = knots_tau1;
    sigma_tauP1 = sigma_tau1;
    knotsP1     = knots1;
    sigmaP1     = sigma1;
    for (k1 = 0; k1 < nknots[1]; k1++){
      *knots_tauP1 = *knotsP1 * tau[1];
      *sigma_tauP1 = *sigmaP1 * tau[1];
      knots_tauP1++;
      sigma_tauP1++;
      knotsP1++;
      sigmaP1++;
    }   

    /*** Compute the value of the G-spline in the grid ***/
    /*** ============================================= ***/
    valP   = value;
    averP  = average;

    gridP1 = grid1;
    for (ix1 = 0; ix1 < ngrid[1]; ix1++){

      gridP0 = grid0;
      for (ix0 = 0; ix0 < ngrid[0]; ix0++){
        *valP = 0.0;
        wP    = weight;
        indP  = ind_w_effect;

        knots_tauP1 = knots_tau1;
        sigma_tauP1 = sigma_tau1;
        for (k1 = 0; k1 < nknots[1]; k1++){

          knots_tauP0 = knots_tau0;
          sigma_tauP0 = sigma_tau0;
          for (k0 = 0; k0 < nknots[0]; k0++){
            if (*indP){
              val_gk0 = ((*gridP0) - intcpt[0] - (*knots_tauP0))/(*sigma_tauP0);
              val_gk0 = -0.5 * val_gk0 * val_gk0;

              val_gk1 = ((*gridP1) - intcpt[1] - (*knots_tauP1))/(*sigma_tauP1);
              val_gk1 = -0.5 * val_gk1 * val_gk1;

              val_gk  = val_gk0 + val_gk1;
              val_gk  = val_gk < _AK_EMIN ? 0.0 : (*wP)*exp(val_gk)/(*sigma_tauP0 * *sigma_tauP1);
              *valP  += val_gk;
            }
 
            knots_tauP0++;
            sigma_tauP0++;
            wP++;
            indP++;
          }
          
          knots_tauP1++;
          sigma_tauP1++;
        }

        *valP /= sumexpa;
        *valP *= _AK_invTWO_PI;
        *averP += *valP;

        valP++;
        averP++;
        gridP0++;
      }
      gridP1++;
    }

    /*** Compute the value of the marginal G-spline in margin 0 ***/
    /*** ====================================================== ***/
    valP0  = value0;
    averP0 = average0;
    gridP0 = grid0;
    for (ix0 = 0; ix0 < ngrid[0]; ix0++){
      *valP0  = 0.0;
      wP0     = weight0;
      knots_tauP0 = knots_tau0;
      sigma_tauP0 = sigma_tau0;
      for (k0 = 0; k0 < nknots[0]; k0++){
        val_gk0 = ((*gridP0) - intcpt[0] - (*knots_tauP0))/(*sigma_tauP0);
        val_gk0 = -0.5 * val_gk0 * val_gk0;

        val_gk0  = val_gk0 < _AK_EMIN ? 0.0 : (*wP0)*exp(val_gk0)/(*sigma_tauP0);
        *valP0  += val_gk0;

        knots_tauP0++;
        sigma_tauP0++;
        wP0++;
      }

      *valP0 /= sumexpa;
      *valP0 *= _AK_invSQRT_TWO_PI;
      *averP0 += *valP0;

      valP0++;
      averP0++;
      gridP0++;
    }

    /*** Compute the value of the marginal G-spline in margin 1 ***/
    /*** ====================================================== ***/
    valP1  = value1;
    averP1 = average1;
    gridP1 = grid1;
    for (ix1 = 0; ix1 < ngrid[1]; ix1++){
      *valP1  = 0.0;
      wP1     = weight1;
      knots_tauP1 = knots_tau1;
      sigma_tauP1 = sigma_tau1;
      for (k1 = 0; k1 < nknots[1]; k1++){
        val_gk1 = ((*gridP1) - intcpt[1] - (*knots_tauP1))/(*sigma_tauP1);
        val_gk1 = -0.5 * val_gk1 * val_gk1;

        val_gk1  = val_gk1 < _AK_EMIN ? 0.0 : (*wP1)*exp(val_gk1)/(*sigma_tauP1);
        *valP1  += val_gk1;

        knots_tauP1++;
        sigma_tauP1++;
        wP1++;
      }

      *valP1 /= sumexpa;
      *valP1 *= _AK_invSQRT_TWO_PI;
      *averP1 += *valP1;

      valP1++;
      averP1++;
      gridP1++;
    }
  }

  return;
}

}   /*** end of extern "C" ***/


/*** readWeightsFromFilesBi:  Read (log-)weights from the file (bivariate G-spline)                            ***/
/***                                                                                                           ***/
/***    It is assumed that the file with weights contains only the weights on "non-zero" components.           ***/
/***    Number and indeces of non-zero components are read from another file.                                  ***/
/***                                                                                                           ***/
/*** ========================================================================================================= ***/
//
// k_effect:                     OUTPUT: Read number of "non-zero" components
// ind_w_effect[total_length]:   OUTPUT: Derived 0/1 indicators of "non-zero" components
//                               ASSUMPTION: Indeces of "non-zero" components in indfile are on the scale
//                                           0,..., total_length-1
//
// weight[total_length]:         OUTPUT: Read (log-)weights
//                                       Weights of non-zero components are filled by zero
//
// skip:   Number of rows that should be skipped before reading real data
// row:    How many rows of the data have already been read before the function call?
//         (only used to generate possible error messages)
//
void
readWeightsFromFilesBi(int *k_effect,           int *ind_w_effect,           double *weight,  
                       const int &skip,         const int &row,              const int *total_length,
                       std::ifstream &indfile,  const std::string& sindpath,
                       std::ifstream &wfile,    const std::string& swpath)
{
  try{
    static int j, k, ihelp, index;
    static char ch;
    static std::string errmes;
    static double *wP;
    static int *indP;
    int iihelp;

    /**  Skip rows that are to be skipped  **/
    for (j = 0; j < skip; j++){
      wfile.get(ch);        
      while (ch != '\n') wfile.get(ch);

      indfile.get(ch);
      while (ch != '\n') indfile.get(ch);
    }


    /** Read number and indeces of non-zero weights **/
    if (indfile.eof()){
      ihelp = row + 1;
      errmes = std::string("Error in summary_BiGspline.cpp: readWeightsFromFilesBi(): Reached end of file ") + sindpath + " before "
               + char(ihelp) + std::string(" sets of G-spline indeces were read.");
      throw returnR(errmes, 99);
    }

    indfile >> iihelp;
    if (!iihelp){
      REprintf("Iteration %d,  G-spline with 0 non-zero components.\n", row+1);
      throw returnR("Error in summary_BiGspline.cpp: readWeightsFromFilesBi()", 99);
    }
    *k_effect = iihelp;

    indP = ind_w_effect;
    k    = 0;
    for (j = 0; j < *k_effect; j++){
      indfile >> index;
      while (k < index){
        *indP = 0;
        indP++;      
        k++;
      }
      *indP = 1;
      indP++;
      k++;
    }
    while (k < *total_length){
      *indP = 0;
      indP++;
      k++;
    }
    indfile.get(ch);                 
    while (ch != '\n') indfile.get(ch);


    /** Read G-spline weights **/
    if (wfile.eof()){
      ihelp = row + 1;
      errmes = std::string("Error in summary_BiGspline1.cpp: readWeightsFromFilesBi(): Reached end of file ") + swpath + " before "
               + char(ihelp) + std::string(" sets of G-spline weights were read.");
      throw returnR(errmes, 99);
    }

    wP   = weight;
    indP = ind_w_effect; 
    for (k = 0; k < *total_length; k++){
      if (*indP) wfile >> *wP;
      else       *wP = 0;
      wP++;
      indP++;
    }
    wfile.get(ch);                 
    while (ch != '\n') wfile.get(ch);

    return;
  }
  catch(returnR){
    throw;
  }  
}


}   /*** end of the namespace summary_BiGsplineA ***/
