/*** Random2.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  24/10/2006
//
//        discreteSampler2:   24/10/2006
//        findUniformIndex:   24/10/2006
//               findIndex:   24/10/2006
//            rltruncGamma:   27/10/2006
//               rGspline1:   29/01/2007
//              rGspline1R:   29/01/2007
//              rBiGspline:   12/04/2007
//             rBiGsplineR:   12/04/2007
//
// PURPOSE: Random numbers generator (several random number generation)
//
//   Partially taken from 'random.cpp' of 'bayesSurv' package
//
/* ********************************************************************************* */

#include "Random2.h"

namespace Random {

// ********** rltruncGamma ************
// Sample from left-truncated gamma distribution
//
// PARAMETERS:
//   shape .......... shape of the gamma distribution
//   rate ........... rate (= 1/scale) of the gamma distribution
//   minx ........... a value > 0 indicating where the gamma distribution is to be truncated
//   n .............. number of random numbers that should be returned in the vector x
//   
void
rltruncGamma(double* x,  
             const double* shape,  const double* rate,   const double* minx,
             const int* n)
{
  int i;
  double u, *xP;  

  /* Rprintf("minx = %e, shape = %e, rate = %e\n", *minx, *shape, *rate);  */
  double scale = 1/(*rate);
  double Flower = pgamma(*minx, *shape, scale, 1, 0);
  if (Flower >= 1 - _AK_NORM_ZERO){   // truncation time irrealistic large => sampled values = minx
    xP = x;
    for (i = 0; i < *n; i++){
      *xP = *minx;
      xP++;
    }
  }
  else{
    if (Flower <= _AK_NORM_ZERO){     // truncation time = 0, sample from gamma distribution
      xP = x;
      for (i = 0; i < *n; i++){
        *xP = rgamma(*shape, scale);
        xP++;
      }
    }
    else{
      xP = x;
      for (i = 0; i < *n; i++){
        u = runif(0, 1) * (1 - Flower) + Flower;
        *xP = qgamma(u, *shape, scale, 1, 0);
        xP++;
      }
    }
  }

  return;
}


// ********** discreteSampler2 **********
//
// Just small modification of 'discreteSampler' from 'random.cpp' of 'bayesSurv' package
//   * it does not check which proportions are zero
//   * to be used primarily by updateAlloc
//
// Sample from a discrete distribution on 0, 1, ..., *kP - 1
//   when either cumulative proportions are given or proportions are given,
//   i.e. if (cumul){
//          P(Y = 0) \propto propA[0]
//          P(Y = j) \propto propA[j] - propA[j-1], j = 1, ..., *kP - 1
//        }
//        else{
//          P(Y = j) \propto  propA[j], j = 0, ..., *kp - 1
//        }
//
// ASSUMPTIONS:
//   * if (cumul)
//       (i)  \exist j propA[j] > 0 
//       (ii) propA[0] <= ... <= propA[*kP - 1]
//   * if (!cumul)
//       (ii) \exist j propA[j] > 0  
//       (i)  propA[j] => 0 for all j  
//
// I do not check these assumptions in the C++ code!!!
//
// PARAMETERS:
//
// propA ................. array of either proportions or cumulative proportions
// kP .................... length of the array propA
// nP .................... number of random variates to be sampled 
// cumul ................. logical, true if the array propA gives cumulative proportions
//
// RETURN:
//
// sampledj .............. array of sampled values
//
void
discreteSampler2(int* sampledj,  double* propA,  const int* kP,  const int* nP,  const int* cumul)
{ 
  int i, j;
  double u;
  double propMax;

  int *jP;
  double *cumwP;

  if (*kP <= 1){
    jP = sampledj;
    for (i = 0; i < *nP; i++){
      *jP = 0;
      jP++;
    }
    return;
  }

  // Compute cumulative proportions if necessary
  if (!(*cumul)){
    cumwP = propA + 1;
    for (j = 1; j < *kP; j++){
      *cumwP += *(cumwP - 1);
      cumwP++; 
    }
  }

  propMax = propA[(*kP) - 1];
  jP = sampledj;
  for (i = 0; i < *nP; i++){
    u = runif(0, propMax);
    *jP = Random::findIndex(u, 0, (*kP) - 1, propA);
    jP++;
  }

  return;  
}   // end of function discreteSampler2


// ********** findUniformIndex **********
// Find the index of the smallest element of ValuesA = (1/k, 2/k, ..., k/k)
//  which is higher than u 
// The search is restricted to ValuesA[startInd], ..., ValuesA[endInd]
//
// This function is to be used by 'discreteUniformSampler' routine.
// That's why it relies on the following assumptions:
//   * startInd < endInd
//   * 0 <= u <= 1
//
// I do not check these assumptions since it would be wasting of CP time.
//
int 
findUniformIndex(const double u,  const int startInd,  const int endInd,  const int k)
{
  if (startInd == endInd - 1){
    if (u <= double(startInd+1)/double(k)) return startInd;
    else                                   return endInd;
  } 
  else{
    int midInd = int(ceil(0.5 * (startInd + endInd)));
    int Index;
    if (u <= double(midInd+1)/double(k)) Index = Random::findUniformIndex(u, startInd, midInd, k);
    else                                 Index = Random::findUniformIndex(u, midInd, endInd, k); 
    return Index;
  }
}   // end of function findUniformIndex



// ********** findIndex **********
// Find the index of the smallest element of ValuesA
//  which is higher than u 
// The search is restricted to ValuesA[startInd], ..., ValuesA[endInd]
//
// This function is to be used by 'discreteSampler' routine.
// That's why it relies on the following assumptions:
//   * startInd < endInd
//   * ValuesA[0] > 0
//   * ValuesA[startInd] < ... < ValuesA[endInd]
//   * 0 <= u <= valuesA[endInd]
//
// I do not check these assumptions since it would be wasting of CP time.
//
int 
findIndex(const double u,  const int startInd,  const int endInd,  const double* ValuesA)
{
  if (startInd == endInd - 1){
    if (u <= ValuesA[startInd]) return startInd;
    else                        return endInd;
  } 
  else{
    int midInd = int(ceil(0.5 * (startInd + endInd)));
    int Index;
    if (u <= ValuesA[midInd]) Index = Random::findIndex(u, startInd, midInd, ValuesA);
    else                      Index = Random::findIndex(u, midInd, endInd, ValuesA); 
    return Index;
  }
}   // end of function findIndex


/*** rGspline1:   Random number generation from the UNIVARIATE G-spline     ***/
/*** rGspline1R:  Interface to R                                            ***/
/***                                                                        ***/
/*** ====================================================================== ***/
//
// x[nx]:            Generated random number(s)
// weight[nknots]:   INPUT:   (Log-)weights (up to a proportionality constant)
//                   OUTPUT:  Cumulative weights (up to a proportionality constant)
// nx:               Number of points to generate
// knots[nknots]:    Knots (mixture means)
// sigma[nknots]:    Standard deviations of mixture components
// nknots:           Number of knots
// intcpt:           Intercept
// tau:              Scale
// logw:             If <> 0 then it is assumed that 'weight' contains log-weights
// 
//
void
rGspline1(double *x,             double *weight, 
          const int *nx,
          const double *knots,   const double *sigma,  const int *nknots,
          const double *intcpt,  const double *tau,    const int *logw)
{
  static int ix, k, sampled_k;
  static double *wP, *xP;
  static double sumw, u;

  /*** Compute cumulative weights ***/
  wP = weight;  
  if (*logw){
    *wP = exp_AK(*wP);
    for (k = 1; k < *nknots; k++){
      wP++;
      *wP = exp_AK(*wP);
      *wP += *(wP-1);
    }
  }
  else{
    for (k = 1; k < *nknots; k++){
      wP++;
      *wP += *(wP-1);
    }
  }
  sumw = *wP;
  
  /*** Loop over the sampled values ***/
  xP = x;
  for (ix = 0; ix < *nx; ix++){

    /*** Sample the mixture component ***/
    u = runif(0, sumw);
    sampled_k = Random::findIndex(u, 0, (*nknots) - 1, weight);

    /*** Sample the value from the mixture component ***/
    *xP = norm_rand();
    *xP *= sigma[sampled_k];
    *xP += knots[sampled_k];
    *xP *= *tau;
    *xP += *intcpt;

    xP++;
  }

  return;
}  /*** end of rGspline1() ***/


extern "C"{

void
rGspline1R(double *x,             double *weight, 
           const int *nx,
           const double *knots,   const double *sigma,  const int *nknots,
           const double *intcpt,  const double *tau,    const int *logw)
{
  try{
    GetRNGstate();

    Random::rGspline1(x, weight, nx, knots, sigma, nknots, intcpt, tau, logw);

    PutRNGstate();
    return;
  }
  catch(returnR rr){
    PutRNGstate();
    return;
  }
}

}  /*** end of extern "C" ***/

/*** rBiGspline:   Random number generation from the BIVARIATE G-spline      ***/
/*** rBiGsplineR:  Interface to R                                            ***/
/***                                                                         ***/
/*** ======================================================================= ***/
//
// x[2*nx]:                Generated random number(s)
// weight[total_length]:   INPUT:   (Log-)weights (up to a proportionality constant)
//                                   Places that correspond to zero weights are filled by zeros and ignored
//                         OUTPUT:  Cumulative weights (up to a proportionality constant)
//                                  for non-zero components on first k_effect places
//
// ind_w_effect[total_length]:   Used only if (is_indfile <> 0)
//                               in that case INPUT:   0/1 indicators of non-zero components
//                                            OUTPUT:  indeces of non-zero components
//                                                     (on scale 0,...,nknots[0]*nknots[1]-1)
//                               
// nx:                   Number of points to generate
//
// knots0[nknots[0]]:    Knots (mixture means) for margin 0
// knots1[nknots[1]]:    Knots (mixture means) for margin 1
// sigma0[nknots[0]]:    Standard deviations of mixture components for margin 0
// sigma1[nknots[1]]:    Standard deviations of mixture components for margin 1
//
// nknots[2]:            Numbers of knots in each margin
// total_length[1]:      = nknots[0] * nknots[1]
//
// intcpt[2]:            Intercepts for each margin
// tau[2]:               Scales for each margin
//
// logw:                 If <> 0 then it is assumed that 'weight' contains log-weights
// is_indfile:           If = 0 then it is assumed that k_effect = total_length and 
//                       ind_w_effect is ignored
// 
void
rBiGspline(double *x,              double *weight,             int *ind_w_effect,
           const int *nx,
           const double *knots0,   const double *knots1,
           const double *sigma0,   const double *sigma1,       const int *nknots,   const int *total_length,
           const double *intcpt,   const double *tau,      
           const int *logw,        const int *is_indfile)
{
  static int ix, k_effect, k, istop, sampled_k, k0, k1;
  static double *wP, *wP1, *xP;
  static int *indP, *indP1;
  static double sumw, u;

  /*** Compute cumulative weights and if needed, re-calculate ind_w_effect ***/
  wP = weight;
  if (*is_indfile){
    k_effect = 0;

    /** Iterate to the first non-zero weight **/
    wP1   = weight;
    indP1 = ind_w_effect;

    indP  = ind_w_effect;
    istop = 0;
    while (*indP1 == 0){
      wP1++;
      indP1++;
      istop++;
      if (istop == *total_length) throw returnR("Error in Random2.cpp: rBiGspline(). All G-spline components have zero weights", 99);
    }

    if (*logw){
      *wP      = exp_AK(*wP1);
      *indP    = istop;
      k_effect = 1;
      for (k = istop + 1; k < *total_length; k++){
        wP1++;
        indP1++;
        if (*indP1){
          wP++;
          indP++;
          k_effect++;
          *wP = exp_AK(*wP1);
          *wP += *(wP-1);
          *indP = k;
        }
      }
    }
    else{
      *wP      = *wP1;
      *indP    = istop;
      k_effect = 1;
      for (k = istop + 1; k < *total_length; k++){
        wP1++;
        indP1++;
        if (*indP1){
          wP++;
          indP++;
          k_effect++;
          *wP = *wP1;
          *wP += *(wP-1);
          *indP = k;
        }
      }
    }
    sumw = *wP;
    //Rprintf("\nk_effect=%d,  sumw=%g\n", k_effect, sumw);

    /*** Loop over the sampled values ***/
    xP = x;
    for (ix = 0; ix < *nx; ix++){

      /*** Sample the mixture component ***/
      u = runif(0, sumw);
      sampled_k = Random::findIndex(u, 0, k_effect - 1, weight);
      sampled_k = ind_w_effect[sampled_k];

      /*** Compute double index of the mixture component  (on the scale 0, ..., nknots) ***/
      k0 = sampled_k % nknots[0];
      k1 = sampled_k / nknots[0];
      
      /*** Sample the value from the mixture component ***/
      *xP = norm_rand();
      *xP *= sigma0[k0];
      *xP += knots0[k0];
      *xP *= tau[0];
      *xP += intcpt[0];
      xP++;

      *xP = norm_rand();
      *xP *= sigma1[k1];
      *xP += knots1[k1];
      *xP *= tau[1];
      *xP += intcpt[1];
      xP++;
    }
  }
  else{  /** ignore ind_w_effect, k_effect = total_length **/
    k_effect = *total_length;
    if (*logw){
      *wP = exp_AK(*wP);
      for (k = 1; k < *total_length; k++){
        wP++;
        *wP = exp_AK(*wP);
        *wP += *(wP-1);
      }
    }
    else{
      for (k = 1; k < *total_length; k++){
        wP++;
        *wP += *(wP-1);
      }
    }
    sumw = *wP;
    //Rprintf("\nk_effect=%d,  sumw=%g\n", k_effect, sumw);

    /*** Loop over the sampled values ***/
    xP = x;
    for (ix = 0; ix < *nx; ix++){

      /*** Sample the mixture component ***/
      u = runif(0, sumw);
      sampled_k = Random::findIndex(u, 0, k_effect - 1, weight);

      /*** Compute double index of the mixture component  (on the scale 0, ..., nknots) ***/
      k0 = sampled_k % nknots[0];
      k1 = sampled_k / nknots[0];
      
      /*** Sample the value from the mixture component ***/
      *xP = norm_rand();
      *xP *= sigma0[k0];
      *xP += knots0[k0];
      *xP *= tau[0];
      *xP += intcpt[0];
      xP++;

      *xP = norm_rand();
      *xP *= sigma1[k1];
      *xP += knots1[k1];
      *xP *= tau[1];
      *xP += intcpt[1];
      xP++;
    }
  }

  return;
}  /*** end of rGsplineBi() ***/

extern "C"{

void
rBiGsplineR(double *x,              double *weight,             int *ind_w_effect,
            const int *nx,
            const double *knots0,   const double *knots1,
            const double *sigma0,   const double *sigma1,       const int *nknots,   const int *total_length,
            const double *intcpt,   const double *tau,      
            const int *logw,        const int *is_indfile)
{
  try{
    GetRNGstate();

    Random::rBiGspline(x, weight, ind_w_effect, nx, knots0, knots1, sigma0, sigma1, nknots, total_length, intcpt, tau, logw, is_indfile);

    PutRNGstate();
    return;
  }
  catch(returnR rr){
    PutRNGstate();
    return;
  }
}

}  /*** end of extern "C" ***/

}  /*** end of namespace Random ***/
