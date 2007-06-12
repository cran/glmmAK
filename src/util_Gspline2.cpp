/*** util_Gspline2.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  20/10/2006
//
//                               max_poster_prob:  20/10/2006
//      array2alloc, alloc2array, alloc2mixtureN:  23/10/2006
//                                  updateAlloc1:  29/10/2006
//
// PURPOSE: G-spline, several smaller functions related to it
//
/* ********************************************************************************* */

#include "util_Gspline2.h"

namespace util_Gspline2{

extern "C"{

/***** max_poster_prob ********************************************************************************************************************/
/*                                                                                                                                        */
/* PURPOSE: For each marginal data point, determine the mixture component having the highest probability for that point                   */
/*                                                                                                                                        */
/* b[j,i] | k ~ alpha[j] + tau[j] N(mu[j,k], sigma2[j])                                                                                   */
/*            ~ N(alpha[j] + tau[j]*mu[j,k], tau2[j]*sigma2[j])                                                                           */
/* p(b[j,i] | k) = phi((b[j,i] - alpha[j])/tau[j] | mu[j,k], sigma2[j])                                                                   */
/*               = phi((b[j,i] - alpha[j] - tau[j]*mu[j,k])/(tau[j]*sigma[j]))                                                            */
/*                                                                                                                                        */
/******************************************************************************************************************************************/
//
// Alloc[Dim x nData]:
//    OUTPUT: allocations (on the scale -K,...,K) for the Data to the mixture components having the highest value of the density
//            that is, Alloc[j,i] = argmax_k p(b[j,i] | k), k=-K[j],...,K[j]
//
// BasisStdDev[Dim]:
//    INPUT:   basis standard deviations for each margin sigma[j], j=1,...,Dim
//    OUTPUT:  sigma[j]*tau[j], j=1,...,Dim
//
// Knots[sum_length]:
//    INPUT:  knots for each margin mu[k,j], k=-K[j],...,K[j], j=1,...,Dim
//    OUTPUT: Knots[k,j] = alpha[j] + tau[j]*mu[k,j], k=-K[j],...,K[j], j=1,...,Dim
//
// K[Dim]
//    number of knots on each side of the reference knot for each margin
//    length[i] = 2*K[i] + 1
//    sum_length = sum(length[i])
//
// Data[Dim x nData]:
//    data points b[j,i], j=1,...,Dim; i=1,...,nData
//
// Dim[1]:
//    dimension of each data point
//
// nData[1]:
//    number of each data point
//
// Intcpt[Dim]:
//    intercepts for each margin alpha[j], j=1,...,Dim
//
// StdDev[Dim]:
//    standard deviations for each margin tau[j], j=1,,,.Dim
//
void
max_poster_prob(int *Alloc,               double *BasisStdDev,   double *Knots,      const int *K,
                const double *Data,       const int *Dim,        const int *nData,
                const double *Intcpt,     const double *StdDev)
{
  int i, j, k;

  const double *intcptP, *stdDevP, *dataP; 
  const int *KP;
  double *basisStdDevP, *knotP;  
  int *allocP;
  double log_pb_max, log_pb_NEW;

  /*** Compute tau*sigma and alpha + tau*mu ***/
  intcptP = Intcpt;
  stdDevP = StdDev;
  basisStdDevP = BasisStdDev;
  KP = K;   
  knotP = Knots;
  for (j = 0; j < *Dim; j++){
    *basisStdDevP *= *stdDevP;
    for (k = -(*KP); k <= *KP; k++){
      *knotP *= *stdDevP;
      *knotP += *intcptP;
      knotP++;
    }
    basisStdDevP++;
    KP++;
    intcptP++;
    stdDevP++;
  }

  /*** Determine allocations that maximize the density of the component ***/
  dataP = Data;
  allocP = Alloc;
  for (i = 0; i < *nData; i++){
    basisStdDevP = BasisStdDev;
    KP = K;
    knotP = Knots;
    for (j = 0; j < *Dim; j++){
      *allocP = -(*KP);
      log_pb_max = dnorm(*dataP, *knotP, *basisStdDevP, 1);
      knotP++;
      for (k = -(*KP)+1; k <= *KP; k++){      
        log_pb_NEW = dnorm(*dataP, *knotP, *basisStdDevP, 1);
        if (log_pb_NEW >= log_pb_max){
          log_pb_max = log_pb_NEW;
          *allocP = k;
        }
        knotP++;
      }
      basisStdDevP++;
      KP++;
      dataP++;
      allocP++;
    }
  }

  return;
}

}    /* end of extern "C" */


/***** array2alloc ****************************************************************************************************************/
/*                                                                                                                                */
/* PURPOSE: Create a matrix of allocations                                                                                        */
/*          and perform input checks                                                                                              */
/*                                                                                                                                */
/**********************************************************************************************************************************/
//
// allocM[dim x nObs]:  OUTPUT: matrix dim x nObs with allocations on scale -K,...,K
// allocA[dim * nObs]:  INPUT: array with allocations on scale -K,...,K, 
//                             = matrix dim x nObs stored in column major order, each column = allocations for each observation
// K[dim]:              numbers of knots on each side of the reference knot for each margin
// dim:                 dimension
// nObs:                number of observations
//
void
array2alloc(MatrixRect<int> *allocM,  const int *allocA,  const int *K,  const int &dim,  const int &nObs)
{
  int i, j;
  int *aMP;
  const int *aAP, *KP;

  if (!dim || !nObs){
    *allocM = MatrixRect<int>(0, 0);
    return;
  }

  *allocM = MatrixRect<int>(dim, nObs);
  aMP = allocM->a();
  aAP = allocA;
  for (i = 0; i < nObs; i++){
    KP = K;
    for (j = 0; j < dim; j++){
      if (*aAP < -(*KP) || *aAP > *KP) throw returnR("Error in util_Gspline2.cpp: array2alloc. Allocation out of the range.", 1);
      *aMP = *aAP;                          /*** If shift from scale -K,...,K to scale 0,...,2*K needed:  *aMP = *aAP + (*KP); ***/
      aAP++;
      aMP++;
      KP++;
    }
  }

  return;
}


/***** alloc2array ****************************************************************************************************************/
/*                                                                                                                                */
/* PURPOSE: Copy a matrix of allocations                                                                                          */
/*                                                                                                                                */
/**********************************************************************************************************************************/
void
alloc2array(int *allocA,  const MatrixRect<int> *allocM,  const int *K,  const int &dim,  const int &nObs)
{
  int i, j;
  int *aAP;
  const int *aMP, *KP;

  if (!dim || !nObs){
    return;
  }

  aMP = allocM->aconst();
  aAP = allocA;
  for (i = 0; i < nObs; i++){
    KP = K;
    for (j = 0; j < dim; j++){
      *aAP = *aMP;        /*** Shift from scale 0,...,2*K to scale -K,...,K:  *aAP = *aMP - (*KP); ***/
      aAP++;
      aMP++;
      KP++;
    }
  }

  return;
}


/***** alloc2mixtureN *************************************************************************************************************/
/*                                                                                                                                */
/* PURPOSE: From allocations, compute numbers of observations in each mixture component                                           */
/*                                                                                                                                */
/* mixtureN[dim]:      array of matrices, each of dimension 1 x 2*K[i]+1                                                          */
/* alloc[dim x nObs]:  allocations (on scale -K,...,K)                                                                            */
/*                                                                                                                                */
/**********************************************************************************************************************************/
void
alloc2mixtureN(MatrixRect<int> *mixtureN,  const int *alloc,  const int *K,  const int &dim,  const int &nObs)
{
  int j;

  /*** Reset mixtureN ***/
  int *mixNP;
  MatrixRect<int> *mixtureNP = mixtureN;
  const int *KP = K;
  for (j = 0; j < dim; j++){
    mixNP = mixtureNP->a();
    for (int k = -(*KP); k <= *KP; k++){
      *mixNP = 0;
      mixNP++;
    }
    mixtureNP++;
    KP++;
  }

  /*** Compute mixtureN ***/ 
  const int *allocP = alloc;
  for (int i = 0; i < nObs; i++){
    mixtureNP = mixtureN;
    KP = K;
    for (j = 0; j < dim; j++){    
      mixtureNP->a()[*allocP + (*KP)]++;
      mixtureNP++;
      KP++;
      allocP++;
    }
  }

  return;
}


/***** updateAlloc1 ************************************************************************************************************/
/*                                                                                                                             */
/* PURPOSE:  Update allocations                                                                                                */
/*           This version assumes that given model parameters, margins are independent                                         */
/*           * useful when multivariate distributions are modelled using copulas                                               */
/*                                                                                                                             */
/*  Alloc[dim x nObs]:                matrix dim x nObs stored in an array with allocations to be updated                      */
/*  MixtureN[dim, 1 x length[i]]:     array of matrices 1 x length[i] with numbers of observations allocated in each component */
/*  Gspl:                             Gspline,                                                                                 */
/*                    NEEDED UNMODIFIED ITEMS: _d_knots, _inv_sigma2_d2                                                        */
/*                    MODIFIED ITEMS:          _workD, _workI (storage space)                                                  */
/*  Obs[dim x nObs]:                  matrix dim x nObs stored in an array (column major order) with observations              */
/*  EObs[dim]:                        vector of intercepts (means) of observations                                             */
/*  nObs:                             number of observations                                                                   */
/*                                                                                                                             */
/* =========================================================================================================================== */
void
updateAlloc1(int *Alloc,         MatrixRect<int> *MixtureN,  Gspline2 *Gspl,
             const double *Obs,  const double *EObs,         const int &nObs)
{
  static int i, j, k;

  /*** Reset MixtureN ***/
  const int *KP = Gspl->KAconst();
  MatrixRect<int> *mixtureNM = MixtureN;
  int *mixtureNA;
  for (j = 0; j < Gspl->dim(); j++){
    mixtureNA = mixtureNM->a();
    for (k = -(*KP); k <= *KP; k++){
      *mixtureNA = 0;
      mixtureNA++;
    }
    KP++;
    mixtureNM++;
  }

  /*** Loop over observations and margins. ***/
  int *allocP      = Alloc;
  const double *bP = Obs;

  const double *EObsP;
  const MatrixRect<double> *d_knotsM, *aM;
  const double *log_null_wP, *d_knotsP, *inv_sigma2_d2P, *aP;

  double *logw, *cumw, max_logw;
  int *Index, n_nonzero_comp, index;

  for (i = 0; i < nObs; i++){           /** loop over observations **/
    KP = Gspl->KAconst();
    d_knotsM = Gspl->d_knotsconst();
    aM = Gspl->aconst();
    inv_sigma2_d2P = Gspl->inv_sigma2_d2Aconst();
    log_null_wP = Gspl->log_null_wAconst();
    EObsP = EObs;
    mixtureNM = MixtureN;

    for (j = 0; j < Gspl->dim(); j++){        /** loop over margins **/

      /*** Compute log-weights. ***/
      d_knotsP = d_knotsM->aconst();
      aP       = aM->aconst();
      logw     = Gspl->workDA();
      max_logw = R_NegInf;
      for (k = -(*KP); k <= *KP; k++){
        *logw = *bP - (*EObsP) - (*d_knotsP);     
        *logw *= (*logw);
        *logw *= -0.5 * (*inv_sigma2_d2P);
        *logw += *aP;
        if (*logw > max_logw) max_logw = *logw;

        d_knotsP++;
        aP++;
        logw++;
      }

      /*** Rescale log-weights such that the highest one will be equal to zero. ***/
      /*** Keep only the weights that are higher than 'null' weight.            ***/
      /*** Compute cumulative sums of nonzero weights.                          ***/
      logw = Gspl->workDA();
      cumw = logw;
      Index = Gspl->workIA();
      n_nonzero_comp = 0;
      index = -(*KP);
      while (n_nonzero_comp == 0){
        *logw -= max_logw;
        if (*logw > *log_null_wP){
          *Index = index;
          *cumw = exp(*logw);

          n_nonzero_comp++;
          Index++;                    
          cumw++;
        }
        logw++;
        index++;
      }
      while (index <= *KP){
        *logw -= max_logw;
        if (*logw > *log_null_wP){
          *Index = index;
          *cumw = *(cumw-1) + exp(*logw);

          n_nonzero_comp++;
          Index++;                    
          cumw++;
        }
        logw++;
        index++;
      }
      if (n_nonzero_comp == 0){                               /* This is impossible when _log_null_w[j] < 0 */
        throw returnR("Trap in util_Gspline2.cpp:updateAlloc1. All weights are zero for some observation", 1);
      }
 
      /*** Sample new allocation ***/
      Random::discreteSampler2(&index, Gspl->workDA(), &n_nonzero_comp, &_AK_ONE_INT, &_AK_ONE_INT);
      *allocP = Gspl->workIA()[index];

      /*** Update mixtureN ***/
      mixtureNA = mixtureNM->a();
      mixtureNA[*allocP + (*KP)]++;
      
      /*** Increase pointers ***/
      KP++;
      d_knotsM++;
      aM++;
      inv_sigma2_d2P++;
      log_null_wP++;
      EObsP++;
      mixtureNM++;
 
      allocP++;
      bP++;
    }                   /*** end of the loop over margins ***/
  }                  /*** end of the loop over observations ***/

  return;
}

}
