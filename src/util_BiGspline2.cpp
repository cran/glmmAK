/*** util_BiGspline2.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  29/03/2007
//
//                      array2allocBi:   29/03/2007
//                      alloc2arrayBi:   29/03/2007
//                   alloc2mixtureNBi:   29/03/2007
//                      updateAllocBi:   29/03/2007
//
// PURPOSE: Bivariate G-spline, several smaller functions related to it
//
/* ********************************************************************************* */

#include "util_BiGspline2.h"

namespace util_BiGspline2{

/***** array2allocBi **************************************************************************************************************/
/*                                                                                                                                */
/* PURPOSE: Create a matrix of allocations                                                                                        */
/*          and perform input checks                                                                                              */
/*          Double index allocations ((-K0,_K1),...) -> single index allocations (0, ...)                                         */
/*                                                                                                                                */
/**********************************************************************************************************************************/
//
// allocM[1 x nObs]:    OUTPUT: matrix 1 x nObs with allocations on scale 0, ..., (2*K[0]+1)*(2*K[1]+1)-1
// allocA[dim * nObs]:  INPUT: array with allocations on scale -K,...,K, 
//                             = matrix BiGspline2A::_dim x nObs stored in column major order, 
//                               each column = allocations for each observation
// K[dim]:              numbers of knots on each side of the reference knot for each margin
// nObs:                number of observations
//
void
array2allocBi(MatrixRect<int> *allocM,  const int *allocA,  const int *K,  const int &nObs)
{
  int i, j;
  int *aMP;
  const int *aAP;

  int length0 = 2*K[0] + 1;

  if (BiGspline2A::_dim != 2){
    REprintf("BiGspline2A::_dim=%d\n");
    throw returnR("Error in util_BiGspline2.cpp: array2allocBi(). Implemented only for bivariate G-splines.", 1);
  }

  if (!nObs){
    *allocM = MatrixRect<int>(0, 0);
    return;
  }

  *allocM = MatrixRect<int>(1, nObs);
  aMP = allocM->a();
  aAP = allocA;
  for (i = 0; i < nObs; i++){
    if (aAP[0] < -K[0] || aAP[0] > K[0] || aAP[1] < -K[1] || aAP[1] > K[1]){
      REprintf("K=(%d, %d),  alloc[%d]=(%d, %d)\n", K[0], K[1], i, aAP[0], aAP[1]);
      throw returnR("Error in util_BiGspline2.cpp: array2allocBi. Allocation out of the range.", 1);
    }
    *aMP = (aAP[1] + K[1])*length0 + (aAP[0] + K[0]);
    aAP += BiGspline2A::_dim;
    aMP++;
  }

  return;
}


/***** alloc2arrayBi **************************************************************************************************************/
/*                                                                                                                                */
/* PURPOSE: Single index allocations (0, ...) -> Double index allocations ((-K0, -K1), ...)                                       */
/*                                                                                                                                */
/**********************************************************************************************************************************/
//
// allocA[dim * nObs]:  OUTPUT: array with allocations on scale -K,...,K, 
//                              = matrix BiGspline2A::_dim x nObs stored in column major order, 
//                                each column = allocations for each observation
// allocM[1 x nObs]:    INPUT: matrix 1 x nObs with allocations on scale 0, ..., (2*K[0]+1)*(2*K[1]+1)-1
//
void
alloc2arrayBi(int *allocA,  const MatrixRect<int> *allocM,  const int *K,  const int &nObs)
{
  int i, j;
  int *aAP;
  const int *aMP;

  int length0 = 2*K[0] + 1;

  if (BiGspline2A::_dim != 2){
    REprintf("BiGspline2A::_dim=%d\n");
    throw returnR("Error in util_BiGspline2.cpp: array2allocBi(). Implemented only for bivariate G-splines.", 1);
  }

  if (!nObs){
    return;
  }

  aMP = allocM->aconst();
  aAP = allocA;
  for (i = 0; i < nObs; i++){
    aAP[0] = *aMP % length0 - K[0];
    aAP[1] = *aMP / length0 - K[1];
    aAP += BiGspline2A::_dim;
    aMP++;    
  }

  return;
}


/***** alloc2mixtureNBi ***************************************************************************************************************/
/*                                                                                                                                    */
/* PURPOSE: From allocations, compute numbers of observations in each mixture component                                               */
/*                                                                                                                                    */
/**************************************************************************************************************************************/
//
// mixtureN[length0 x length1]:    matrix with numbers of observations per mixture component, length0 = 2*K[0]+1, length1 = 2*K[1]+1
// allocS[nObs]:                   allocations (single indeces 0, ...)
// total_length:                   = length0 * length1
//
void
alloc2mixtureNBi(int *mixtureN,  const int *allocS,  const int &total_length,  const int &nObs)
{
  static int j;
  static int *mixNP;
  static const int *allocP;

  /*** Reset mixtureN ***/
  mixNP = mixtureN;
  for (j = 0; j < total_length; j++){
    mixNP = 0;
    mixNP++;
  }

  /*** Compute mixtureN ***/ 
  allocP = allocS;
  mixNP  = mixtureN;
  for (j = 0; j < nObs; j++){
    mixNP[*allocP]++;
    allocP++;
  }

  return;
}


/***** updateAllocBi ***********************************************************************************************************/
/*                                                                                                                             */
/* PURPOSE:  Update allocations in the model with BIVARIATE G-spline                                                           */
/*                                                                                                                             */
/* =========================================================================================================================== */
//
//  ASSUMPTION: dim = 2
//  -------------------
//
//  allocS[nObs]:                 Single index (0,...) allocations to be updated
//  MixtureN[length0*length1]:    Array from the matrix length0 x length1 with numbers of observations allocated in each component
//                                (length0=2*K[0]+1, length1=2*K[1]+1, length0*length1 = Gspl->total_length())
//  Gspl:                         Gspline,
//                    NEEDED UNMODIFIED ITEMS: _d_knots, _inv_sigma2_d2
//                    MODIFIED ITEMS:          _workD, _workI (storage space)
//  Obs[dim*nObs]:                Array from the matrix dim x nObs (column major order) with observations
//  EObs[dim]:                    Vector of intercepts (alpha parameters) of observations
//  nObs:                         Number of observations
//
void
updateAllocBi(int *allocS,        int *mixtureN,       BiGspline2 *Gspl,
              const double *Obs,  const double *EObs,  const int &nObs)
{
  static int i, k0, k1, K0, K1, total_length;
  static double b_Eb0, b_Eb1, inv_sigma2_d20, inv_sigma2_d21, z0, z1, log_null_w;  
  static int *mixNP;
  static int *allocP;
  static const MatrixRect<double>* d_knotsM;
  static const double *d_knots0, *d_knots1, *bP, *EbP, *aP;  
  static double *logw, *cumw, max_logw;
  static int *Index, n_nonzero_comp, index;

  K0 = Gspl->Kconst()[0];
  K1 = Gspl->Kconst()[1];
  d_knotsM = Gspl->d_knotsconst();
  inv_sigma2_d20 = Gspl->inv_sigma2_d2const()[0];
  inv_sigma2_d21 = Gspl->inv_sigma2_d2const()[1];
  log_null_w     = Gspl->log_null_w();
  total_length   = Gspl->total_length();

  /*** Reset mixtureN ***/
  mixNP = mixtureN;
  for (k0 = 0; k0 < total_length; k0++){
    *mixNP = 0;
    mixNP++;
  }

  /*** Loop over observations ***/
  allocP = allocS;
  bP     = Obs;
  for (i = 0; i < nObs; i++){   

    /*** Compute log-weights ***/ 
    b_Eb0 = bP[0] - EObs[0];
    b_Eb1 = bP[1] - EObs[1];

    aP       = Gspl->aAconst();
    logw     = Gspl->workDA();
    max_logw = R_NegInf;
    d_knots1 = d_knotsM[1].aconst();

    for (k1 = -K1; k1 <= K1; k1++){
      d_knots0 = d_knotsM[0].aconst();
      for (k0 = -K0; k0 <= K0; k0++){
        z0 = b_Eb0 - (*d_knots0);
        z0 *= z0;
        z0 *= inv_sigma2_d20;

        z1 = b_Eb1 - (*d_knots1);
        z1 *= z1;
        z1 *= inv_sigma2_d21;
        
        *logw = -0.5*(z0 + z1) + (*aP);        
        if (*logw > max_logw) max_logw = *logw;

        d_knots0++;
        aP++;
        logw++;
      }
      d_knots1++;
    }

    /*** Rescale log-weights such that the highest one will be equal to zero. ***/
    /*** Keep only the weights that are higher than 'null' weight.            ***/
    /*** Compute cumulative sums of nonzero weights.                          ***/
    logw = Gspl->workDA();
    cumw = logw;
    Index = Gspl->workIA();
    n_nonzero_comp = 0;
    index = 0;
    while (n_nonzero_comp == 0){
      *logw -= max_logw;
      if (*logw > log_null_w){
        *Index = index;
        *cumw = exp(*logw);

        n_nonzero_comp++;
        Index++;                    
        cumw++;
      }
      logw++;
      index++;
    }
    while (index < total_length){
      *logw -= max_logw;
      if (*logw > log_null_w){
        *Index = index;
        *cumw = *(cumw-1) + exp(*logw);

        n_nonzero_comp++;
        Index++;                    
        cumw++;
      }
      logw++;
      index++;
    }
    if (n_nonzero_comp == 0){                               /* This is impossible when _log_null_w < 0 */
      throw returnR("Trap in util_BiGspline2.cpp:updateAllocBi. All weights are zero for some observation", 1);
    }

    /*** Sample new allocation ***/
    Random::discreteSampler2(&index, Gspl->workDA(), &n_nonzero_comp, &_AK_ONE_INT, &_AK_ONE_INT);
    *allocP = Gspl->workIA()[index];

    /*** Update mixtureN ***/
    mixtureN[*allocP]++;

    /*** Increase pointers ***/
    allocP++;
    bP += BiGspline2A::_dim;
  }

  return;
}

}   /*** end of the namespace util_BiGspline2 ***/

