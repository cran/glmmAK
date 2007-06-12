/*** util_poisson.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  16/02/2007
//
//   PURPOSE:  Several smaller functions related to the Poisson log-linear regression model
//
//
/* ********************************************************************************* */

#include "util_poisson.h"

namespace util_poisson{


/******************* ecount_etas_poisson ***********************************/
/***                                                                     ***/
/***                                                                     ***/
/***************************************************************************/
void
ecount_etas_poisson(double *ecount,   double *eta,   const double *offset,  const double *theta,  const double *x,
                     const int &nObs,  const int &p)
{
  static int i, j;
  const double *thetaP;

  const double *offsetP = offset;
  const double *xP      = x;
  double *ecountP = ecount;
  double *etaP    = eta;

  for (i = 0; i < nObs; i++){
    thetaP = theta;
    *etaP  = 0.0;

    for (j = 0; j < p; j++){
      *etaP += (*xP) * (*thetaP);
      xP++;
      thetaP++;
    }

    *ecountP  = exp_AK(*etaP + (*offsetP));
    offsetP++;
    ecountP++;
    etaP++;      
  }

  return;  
}


/******************* ecount_etas_etaREs_poisson ***********************************/
/***                                                                            ***/
/***                                                                            ***/
/**********************************************************************************/
void
ecount_etas_etaREs_poisson(double *ecount,        double *etaF,          double *etaRE,          double *eta, 
                           const double *offset,  const double *theta,   const double *thetaRE,  
                           const double *x,       const double *xRE,     const int &nCluster,    const int *ni,
                           const int &p,          const int &pRE)
{
  static int i, j, k;
  const double *thetaP, *thetaREP;

  const double *thetaRE_now = thetaRE;
  const double *offsetP     = offset;
  const double *xP          = x;
  const double *xREP        = xRE;
  const int *niP            = ni;
  double *ecountP = ecount;
  double *etaFP   = etaF;
  double *etaREP  = etaRE;
  double *etaP    = eta;

  for (i = 0; i < nCluster; i++){

    for (k = 0; k < *niP; k++){

      /*** etaF ***/
      thetaP = theta;
      *etaFP  = 0.0;
      for (j = 0; j < p; j++){
        *etaFP += (*xP) * (*thetaP);
        xP++;
        thetaP++;
      }

      /*** etaRE ***/
      thetaREP = thetaRE_now;
      *etaREP  = 0.0;
      for (j = 0; j < pRE; j++){
        *etaREP += (*xREP) * (*thetaREP);
        xREP++;
        thetaREP++;
      }

      /*** eta and ecount ***/      
      *etaP     = *etaFP + (*etaREP);
      *ecountP  = exp_AK(*etaP + (*offsetP));
      offsetP++;
      ecountP++;
      etaP++;
      etaFP++;
      etaREP++;
    }

    thetaRE_now = thetaREP;
    niP++;
  }

  return;
}

}  /*** end of the namespace util_poisson ***/

