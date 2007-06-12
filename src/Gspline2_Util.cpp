/*** Gspline2_Util.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  15/01/2007
//                            
//                                     gMoments:  ??/10/2006
//                                      penalty:  17/10/2006
//                        write_logweights2file:  ??/10/2006
//                                        print:  ??/10/2006
//
// PURPOSE: Utilities for class Gspline2
//
/* ********************************************************************************* */

#include "Gspline2.h"


/***** gMoments:  Compute mean and variance of the G-spline from the scratch                *****/
/***                                                                                          ***/
/***    * Overall intercept (alpha) and scale (tau) parameters are not taken into account     ***/
/***                                                                                          ***/
/*** ======================================================================================== ***/
void
Gspline2::gMoments()
{
  int k;
  double mu_EY;

  double *gMeanP = _gMean.a();
  double *gVarP  = _gVar.a();
  const int *KP                     = _K.aconst();  
  const double *sigmaP              = _sigma.aconst();
  const MatrixRect<double> *knotsM  = _knots;
  const MatrixRect<double> *expaM   = _expa;
  const double *sumexpaP            = _sumexpa.aconst();
  const double *knotsP, *expaP;

  for (int j = 0; j < _dim; j++){

    /*** Expectation[j] ***/
    knotsP  = knotsM->aconst();
    expaP   = expaM->aconst();
    *gMeanP = 0.0;
    for (k = -(*KP); k <= (*KP); k++){
      *gMeanP += (*expaP) * (*knotsP);
      knotsP++;
      expaP++;
    }
    *gMeanP /= *sumexpaP;

    /*** Variance [j] ***/
    knotsP  = knotsM->aconst();
    expaP   = expaM->aconst();
    *gVarP  = 0.0;
    for (k = -(*KP); k <= (*KP); k++){
      mu_EY = *knotsP - (*gMeanP);
      *gVarP += (*expaP) * mu_EY * mu_EY;
      knotsP++;
      expaP++;
    }
    *gVarP  /= *sumexpaP;
    *gVarP  += (*sigmaP) * (*sigmaP);

    gMeanP++;
    gVarP++;
    KP++;
    sigmaP++;
    knotsM++;
    expaM++;
    sumexpaP++;
  }

  return;
}


/***** penalty                                                               *****/
/***                                                                           ***/
/*** Partially taken from Gspline_update_lambda.cpp of bayesSurv package       ***/
/***                                                                           ***/
/* ============================================================================= */
void
Gspline2::penalty(double *value, const double *a, const int *order, const int *na)
{
  int i;
  double *DaP;
  const double *aP;

  aP  = a;
  DaP = _workD.a();
  for (i = 0; i < *na; i++){
    *DaP = *aP;
    aP++;
    DaP++;
  }

  GMRF::diff(_workD.a(), order, na);                /*** from GMRF.h ***/
  *value = 0.0;
  DaP = _workD.a();
  for (i = 0; i < *na - (*order); i++){
    *value += (*DaP) * (*DaP);
    DaP++;
  }
  *value *= -0.5;

  return;
}


/***** write_logweights2file: Write log-weights to the file                  *****/
/***                                                                           ***/
/*** ========================================================================= ***/
void
Gspline2::write_logweights2file(std::ofstream& afile,  const int& prec,  const int& width) const
{
  static int k;

  const MatrixRect<double> *aM = _a;
  const int *KP = _K.aconst();
  const double *aP;  
  for (int j = 0; j < _dim; j++){
    aP = aM->aconst();
    for (k = -(*KP); k <= *KP; k++){
      if (*aP < 1 && *aP > -1 && *aP != 0){
        afile << std::scientific << std::setw(width) << std::setprecision(prec) << *aP;
        afile << "   ";
      }
      else{
        afile << std::fixed << std::setw(width) << std::setprecision(prec) << *aP;
        afile << "   ";
      }      
      aP++;
    }
    aM++;
    KP++;
  }
  afile << std::endl;

  return;
}


/***** print                                *****/
/***                                          ***/
/*** ======================================== ***/
void
Gspline2::print() const
{
  int j;

  Rprintf("\nObject of class Gspline2:\n");
  Rprintf("===========================\n");
  Rprintf("dim=%d,  sum_n_knots=%d,  max_n_knots=%d\n", _dim, _sum_length, _max_length);

  for (j = 0; j < _dim; j++){
    Rprintf("\nMargin %d (K=%d, n_knots=%d, sigma=%g, invsigma2=%g):\n", j, _K.a(j), _length.a(j), _sigma.a(j), _invsigma2.a(j));
    Rprintf("          knots: ");
    _knots[j].print(0);
    Rprintf("    log-weights: ");
    _a[j].print(0);
    Rprintf("\n    g-Mean=%g,  g-Var=%g", _gMean.a(j), _gVar.a(j));   
    Rprintf("\n    sum(exp(a))=%g,  log(null_w)=%g", _sumexpa.a(j), _log_null_w.a(j));
    Rprintf("\n");
  }

  return;
}


