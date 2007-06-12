/*** BiGspline2_Util.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  21/03/2007
//                            
//                                     gMoments:  22/03/2007
//                                      penalty:  22/03/2007
//                        write_logweights2file:  29/03/2007
//                                        print:  21/03/2007
//
// PURPOSE: Utilities for class BiGspline2
//
/* ********************************************************************************* */

#include "BiGspline2.h"


/***** gMoments:  Compute mean and (co-)variance of the G-spline from the scratch           *****/
/***                                                                                          ***/
/***    * Overall intercept (alpha) and scale (tau) parameters are not taken into account     ***/
/***                                                                                          ***/
/*** ======================================================================================== ***/
void
BiGspline2::gMoments()
{
  int j, k, k0, k1;
  double mu_EY;

  double *gMeanP = _gMean;
  double *gVarP  = _gVar;
  const int *KP                             = _K;
  const double *sigmaP                      = _sigma;
  const MatrixRect<double> *knotsM          = _knots;
  const MatrixRect<double> *sumexpa_marginM = _sumexpa_margin;
  const double *knotsP, *knots0P, *knots1P, *expaP;

  for (j = 0; j < BiGspline2A::_dim; j++){

    /*** Expectation[j] ***/
    knotsP  = knotsM->aconst();
    expaP   = sumexpa_marginM->aconst();
    *gMeanP = 0.0; 
    for (k = -(*KP); k <= (*KP); k++){
      *gMeanP += (*expaP) * (*knotsP);
      knotsP++;
      expaP++;
    }
    *gMeanP /= _sumexpa;

    /*** Variance [j] ***/
    knotsP  = knotsM->aconst();
    expaP   = sumexpa_marginM->aconst();
    *gVarP  = 0.0; 
    for (k = -(*KP); k <= (*KP); k++){
      mu_EY = *knotsP - (*gMeanP);
      *gVarP += (*expaP) * mu_EY * mu_EY;
      knotsP++;
      expaP++;
    }
    *gVarP /= _sumexpa;
    *gVarP  += (*sigmaP) * (*sigmaP);

    gMeanP++;
    gVarP += (BiGspline2A::_dim - j);
    KP++;
    sigmaP++;
    knotsM++;
    sumexpa_marginM++;
  }

  /*** Covariance ***/
  gVarP   = _gVar + 1;
  knots1P = _knots[1].aconst();
  expaP   = _expa.aconst();
  *gVarP  = 0.0;
  for (k1 = -_K[1]; k1 <= _K[1]; k1++){  
    knots0P = _knots[0].aconst();
    for (k0 = -_K[0]; k0 <= _K[0]; k0++){  
      *gVarP += (*expaP) * (*knots0P - _gMean[0]) * (*knots1P - _gMean[1]);
      knots0P++;
      expaP++;
    }
    knots1P++;
  }
  *gVarP /= _sumexpa;

  //Rprintf("_gVar is (%g, %g, %g)\n", _gVar[0], _gVar[1], _gVar[2]);

  return;
}


/***** penalty                                                               *****/
/***                                                                           ***/
/*** Partially taken from Gspline_update_lambda.cpp of bayesSurv package       ***/
/***                                                                           ***/
/* ============================================================================= */
//
// If _equal_lambda: _penalty[0] = overall penalty, in the case of uniCAR prior, it is
//   the sum of the penalty over rows and over columns
//
// If !_equal_lambda and uniCAR prior: _penalty[0] = penalty over rows with fixed columns
//                                     _penalty[1] = penalty over columns with fixed rows
//
void
BiGspline2::penalty()                                       
{
  static int i, j, nr, nc;                  /*** nr and nc should be at least 5 ***/                             
  static double *DaP;
  static const double *a0, *a1, *a2, *a3, *a4, *a5, *a6, *a7, *a8;     /***  a5 a4 a7  ***/
  nr = _length[0];                                                     /***  a1 a0 a3  ***/   
  nc = _length[1];                                                     /***  a6 a2 a8  ***/

  switch (_neighbor_system){
  case BiGspline2A::uniCAR:

    /*** Penalty over rows with fixed columns ***/
    _penalty[0] = 0.0;
    a0 = _a.aconst();
    for (j = 0; j < nc; j++){
      DaP = _workD.a();
      for (i = 0; i < nr; i++){
        *DaP = *a0;
        a0++;
        DaP++;
      }

      GMRF::diff(_workD.a(), &_order, &nr);
      DaP = _workD.a();
      for (i = 0; i < nr - _order; i++){
        _penalty[0] += (*DaP) * (*DaP);
        DaP++;
      }
    }
    _penalty[0] *= (-0.5);

    /*** Penalty over columns with fixed rows ***/
    _penalty[1] = 0.0;
    a0 = _a.aconst();
    for (i = 0; i < nr; i++){
      DaP = _workD.a();
      a1 = a0;
      for (j = 0; j < nc; j++){
        *DaP = *a1;
        a1 += nr;
        DaP++;
      }
      a0++;

      GMRF::diff(_workD.a(), &_order, &nc);
      DaP = _workD.a();
      for (j = 0; j < nc - _order; j++){
        _penalty[1] += (*DaP) * (*DaP);
        DaP++;
      }
    }
    _penalty[1] *= (-0.5);

    if (_equal_lambda)  _penalty[0] += _penalty[1];
    return;

  case BiGspline2A::eight_neighbors:  /*** penalty is computed via sum of weighted squared pairwise differences ***/
    _penalty[0] = 0.0;

    // Neighbors of sites (i, 0)
    // ===========================
    /* neighbors of site (0, 0) */
    a0 = _a.aconst();
    a2 = a0 + 1;        
    a3 = a2 + (nr-1);
    a8 = a3 + 1;
    _penalty[0] += (*a0 - *a2)*(*a0 - *a2) + (*a0 - *a3)*(*a0 - *a3) - (*a0 - *a8)*(*a0 - *a8);
    a0++;
    a2++;
    a3++;
    a7 = a3 - 1;
    a8++;

    /* neighbors (not yet included) of sites (i, 0), i=1,...,nr-2 */
    for (i = 1; i <= nr-2; i++){
      _penalty[0] += 2*(*a0 - *a3)*(*a0 - *a3) + (*a0 - *a2)*(*a0 - *a2) - (*a0 - *a8)*(*a0 - *a8) - (*a0 - *a7)*(*a0 - *a7);
      a0++;
      a2++;
      a3++;
      a7++;
      a8++;
    }

    /* neighbors (not yet included) of site (nr-1, 0) */
    _penalty[0] += (*a0 - *a3)*(*a0 - *a3) - (*a0 - *a7)*(*a0 - *a7);
    a0++;
    a2++;
    a3++;
    a7++;
    a8++;


    // Neighbors (not yet included) of sites (i, j), j=1,...,nc-2
    // ==========================================================
    for (j = 1; j <= nc-2; j++){
      /* neighbors (not yet included) of site (0, j) */
      _penalty[0] += 2*(*a0 - *a2)*(*a0 - *a2) + (*a0 - *a3)*(*a0 - *a3) - (*a0 - *a8)*(*a0 - *a8);
      a0++;
      a2++;
      a3++;
      a7++;
      a8++;

      /* neighbors (not yet included) of sites (i, j), i=1,...,nr-2 */
      for (i = 1; i <= nr-2; i++){
        _penalty[0] += 2*((*a0 - *a3)*(*a0 - *a3) + (*a0 - *a2)*(*a0 - *a2)) - (*a0 - *a7)*(*a0 - *a7) - (*a0 - *a8)*(*a0 - *a8);
        a0++;
        a2++;
        a3++;
        a7++;
        a8++;
      }

      /* neighbors (not yet included) of site (nr-1, j) */
      _penalty[0] += (*a0 - *a3)*(*a0 - *a3) - (*a0 - *a7)*(*a0 - *a7);
      a0++;
      a2++;
      a3++;
      a7++;
      a8++;
    }


    // Neighbors (not yet included) of sites (i, nc-1)
    // ===============================================
    /* neighbors (not yet included) of sites (i, nc-1), i=0,...,nr-2 */
    for (i = 0; i <= nr-2; i++){
      _penalty[0] += (*a0 - *a2)*(*a0 - *a2);
      a0++;
      a2++;
    }

    _penalty[0] *= (-0.5);
    _penalty[1] = 0;
    return;

  case BiGspline2A::twelve_neighbors:
    REprintf("_neighbor_system=%d (twelve_neighbors)\n", _neighbor_system);  
    throw returnR("Error in BiGspline2_Util.cpp: BiGspline2::penalty(...), not (yet) implemented _neighbor_system", 1);        
    return;

  default:
    REprintf("_neighbor_system=%d\n", _neighbor_system);
    throw returnR("Error in BiGspline2_Util.cpp: BiGspline2::penalty(...), incorrect _neighbor_system argument", 1);        
  }
}


/***** write_logweights2file: Write weights and log-weights to the file                  *****/
/***                                                                                       ***/
/*** ===================================================================================== ***/
//
// afile:     File with log-weights (written all and only if writeAll)
// wfile:     File with weights (only weights corresponding to _ind_w_effect are written)
// indfile:   File with _k_effect, _ind_w_effect
// writea:    If <> 0 then log-weights are written to the file
//
void
BiGspline2::write_logweights2file(std::ofstream& afile,  std::ofstream& wfile,    std::ofstream& indfile,
                                  const int *writea,     const int& prec,         const int& width) const
{
  static int ind;
  static double gewicht;
  static const double *dP;
  static const int *iP;

  /*** log-weights -> afile ***/
  if (*writea){
    dP = _a.aconst();
    for (ind = 0; ind < _total_length; ind++){
      if (*dP < 1 && *dP > -1 && *dP != 0){
        afile << std::scientific << std::setw(width) << std::setprecision(prec) << *dP;
        afile << "  ";
      }
      else{
        afile << std::fixed << std::setw(width) << std::setprecision(prec) << *dP;
        afile << "  ";
      }      
      dP++;        
    }
    afile << std::endl;    
  }

  /*** weights                  -> wfile   ***/
  /*** _k_effect, _ind_w_effect -> indfile ***/
  indfile << _k_effect;
  indfile << "  ";

  dP = _expa.aconst();
  iP = _ind_w_effect.aconst();
  for (ind = 0; ind < _total_length; ind++){
    if (*iP){
      indfile << ind;                                         // index on the scale 0, ..., _total_length
      indfile << "  ";
       
      gewicht = *dP / _sumexpa;
      wfile << std::scientific << std::setw(width) << std::setprecision(prec) << gewicht;
      wfile << "  ";        
    }
    dP++;
    iP++;
  }
  indfile << std::endl;
  wfile << std::endl;    

  return;
}


/***** print                                *****/
/***                                          ***/
/*** ======================================== ***/
void
BiGspline2::print() const
{
  int j, iv;
  double gCorr = _gVar[1]/sqrt(_gVar[0]*_gVar[2]);

  Rprintf("\nObject of class BiGspline2:\n");
  Rprintf("=============================\n");
  Rprintf("dim=%d,  total_n_knots=%d,  max_n_knots (over margins)=%d\n", BiGspline2A::_dim, _total_length, _max_length);
  Rprintf("sum(exp(a))=%g,  max(a)=%g,  log(w[NULL])=%g\n", _sumexpa, _a_max, _log_null_w);
  Rprintf("equal_lambda=%s\n", _equal_lambda ? "true" : "false");

  for (j = 0; j < BiGspline2A::_dim; j++){
    Rprintf("\nMargin %d (K=%d, n_knots=%d, sigma=%g, invsigma2=%g):\n", j, _K[j], _length[j], _sigma[j], _invsigma2[j]);
    Rprintf("          knots: ");
    _knots[j].print(0);
    Rprintf("          knots^2: ");
    _knots2[j].print(0);
    iv = (j == 0 ? 0 : 2);
    Rprintf("\n    g-Mean=%g,  g-Var=%g", _gMean[j], _gVar[iv]);   
    Rprintf("\n    lambda=%g,  penalty=%g", _lambda[j], _penalty[j]);
    Rprintf("\n    shape(lambda)=%g,  rate(lambda)=%g", _lambdaPriorShape[j], _lambdaPriorInvScale[j]);

    Rprintf("\nSum(exp(a)) for this margin:");
    _sumexpa_margin[j].print(0);
    Rprintf("\n");
  }
  Rprintf("\n    g-Covariance=%g,  g-Correlation=%g\n", _gVar[1], gCorr);

  Rprintf("\n    k_effect=%d", _k_effect);
  //Rprintf("\n    Indices of non-zero components:\n");
  //_ind_w_effect.printI(0);

  //Rprintf("\nLog-weights:\n");
  //_a.print(0);
  //Rprintf("\nWeights (exp(a)):\n");
  //_expa.print(0);

  //Rprintf("\nStarting abscissae:\n");
  //_abscis.print(0);

  Rprintf("\n");

  return;
}

