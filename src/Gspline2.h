/*** Gspline2.h ***/

#ifndef _G_SPLINE_TWEE_H_
#define _G_SPLINE_TWEE_H_

#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_Optim.h"

#include "MatrixLT.h"
#include "MatrixRect.h"

#include "mcmc_Gspline2.h"

#include "GMRF.h"
#include "GMRF_Gspline.h"
#include "GMRF_Gspline_Util.h"

class Gspline2
{
  private:
  /** Dimensionality parameters **/
  int  _dim;                  // dimension of the G-spline                                                        1 
  MatrixRect<int> _K;         // (_length[j] - 1)/2 for each margin                                               1 x _dim
  MatrixRect<int> _length;    // lengths of the G-spline in each margin (it should be odd)                        1 x _dim
  int _sum_length;            // sum(_length)                                                                     1
  int _max_length;            // max(_length)                                                                     1


  /** G-spline basis parameters **/
  MatrixRect<double> *_knots;      // knots (means) for each margin                                               _dim, 1 x _length[i]
  MatrixRect<double> *_knots2;     // knots^2 (means^2) for each margin                                           _dim, 1 x _length[i]
  MatrixRect<double> _invsigma2;   // basis sigma^{-2} for each margin                                            1 x _dim
  MatrixRect<double> _sigma;       // basis sigma for each margin                                                 1 x _dim


  /** G-spline parameters in connection with overall standard deviations           **/
  /**   NOTATION: overall covariance matrix = D = diag(d[i]^2)                     **/
  MatrixRect<double> *_d_knots;      // d[i] * knots[i,...] for each margin                                         _dim, 1 x _length[i]
  MatrixRect<double> _sigma_d;       // sigma[i]*d[i]   (total standard deviation in each margin)                   _dim
  MatrixRect<double> _inv_sigma2_d2; // (sigma[i]*d[i])^{-2}   (total inverse variance in each margin)              _dim


  /** 'a' coefficients and related quantities **/
  MatrixRect<double> *_a;          // a coefficients for each margin                                              _dim, 1 x _length[i]
  MatrixRect<double> *_expa;       // exp(a) coefficients for each margin                                         _dim, 1 x _length[i]
  MatrixRect<double> _sumexpa;     // sum(exp(a)) for each margin                                                 1 x _dim

  MatrixRect<double> _log_null_w;  // 'a' coefficient which satisfies (a - a_max < _log_null_w) will be assumed   1 x _dim
                                   //   to lead to zero weight
                                   // it is defined such that 
                                   //   (_length[i]*exp(_log_null_w[i]))/(1+_length[i]*exp(_log_null_w[i])) = _null_mass


  /** Identifiability of 'a' coefficients **/
  int _aIdent;                       // type of identifiability constraint                                             1
  MatrixRect<int> _aReference;       // indeces of reference knots in each margin (used when _aIdent=_Reference_)      1 x _dim
                                     //  indeces are stored on scale -K[i],...,K[i]
  MatrixRect<double> *_aContrast;    // identifiability contrasts for each margin                                      _dim, 1 x _length[i]
                                     //  these are vectors whose components sum up to 1


  /** Quantities related to the update of a coefficients **/
  int _type_update_a;                //  value from the enum 'Gspline2::a_sampling_scheme'                            1

  MatrixRect<double> *_abscis;       //  matrix with starting abscisae for ARS                                        _dim, Gspline2A::_nabscis x _length[i]
                                     //  it is used to store intervals defining the slice when the slice sampler is used
  MatrixRect<int> _iwv;              // integer working place for ARS                                                 1 x (7 + Gspline2A::_ns)
  MatrixRect<double> _rwv;           // double working place for ARS                                                  1 x (9 + 6*(1+Gspline2A::_ns))
  MatrixRect<double> _hx;            // working array for ARS/slice sampler/update_sigma                              1 x Gspline2A::_nabscis
  MatrixRect<double> _hpx;           // working array for ARS/slice sampler/update_sigma                              1 x Gspline2A::_nabscis


  /** Quantities needed for joint update of a coefficients **/
  /** !!! Currently only implemented for _dim = 1 !!!      **/
  /** !!! Set to NULL for _dim > 1 !!!                     **/
  MatrixRect<double> _w;            // weights (in a long array in column major order)                          1 x _sum_length
  double _minw;                     // minimal weight                                                           1
  MatrixLT<double>   _Q;            // Q = t(D)*D                                                                LT(_sum_length)
  MatrixRect<double> _Da;           // D*a in first _total_length - _order places                                1 x _sum_length
  MatrixRect<double> _Qa;           // Q*a = t(D)*D*a                                                            1 x _sum_length
  MatrixRect<int>    _diffOper;     // vector defining the difference operator                                   1 x (_order + 1) (if _dim = 1)

  double _par_rscale[6];            // parameters for the lambda proposal according to Rue and Held 
                                    // actually not used in this version

  int _LTna;                        // = LT(_sum_length)
  int _LTna_1;                      // = LT(_sum_length - 1)
  int _nworkML;
  int _nworka;
  int _nworkGMRF;
  MatrixRect<double> _workML;                  // working arrays, see GMRF_Gspline.cpp for how long they should be
  MatrixRect<double> _worka;
  MatrixRect<double> _workGMRF;


  /** Penalty related quantities **/
  MatrixRect<double> _lambda;               // lambda precision parameters for each margin                                          1 x _dim
  MatrixRect<double> _penalty;              // -0.5*sum(Delta a)^2 for each margin                                                  1 x _dim

  MatrixRect<int> _order;                   // order of autoregression in the penalty for each margin                               1 x _dim
                                            //    if equal to 0, a coefficients are assumed to be fixed
  bool _equal_lambda;                       // are lambda's for all margins (if _dim > 1) equal?                                    1
  int _prior_for_lambda;                    // type of prior for lambda(s)                                                          1
  MatrixRect<double> _lambdaPriorShape;     // if _prior_for_lambda = Gamma_: shape parameters of gamma priors                      _dim
                                            // if _prior_for_lambda = SDUnif_: 1/S^2, where S is the upper limit 
                                            //                                        of the uniform prior for lambda^{-1/2}
  MatrixRect<double> _lambdaPriorInvScale;  // if _prior_for_lambda = Gamma_: inverse scale (= rate) parameters of gamma priors     1 x _dim
  MatrixRect<double> _lambdaPriorPar;       // prior shape and inverse scale or 1/S^2 in the format as required by the argument
                                            //   par_lambda of function GMRF_Gspline::update(), stored in columns                   2 x _dim


  /** Mean and variance of the G-spline (not taking into account general shift and general standard deviation),           **/
  /** separately for each margin                                                                                          **/
  /** _gMean[j]   = sum(w[j,k]*mu[j,k])                                                                                   **/
  /** _gVar[j]    = sum(w[j,k]*mu[j,k]^2) - _gMean[j]^2 + sigma[j]^2                                                      **/
  /**             = sum(w[j,k]*(mu[j,k] - _gMean[j])^2) + sigma[j]^2                                                      **/
  /**                                                                                                                     **/
  /** UPDATED ONLY WHEN WRITING TO FILES, NOT AFTER EACH WEIGHTS UPDATE                                                   **/
  MatrixRect<double> _gMean;                  // 1 x _dim
  MatrixRect<double> _gVar;                   // 1 x _dim


  /** Working spaces **/
  MatrixRect<double> _workD;               //        1 x _max_length
  MatrixRect<int> _workI;                  //        1 x _max_length


  public:

  /***** Constructors and destructors *****/
  Gspline2();
  Gspline2(const int *dimPar,  
           const int *penMix_iPar,     const double *penMix_dPar,  
           const int *aPar,            const double *aContrast,
           const double *lambdaPrior,  const double *lambda_a,
           const double *invVar);
  Gspline2(const Gspline2 &Gs);
  Gspline2& operator=(const Gspline2 &Gs);
  ~Gspline2();

  /*****  Access functions to various components *****/
  inline int
  dim() const { return _dim;}

  inline int
  length(const int &i) const { return _length.a(i);}

  inline const int*
  KAconst() const { return _K.aconst();}

  inline const double*
  log_null_wAconst() const { return _log_null_w.aconst();}

  inline const MatrixRect<double>*
  aconst() const { return _a;}

  inline const MatrixRect<double>*
  knotsconst() const { return _knots;}

  inline const MatrixRect<double>*
  d_knotsconst() const { return _d_knots;}

  inline MatrixRect<double>*
  d_knots() { return _d_knots;}

  inline const double*
  sigmaAconst() const { return _sigma.aconst();}

  inline const double*
  invsigma2Aconst() const { return _invsigma2.aconst();}

  inline const double*
  sigma_dAconst() const { return _sigma_d.aconst();}

  inline double*
  sigma_dA() { return _sigma_d.a();}

  inline const double*
  inv_sigma2_d2Aconst() const { return _inv_sigma2_d2.aconst();}

  inline double*
  inv_sigma2_d2A() { return _inv_sigma2_d2.a();}

  inline const double*
  lambdaAconst() const { return _lambda.aconst();}

  inline const double*
  gMeanAconst() const { return _gMean.aconst();}

  inline const double*
  gVarAconst() const { return _gVar.aconst();}

  inline double*
  workDA() { return _workD.a();}

  inline int*
  workIA() { return _workI.a();}

  /***** MCMC update of parameters                               *****/
  /***** ======================================================= *****/

  /*** Gspline2_updateLambda.cpp: ***/
  void
  updateLambda();

  /*** Gspline2_updateWeights.cpp: ***/
  void
  update_a1(double *aa,     double *expaa,       double *Abscis,
            const int *ia,  const int *a_ipars,  const int *overrelax);

  void
  update_alla_lambda1(int *a_ipars,  const int *mixtureNM,  const int *overrelax);

  void
  full_a_pars1(double* mean,  double* invvar,  const int *ia,  const double *aa) const;

  void
  find_start_abscis1();


  /***** Utilities                                             *****/
  /***** ===================================================== *****/

  /*** Gspline2_Util.cpp: ***/
  void
  gMoments();

  void
  penalty(double *value, const double *a, const int *order, const int *na);

  void
  write_logweights2file(std::ofstream& afile,  const int& prec,  const int& width) const;

  void
  print() const;
};

#endif
