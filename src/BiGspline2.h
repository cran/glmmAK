/*** BiGspline2.h ***/

#ifndef _BI_G_SPLINE_TWEE_H_
#define _BI_G_SPLINE_TWEE_H_

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
#include "GMRF_Gspline_Util.h"

namespace BiGspline2A {
  enum neighborSystem {uniCAR, eight_neighbors, twelve_neighbors};            // types of neighboring systems for 'a' coefficients

  const int _dim   = 2;
  const int _LTdim = 3;        // = (_dim*(_dim+1))/2
}

class BiGspline2
{
  private:

  /** Dimensionality parameters **/
  int _K[BiGspline2A::_dim];          // (_length[j] - 1)/2 for each margin                                              _dim
  int _length[BiGspline2A::_dim];     // lengths of the G-spline in each margin (it should be odd)                       _dim
  int _max_length;                    // max(_length)                                                                    1
  int _total_length;                  // _length[0]*_length[1]                                                            1


  /** G-spline basis parameters **/
  MatrixRect<double> _knots[BiGspline2A::_dim];      // knots (means) for each margin                             _dim, 1 x _length[i]
  MatrixRect<double> _knots2[BiGspline2A::_dim];     // knots^2 (means^2) for each margin                         _dim, 1 x _length[i]
  double _invsigma2[BiGspline2A::_dim];              // basis sigma^{-2} for each margin                          _dim
  double _sigma[BiGspline2A::_dim];                  // basis sigma for each margin                               _dim


  /** G-spline parameters in connection with overall standard deviations           **/
  /**   NOTATION: overall covariance matrix = D = diag(d[i]^2)                     **/
  MatrixRect<double> _d_knots[BiGspline2A::_dim];        // d[i] * knots[i,...] for each margin                                _dim, 1 x _length[i]
  double _sigma_d[BiGspline2A::_dim];                    // sigma[i]*d[i]   (total standard deviation in each margin)          _dim
  double _inv_sigma2_d2[BiGspline2A::_dim];              // (sigma[i]*d[i])^{-2}   (total inverse variance in each margin)     _dim


  /** 'a' coefficients and related quantities **/
  MatrixRect<double> _a;                                  // a coefficients                  _length[0] x _length[1]
  MatrixRect<double> _expa;                               // exp(a) coefficients             _length[0] x _length[1]
  MatrixRect<double> _sumexpa_margin[BiGspline2A::_dim];  // sum of exp(a) over all indeces up to one index   _dim, 1 x _length[i]
  double _sumexpa;                                        // sum(exp(a))                     1
  double _a_max;                                          // max(a)                          1

  double _log_null_w;              // 'a' coefficient which satisfies (a - a_max < _log_null_w) will be assumed   1 x _dim
                                   //   to lead to zero weight
                                   // it is defined such that 
                                   //   (_total_length*exp(_log_null_w))/(1+_total_length*exp(_log_null_w)) = _null_mass


  int _k_effect;                      // number of weights that correspond to non-zero mass                            1
  MatrixRect<int> _ind_w_effect;      // binary indicators of components with non-zero mass      _length[0] x _length[1]


  /** Identifiability of 'a' coefficients **/
  int _aReference[BiGspline2A::_dim];    // indeces of reference knots in each margin (ignored during the MCMC)      _dim
                                         //  indeces are stored on scale -K[i],...,K[i]


  /** Quantities related to the update of a coefficients **/
  int _type_update_a;                //  value from the enum 'BiGspline2A::a_sampling_scheme'                         1

  MatrixRect<double> _abscis;        //  matrix with starting abscisae for ARS                         BiGspline2A::_nabscis x _total_length
                                                     //  it is used to store intervals defining the slice when the slice sampler is used
  int _iwv[mcmc_Gspline2::_liwv];                    // integer working place for ARS
  double _rwv[mcmc_Gspline2::_lrwv];                 // double working place for ARS
  double _hx[mcmc_Gspline2::_nabscis];               // working array for ARS/slice sampler/update_sigma
  double _hpx[mcmc_Gspline2::_nabscis];              // working array for ARS/slice sampler/update_sigma


  /** Penalty related quantities **/
  double _lambda[BiGspline2A::_dim];               // lambda precision parameters for each margin                          _dim
  double _penalty[BiGspline2A::_dim];              // -0.5*sum(Delta a)^2 for each margin if uniCAR                        _dim
                                                   //  * if _equal_lambda then _penalty[0] gives the sum of row- and column-penalties
                                                   //  * if _equal_lambda then _penalty[1], ... are sometimes used as working space
  int _neighbor_system;                            // neighboring system for the penalty
  int _order;                                      // order of autoregression in the penalty for each margin                               1
                                                   //    if equal to 0, a coefficients are assumed to be fixed
  bool _equal_lambda;                              // are lambda's for all margins (if _dim > 1) equal?                                    1
  int _prior_for_lambda;                           // type of prior for lambda(s)                                                          1
  double _lambdaPriorShape[BiGspline2A::_dim];     // if _prior_for_lambda = Gamma_: shape parameters of gamma priors                      _dim
                                                   // if _prior_for_lambda = SDUnif_: 1/S^2, where S is the upper limit 
                                                   //                                        of the uniform prior for lambda^{-1/2}
  double _lambdaPriorInvScale[BiGspline2A::_dim];  // if _prior_for_lambda = Gamma_: inverse scale (= rate) parameters of gamma priors     _dim
 

  /** Mean and variance of the G-spline (not taking into account general shift and general standard deviation),           **/
  /** separately for each margin                                                                                          **/
  /** _gMean[j]   = sum(w[j,k]*mu[j,k])                                                                                   **/
  /** _gVar[j]    = sum(w[j,k]*mu[j,k]^2) - _gMean[j]^2 + sigma[j]^2                                                      **/
  /**             = sum(w[j,k]*(mu[j,k] - _gMean[j])^2) + sigma[j]^2                                                      **/
  /**                                                                                                                     **/
  /** UPDATED ONLY WHEN WRITING TO FILES, NOT AFTER EACH WEIGHTS UPDATE                                                   **/
  double _gMean[BiGspline2A::_dim];                    // _dim
  double _gVar[BiGspline2A::_LTdim];                   // _LTdim (lower triangle of a symmetric matrix stored in column major order)


  /** Working spaces **/
  MatrixRect<double> _workD;                 //        1 x _total_length
  MatrixRect<int> _workI;                    //        1 x _total_length

  public:

  /***** Constructors and destructors *****/
  BiGspline2();
  BiGspline2(const int *dimPar,  
             const int *penMix_iPar,     const double *penMix_dPar,  
             const int *aPar,
             const double *lambdaPrior,  const double *lambda_a,
             const double *invVar);
  BiGspline2(const BiGspline2 &Gs);
  BiGspline2& operator=(const BiGspline2 &Gs);
  ~BiGspline2();


  /*****  Access functions to various components *****/
  inline const int*
  lengthconst() const { return _length;}

  inline int
  length(const int &i) const { return _length[i];}

  inline int
  total_length() const { return _total_length;}

  inline bool
  equal_lambda() const { return _equal_lambda;}

  inline double
  log_null_w() const { return _log_null_w;}

  inline const int*
  Kconst() const { return _K;}

  inline const double*
  aAconst() const { return _a.aconst();}

  inline const MatrixRect<double>*
  knotsconst() const { return _knots;}

  inline const MatrixRect<double>*
  d_knotsconst() const { return _d_knots;}

  inline MatrixRect<double>*
  d_knots() { return _d_knots;}

  inline const double*
  sigmaconst() const { return _sigma;}

  inline const double*
  invsigma2const() const { return _invsigma2;}

  inline double*
  invsigma2() { return _invsigma2;}

  inline double*
  sigma_d() { return _sigma_d;}

  inline const double*
  inv_sigma2_d2const() const { return _inv_sigma2_d2;}

  inline double*
  inv_sigma2_d2() { return _inv_sigma2_d2;}

  inline const double*
  lambdaconst() const { return _lambda;}

  inline const double*
  gMeanconst() const { return _gMean;}

  inline const double*
  gVarconst() const { return _gVar;}

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
  update_aBi(double *aa,     int *a_log_inf,      double *expaa,       double *sexpaam0,      double *sexpaam1,   double *Abscis,
             const int *ia,  const int *a_ipars,  const int *overrelax);

  void
  update_alla_lambdaBi(int *a_ipars,  const int *mixtureNM,  const int *overrelax);

  void
  full_a_parsBi(double* mean,  double* invvar,  const int *ia,  const double *aa) const;

  void
  find_start_abscisBi();


  /***** Utilities                                             *****/
  /***** ===================================================== *****/

  /*** BiGspline2_Util.cpp: ***/
  void
  gMoments();

  void
  penalty();

  void
  write_logweights2file(std::ofstream& afile,  std::ofstream& wfile,  std::ofstream& indfile,
                        const int *writea,     const int& prec,       const int& width) const;

  void
  print() const;
};

#endif
