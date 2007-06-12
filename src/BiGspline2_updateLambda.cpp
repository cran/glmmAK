/*** BiGspline2_updateLambda.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  27/03/2007
//                            
//                          updateLambda:  28/03/2007
//
// PURPOSE: Update of the smoothing hyperparameters
//
/* ********************************************************************************* */

#include "BiGspline2.h"

/***** updateLambda:  Update of the smoothing hyperparameter *****/
/***                                                           ***/
/*** ========================================================= ***/
void
BiGspline2::updateLambda()
{
  static int k;
  static double shape, rate, scale;

  if (_order == 0) return;                     // a's are fixed -> there is no random lambda

  switch (_prior_for_lambda){
  case GMRF_Gspline_Util::_Fixed_:
    return;

  case GMRF_Gspline_Util::_Gamma_:
    if (_equal_lambda){
      shape = _lambdaPriorShape[0] + 0.5*(_total_length - _order + 1);
      rate  = _lambdaPriorInvScale[0] - _penalty[0]; 
      if (rate <= 0) throw returnR("Trap in BiGspline2_updateLambda.cpp: BiGspline2::updateLambda(). Non-positive rate parameter", 1);
      scale = 1/rate;
      _lambda[0] = rgamma(shape, scale);
      for (k = 1; k < BiGspline2A::_dim; k++) _lambda[k] = _lambda[0];
    }
    else{
      //Rprintf("penalty=%g, %g;  ", _penalty[0], _penalty[1]);
      for (k = 0; k < BiGspline2A::_dim; k++){
        shape = _lambdaPriorShape[k] + 0.5*(0.5*_total_length - _order + 1);
        rate  = _lambdaPriorInvScale[k] - _penalty[k]; 
        if (rate <= 0) throw returnR("Trap in BiGspline2_updateLambda.cpp: BiGspline2::updateLambda(). Non-positive rate parameter", 1);
        scale = 1/rate;
        _lambda[k] = rgamma(shape, scale);
        //Rprintf("shape[%d]=%g, rate[%d]=%g, lambda[%d]=%g;  ", k, shape, k, rate, k, _lambda[k]);
      }
      //Rprintf("\n");
    }
    return;

  case GMRF_Gspline_Util::_SDUnif_:
    if (_equal_lambda){
      shape = 0.5*(_total_length - _order);
      rate  = - _penalty[0]; 
      if (rate <= 0) throw returnR("Trap in BiGspline2_updateLambda.cpp: BiGspline2::updateLambda(). Non-positive rate parameter", 1);
      scale = 1/rate;
      GMRF_Gspline_Util::rltruncGamma(_lambda, &shape, &scale, _lambdaPriorShape);
      for (k = 1; k < BiGspline2A::_dim; k++) _lambda[k] = _lambda[0];
    }
    else{
      for (k = 0; k < BiGspline2A::_dim; k++){
        shape = 0.5*(0.5*_total_length - _order);
        rate  = - _penalty[k]; 
        if (rate <= 0) throw returnR("Trap in BiGspline2_updateLambda.cpp: BiGspline2::updateLambda(). Non-positive rate parameter", 1);
        scale = 1/rate;
        GMRF_Gspline_Util::rltruncGamma(_lambda + k, &shape, &scale, _lambdaPriorShape + k);
      }
    }
    return;

  default:     
    throw returnR("Error in Gspline2_updateLambda.cpp: Gspline2::updateLambda(). Unimplemented prior", 1);
  }
}
