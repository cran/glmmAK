/*** Gspline2_updateLambda.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//                 CREATED:  15/01/2007
//                            
//                          updateLambda:  15/01/2007
//
// PURPOSE: Update of the smoothing hyperparameters
//
/* ********************************************************************************* */

#include "Gspline2.h"


/***** updateLambda:  Update of the smoothing hyperparameter *****/
/***                                                           ***/
/*** ========================================================= ***/
void
Gspline2::updateLambda()
{
  static double shape, rate, scale;

  if (_dim == 1){
    if (_order.aconst()[0] == 0) return;     // a's are fixed -> there is no random lambda

    switch (_prior_for_lambda){
    case GMRF_Gspline_Util::_Fixed_:
      return;

    case GMRF_Gspline_Util::_Gamma_:
      shape = _lambdaPriorShape.aconst()[0] + 0.5*(_sum_length - _order.aconst()[0] + 1);
      rate  = _lambdaPriorInvScale.aconst()[0] - _penalty.aconst()[0]; 
      if (rate <= 0) throw returnR("Trap in Gspline2_updateLambda.cpp: Gspline2::updateLambda(). Non-positive rate parameter", 1);
      scale = 1/rate;
      _lambda.a()[0] = rgamma(shape, scale);
      return;

    case GMRF_Gspline_Util::_SDUnif_:
      shape = 0.5*(_sum_length - _order.aconst()[0]);
      rate  = - _penalty.aconst()[0]; 
      if (rate <= 0) throw returnR("Trap in Gspline2_updateLambda.cpp: Gspline2::updateLambda(). Non-positive rate parameter", 1);
      scale = 1/rate;
      GMRF_Gspline_Util::rltruncGamma(_lambda.a(), &shape, &scale, _lambdaPriorShape.aconst());
      return;

    default:     
      throw returnR("Error in Gspline2_updateLambda.cpp: Gspline2::updateLambda(). Unimplemented prior", 1);
    }
  }
  else{
    throw returnR("Error in Gspline2_updateLambda.cpp: Gspline2::updateLambda(). Not (yet) implemented for _dim <> 1", 1);
  }
}
