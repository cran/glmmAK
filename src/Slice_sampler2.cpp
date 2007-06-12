/*** Slice_sampler2.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  30/10/2006
//              as a copy of slice_sampler.cpp from the bayesSurv package
//
// PURPOSE: Univariate slice sampling and some additional supporting functions
//
/* ********************************************************************************* */
//
// History of the original slice_sampler.cpp:
//
// 27/11/2004: ss_stepping_out
//             ss_doubling
// 28/11/2004: ss_accept_doubling            
//             ss_shrinkage_sample
// 30/11/2004: ss_bisection_overrelax
// 01/12/2004: ss_exact_sample
//             ss_exact_overrelax
//
// Reference: Neal, R.M. (2003). Slice sampling. The Annals of Statistics 31, 705-767.
//
// REMARK: routine void (*eval)(const double* double*, const double*, const int*)
//                  is supposed to return g(x) = log(f(x)) in its second argument
//                  and further it is supposed to return -FLT_MAX (or any other number representing -Inf)
//                  in the case when x is out of the range of the distribution (e.g. when there are bounds)
//

#include "Slice_sampler2.h"

#ifdef __cplusplus
extern "C" {
#endif


namespace Slice_sampler{

// ***** ss_stepping_out: Function for finding an interval around the current point 
//                          using the "stepping out" procedure
//                         (Figure 3, Neal (2003))
// INPUT:
//   x0[1] ........... current point
//   z[1] ............ logarithm of the vertical level defining the slice
//   w[1] ............ estimate of the typical size of slice
//   m[1] ............ integer limiting the size of a slice to m*w
//   bound[2] ........ bounds for the distribution of x (if any)
//   is_bound[2] ..... indicators whether there are bounds on left and right
//   (*eval) ......... routine to compute g(x) = log(f(x))
//   deval[] ......... double paramaters for (*eval)
//   ieval[] ......... integer parameters for (*eval)
// 
// OUTPUT:
//   interv[2] ....... found interval
//   g_interv[2] ..... log(f(x)) evaluated in the limits of the interval (if there are no bounds)
//
void
ss_stepping_out(double* interv,       double* g_interv,
                const double* x0,     const double* z,     const double* w,  const int* m,
                const double* bound,  const int* is_bound, 
                void (*eval)(const double*, double*, const double*, const int*),
                const double* deval,  const int* ieval)
{
  static double u;
  static int n_left, n_right;
  
  /* Initial guess for the interval */
  u = unif_rand();
  interv[0] = (*x0) - (*w)*u;
  interv[1] = interv[0] + (*w);

  /* Get numbers of steps tried to left and to right */
  u = unif_rand();
  n_left = int(floor((*m)*u));
  n_right = (*m - 1) - n_left;

  /* Initial evaluation of g in the left and right limits of the interval */
  eval(interv,   g_interv,   deval, ieval);
  eval(interv+1, g_interv+1, deval, ieval);

  /* Step to left until leaving the slice */
  while (n_left > 0 && g_interv[0] > *z){
    interv[0] -= (*w);
    eval(interv, g_interv, deval, ieval);
    n_left--;
  }

  /* Step to right until leaving the slice */
  while (n_right > 0 && g_interv[1] > *z){
    interv[1] += (*w);
    eval(interv+1, g_interv+1, deval, ieval);
    n_right--;
  }

  if (is_bound[0] && interv[0] <= bound[0]) interv[0] = bound[0];     /* g_interv[0] is already equal to -FLT_MAX  */
  if (is_bound[1] && interv[1] >= bound[1]) interv[1] = bound[1];     /* g_interv[1] is already equal to -FLT_MAX  */

  return;
 }    /* end of the function stepping_out */


// ***** ss_doubling: Function for finding an interval around the current point
//                      using the "doubling" procedure
//                      (Figure 4, Neal (2003))
//
// INPUT:
//   x0[1] ........... current point
//   z[1] ............ logarithm of the vertical level defining the slice
//   w[1] ............ estimate of the typical size of slice
//   p[1] ............ integer limiting the size of a slice to 2^p*w
//   bound[2] ........ bounds for the distribution of x (if any)
//   is_bound[2] ..... is_bound[0:1] = indicators whether there are bounds on left and right
//   unimodal[1] ..... indicator whether the distribution at question is unimodal or not
//   (*eval) ......... routine to compute g(x) = log(f(x))
//   deval[] ......... double paramaters for (*eval)
//   ieval[] ......... integer parameters for (*eval)
//
// OUTPUT:
//   interv[2] ....... found interval
//   g_interv[2] ..... log(f(x)) evaluated in the limits of the interval (if there are no bounds)
//
void
ss_doubling(double* interv,       double* g_interv,
            const double* x0,     const double* z,      const double* w,     const int* p,
            const double* bound,  const int* is_bound,  const int* unimodal,
            void (*eval)(const double*, double*, const double*, const int*),
            const double* deval,  const int* ieval)
{
  static double u;
  static int n_step;
  static bool go_left, go_right, now_left;

  /* Initial guess for the interval */
  u = unif_rand();
  interv[0] = (*x0) - (*w)*u;
  interv[1] = interv[0] + (*w);

  n_step = *p;
  go_left = go_right = true;

  /* Initial evaluation of g in the left and right limits of the interval */
  eval(interv,   g_interv,   deval, ieval);
  eval(interv+1, g_interv+1, deval, ieval);
  if (is_bound[0] && interv[0] <= bound[0]) go_left = false;
  if (is_bound[1] && interv[1] >= bound[1]) go_right = false;
                    /* leave the value of interv[0] or interv[1] outside the range to be able to perform back-tracking */

  if (*unimodal){    
    if (g_interv[0] <= *z) go_left = false;            /* left limit is already outside the slice  */
    if (g_interv[1] <= *z) go_right = false;           /* right limit is already outside the slice */
  }
  if ((!go_left) && (!go_right)) n_step = 0;
  
  /* Perform doubling until both ends are outside the slice */
  while (n_step > 0 && (g_interv[0] > *z || g_interv[1] > *z)){
    if (go_right && go_left){    /* we have to decide where to go in this step */
      u = unif_rand();
      now_left = (u < 0.5);
    }
    else{
      if (!go_right) now_left = true;
      else           now_left = false;
    }

    if (now_left){
      interv[0] -= (interv[1] - interv[0]);
      eval(interv, g_interv, deval, ieval);
      if (is_bound[0] && interv[0] <= bound[0]) go_left = false;
      if ((*unimodal) && g_interv[0] <= *z)     go_left = false;  /* if unimodal distribution */
                                                                  /* left limit is already outside the slice, no slice points to the left */
    }    /* now_left = false */
    else{
      interv[1] += (interv[1] - interv[0]);
      eval(interv+1, g_interv+1, deval, ieval);
      if (is_bound[1] && interv[1] >= bound[1]) go_right = false;
      if ((*unimodal) && g_interv[1] <= *z)     go_right = false; /* if unimodal distribution */
                                                                /* right limit is already outside the slice, no slice points to the right */
    }

    n_step--;
    if ((!go_left) && (!go_right)) n_step = 0;
  }

  return;
}    /* end of the function ss_doubling */


// ***** ss_accept_doubling: Acceptance test of newly sampled point 
//                           when the "doubling" procedure has been used to find an interval to sample from
//                           (Figure 6, Neal (2003))
//
// REMARK: This procedure can be omitted for unimodal distributions
//                          
// INPUT:
//   x0[1] ........... the current point
//   x1[1] ........... the possible next point (candidate)
//   z[1] ............ logarithm of the vertical level defining the slice
//   w[1] ............ estimate of the typical size of slice
//   interv[2] ....... interval around the current point x0
//   z[1] ............ logarithm of the vertical level defining the slice
//   w[1] ............ estimate of the typical size of slice
//   p[1] ............ integer limiting the size of a slice to 2^p*w
//   (*eval) ......... routine to compute g(x) = log(f(x))
//   deval[] ......... double paramaters for (*eval)
//   ieval[] ......... integer parameters for (*eval)
//
// OUTPUT:
//   accept .......... 1/0 indicating whether the point is acceptable or not
//
void
ss_accept_doubling(int* accept,           const double* x0,        const double* x1,     
                   const double* z,       const double* w,         const double* interv,
                   void (*eval)(const double*, double*, const double*, const int*),
                   const double* deval,   const int* ieval)
{
  static double interv1[2];
  static double g_interv1[2];
  static bool diff__;
  static double w11, mid;

  w11 = 1.1*(*w);
  interv1[0] = interv[0];
  interv1[1] = interv[1];
  diff__ = false;

  *accept = 1;
  while ((*accept) && (interv1[1] - interv1[0] > w11)){
    mid = 0.5*(interv1[0] + interv[1]);
    if ((*x1) < mid){
      if ((*x0) >= mid) diff__ = true;
      interv1[1] = mid;
      eval(interv1+1, g_interv1+1, deval, ieval);
    }
    else{
      if ((*x0) < mid) diff__ = true;
      interv1[0] = mid;
      eval(interv1, g_interv1, deval, ieval);
    }
    if (diff__ &&  g_interv1[0] <= (*z) && g_interv1[1] <= (*z)) *accept = 0;
  }

  return;
}   /* end of the function ss_accept_doubling */


// ***** ss_shrinkage_sample: Function to sample a point from the interval while skrinking the interval 
//                            when the sampled point is not acceptable
//                            (Figure 5, Neal (2003))
//
// INPUT:
//   interv .......... interval to sample from
//   g_interv[2] ..... log(f(x)) evaluated in the interval limits
//   x0[1] ........... the current point
//   z[1] ............ logarithm of the vertical level defining the slice
//   w[1] ............ estimate of the typical size of slice
//   unimodal[1] ..... indicator whether the distribution at question is unimodal or not
//   doubling ........ 0/1 indicating whether doubling was used to find an interval
//   (*eval) ......... routine to compute g(x) = log(f(x))
//   deval[] ......... double paramaters for (*eval)
//   ieval[] ......... integer parameters for (*eval)
//
// OUTPUT:
//    x1[1] ...... newly sampled point
//
void
ss_shrinkage_sample(double* x1,           double* interv,       double* g_interv,     const double* x0,     
                    const double* z,      const double* w,      const int* doubling,  const int* unimodal,
                    void (*eval)(const double*, double*, const double*, const int*),
                    const double* deval,  const int* ieval)
{
  static double u, gx1;
  static int accept;

  accept = 0;
  do {
    u = unif_rand();
    (*x1) = interv[0] + u*(interv[1] - interv[0]);
    eval(x1, &gx1, deval, ieval);
    if (gx1 > (*z)){
      if (*doubling && (!(*unimodal))){
        ss_accept_doubling(&accept, x0, x1, z, w, interv, eval, deval, ieval);
        if (!accept){    /* do shrinkage */
          if ((*x1) < (*x0)){
            interv[0] = *x1;
            g_interv[0] = gx1;
          }
          else{
            interv[1] = *x1;
            g_interv[1] = gx1;
          }          
        }
      }
      else{
        accept = 1;
      }
    }
    else{   /* do shrinkage */
      if ((*x1) < (*x0)){
        interv[0] = *x1;
        g_interv[0] = gx1;
      }
      else{
        interv[1] = *x1;
        g_interv[1] = gx1;
      }
    }
  } while (!accept);

  return;
}    /* end of the function ss_shrinkage_sample */


// ***** ss_bisection_overrelax: Overrelaxed update using a bisection method
//                               (Figure 10, Neal (2003))
//
// This procedure will work only with UNIMODAL densities!!!
//  to be used also with multimodal densities, variables 'go_left' and 'go_right' would have to be removed
//  and test for acceptance at the end would have to be extended
//
// INPUT:
//   interv .......... interval to sample from
//   x0[1] ........... the current point
//   z[1] ............ logarithm of the vertical level defining the slice
//   w[1] ............ estimate of the typical size of slice
//   a[1] ............ integer limiting endpoint accuracy to 2^{-a}*w
//   doubling ........ 0/1 indicating whether doubling was used to find an interval
//   (*eval) ......... routine to compute g(x) = log(f(x))
//   deval[] ......... double paramaters for (*eval)
//   ieval[] ......... integer parameters for (*eval)
//
// OUTPUT:
//    x1[1] ...... newly sampled point
//
void
ss_bisection_overrelax(double* x1,           double* interv,      const double* x0,     const double* z,      
                       const double* w,      const int* a,        const int* doubling,
                       void (*eval)(const double*, double*, const double*, const int*),
                       const double* deval,  const int* ieval)
{
  static double mid, w_bar, g_mid;
  static int a_bar;
  static bool go_on, go_left, go_right;
  static double interv_hat[2];

  /* When the interval is only of size w, narrow it until the midpoint is inside the slice (or accuracy limit is reached) */
  w_bar = (*w);
  a_bar = (*a);  
  if (interv[1] - interv[0] < 1.1*(*w)){
    go_on = true;
    while (go_on){
      mid = 0.5*(interv[0] + interv[1]);
      eval(&mid, &g_mid, deval, ieval);
      if (a_bar == 0 || g_mid > (*z)) go_on = false;
      else{
        if (*x0 > mid) interv[0] = mid;
        else           interv[1] = mid;
        a_bar--;
        w_bar *= 0.5;
      }
    }
  }

  /* Redefine endpoint locations by bisection, to the specified accuracy */
  interv_hat[0] = interv[0];
  interv_hat[1] = interv[1];
  go_left = go_right = true;
  while (a_bar > 0 && (go_left || go_right)){
    a_bar--;
    w_bar *= 0.5;

    /* Bisection on the left */
    if (go_left){
      mid = interv_hat[0] + w_bar;
      eval(&mid, &g_mid, deval, ieval);
      if (g_mid <= (*z)) interv_hat[0] = mid;
      else               go_left = false;
    }

    /* Bisection on the right */
    if (go_right){
      mid = interv_hat[1] - w_bar;
      eval(&mid, &g_mid, deval, ieval);
      if (g_mid <= (*z)) interv_hat[1] = mid;
      else               go_right = false;
    }    
  }

  /* Find a candidate point by flipping from the current point to the opposite side of (hat{L}, hat{R}), */
  /* then test for acceptability                                                                         */
  *x1 = interv_hat[0] + interv_hat[1] - (*x0);    /* = (L+R)/2 + ((L+R)/2 - x)  */
  eval(x1, &g_mid, deval, ieval);
  if (g_mid <= (*z)) *x1 = *x0;

  return;
}    /* end of the function ss_bisection_overrelax */


// ***** ss_exact_sample: Function to sample a point from the interval while skrinking the interval 
//                        when the sampled point is not acceptable
//            * version for the case when the density to sample from is UNIMODAL
//              and the slice was found exactly (up to some precision)
//
// INPUT:
//   interv .......... interval to sample from
//   g_interv[2] ..... log(f(x)) evaluated in the interval limits
//   x0[1] ........... the current point
//   z[1] ............ logarithm of the vertical level defining the slice
//   (*eval) ......... routine to compute g(x) = log(f(x))
//   deval[] ......... double paramaters for (*eval)
//   ieval[] ......... integer parameters for (*eval)
//
// OUTPUT:
//    x1[1] ...... newly sampled point
//
void
ss_exact_sample(double* x1,           double* interv,       double* g_interv,     const double* x0,
                const double* z,
                void (*eval)(const double*, double*, const double*, const int*),
                const double* deval,  const int* ieval)
{
  static double u, gx1;
  static int accept;

  accept = 0;
  do {
    u = unif_rand();
    (*x1) = interv[0] + u*(interv[1] - interv[0]);
    eval(x1, &gx1, deval, ieval);
    if (gx1 > (*z)){
        accept = 1;
    }
    else{   /* do shrinkage */
      if ((*x1) < (*x0)){
        interv[0] = *x1;
        g_interv[0] = gx1;
      }
      else{
        interv[1] = *x1;
        g_interv[1] = gx1;
      }
    }
  } while (!accept);

  return;
}    /* end of the function ss_exact_sample */


// ***** ss_exact_overrelax: Overrelaxed update while assuming that the slice was found exactly
//                           (up to some specified  precission)
//
//  For UNIMODAL densities only!!!
//
// INPUT:
//   interv .......... interval to sample from
//   g_interv[2] ..... here it is only working space!!!
//   x0[1] ........... the current point
//   z[1] ............ logarithm of the vertical level defining the slice
//   bound[2] ........ bounds for the distribution of x (if any)
//   is_bound[3] ..... is_bound[0:1] = indicators whether there are bounds on left and right
//                     is_bound[2] = indicator whether the distribution at question is unimodal or not
//   (*eval) ......... routine to compute g(x) = log(f(x))
//   deval[] ......... double paramaters for (*eval)
//   ieval[] ......... integer parameters for (*eval)
//
// OUTPUT:
//    x1[1] ...... newly sampled point
//
void
ss_exact_overrelax(double* x1,           double* interv,  const double* x0, const double* z,
                   void (*eval)(const double*, double*, const double*, const int*),
                   const double* deval,  const int* ieval)
{
  static double g_x;

  /* Find a candidate point by flipping from the current point to the opposite side of (hat{L}, hat{R}), */
  /* then test for acceptability                                                                         */
  *x1 = interv[0] + interv[1] - (*x0);    /* = (L+R)/2 + ((L+R)/2 - x)  */
  eval(x1, &g_x, deval, ieval);
  if (g_x <= (*z)) *x1 = *x0;

  return;
}    /* end of the function ss_exact_overrelax */

}  /*** end of the namespace Slice_sampler ***/

#ifdef __cplusplus
}                          /* end of extern "C" */
#endif


