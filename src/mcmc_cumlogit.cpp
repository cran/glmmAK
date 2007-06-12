/*** mcmc_cumlogit.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  14/09/2006
//                               version allowing fixed effects only:  15/09/2006
//                            version allowing normal random effects:  03/10/2006
//               version allowing univariate G-spline random effects:  16/01/2007
//                version allowing bivariate G-spline random effects:  10/04/2007
//
//    ADDITIONAL SMALLER FUNTIONS:
//              18/10/2006:  print_iter_infoC
//
//
// PURPOSE: Cumulative logit model with random effects
//          * main function to run MCMC
//
/* ********************************************************************************* */

#include "mcmc_cumlogit.h"

namespace mcmc_cumlogit{

extern "C"{

// dirP:                            directory where to store simulated values
//
// Y[n]:                            response vector (values from {0,...,C})
// C[1]:                            number of response categories minus 1
// n[1]:                            number of observations
// N[1]:                            0, if there are no random effects in the model
//                                  equal to number of clusters if there are random effects in the model
// ni[N]:                           numbers of observations for each cluster
//
// X[p*n]:                          covariates for fixed effects whose effect is proportional w.r.t. odds
// V[q*n]:                          covariates for fixed effects whose effect is not proportional w.r.t. odds
// p[1]:                            number of fixed effects whose effect is proportional w.r.t. odds
// q[1]:                            number of fixed effects whose effect is not proportional w.r.t. odds
//
// Beta[C*q+p]:                     INPUT: initial value of the fixed effects                           
//                                  OUTPUT: last sampled value of the fixed effects
// priorBeta[2*(C*q+p)]:            prior means and prior INVERSE variances for fixed effects (beta parameters)
//     priorBeta[0,...,(C*q+p)-1]         = prior means of fixed effects
//     priorBeta[(C*q+p),...,2*(C*q+p)-1] = prior inverse variances of fixed effects
//
// XRE[pRE*n]:                       covariates for random effects whose effect is proportional w.r.t. odds
// VRE[qRE*n]:                       covariates for random effects whose effect is not proportional w.r.t. odds
// pRE[1]:                          number of random effects whose effect is proportional w.r.t. odds
// qRE[1]:                          number of random effects whose effect is not proportional w.r.t. odds
//
// pRE4C[1]:                        used only in G-spline models if !hierarCenter 
//                                  = number of random effects whose effect is originally proportional w.r.t. odds
// qRE4C[1]:                        used only in G-spline models if !hierarCenter
//                                  = number of random effects whose effect is originally proportional w.r.t. odds
//                                  qRE4C + pRE4C should be equal to pRE and qRE should be equal to 0 if !hierarCenter
//        
// bindxb4C[pRE4C]:                  used only in G-spline models if !hierarCenter 
//                                  = indeces (C++ scale starting from 0) of Beta coefficients which correspond to random effect
//                                    whose effect is proportional w.r.t. odds
// bindvb4C[C*qRE4C]:                used only in G-spline models if !hierarCenter 
//                                  = indeces (C++ scale starting from 0) of Beta coefficients which correspond to random effect
//                                    whose effect is not proportional w.r.t. odds
//
// RE[nRandom * N]:                INPUT: initial values of cluster specific random effects
//                                  OUTPUT: last sampled values of cluster specific random effects
// REdist[1]:                       distribution of random effects
//                                  see 'mcmc_Random.h' for the union of possible values
// hierarCenter[1]:                 1 if hierarchical centering of random effects is used
//                                   => nRandom = pRE + C*qRE
//                                  0 otherwise
//                                   => nRandom = pRE + qRE
//
// ParRE[nRandom + LT(nRandom)]:    INPUT: initial values of the mean and inverse variance of random effects
//                                  OUTPUT: last sampled values of the mean and inverse variance of random effects
//     ParRE[0,...,nRandom-1]     = mean of random effects
//     ParRE[nRandom,...]         = lower triangle (in column major order) of the inverse variance of the random effects
//                                  (should be a diagonal matrix when there are G-spline random effects)
//
// priorREMean_InvVar[2]:   types of prior for means and inverse variance of random effects
//     priorREMean_InvVar[0]     = type of prior for means of random effects (see RandomCL.h for possible values)
//     priorREMean_InvVar[1]     = type of prior for inverse variance of random effects (see RandomCL.h for possible values)
//
// priorParREMean[2*nRandom]:          parameters of the prior for the mean of random effects
//                                     can be arbitrary if  priorREMean_InvVar[0] = _Fixed_
//     priorParREMean[0,...,nRandom-1]             = prior means for the means of random effects
//     priorParREMean[pRE+C*qRE,...,2*nRandom-1]   = prior inverse variances for the means of random effects
//
// priorParREInvVar[]:      parameters of the Wishart prior for the inverse variance of random effects
//                          or the uniform prior of the standard deviation of a univariate random effect
//                          or the independent gamma priors
//     if (priorREMean_InvVar[1] = _Fixed)
//         priorParREInvVar                = arbitary
//     if (priorREMean_InvVar[1] = _Wishart)
//         priorParREInvVar[0]                = degrees of freedom
//         priorParREInvVar[1,...]            = inverse scale matrix of the Wishart distribution (lower triangle)
//     if (priorREMean_InvVar[1] = _SDUnif)
//         priorParREInvVar[0,...,nRandom-1]  = upper limits of the uniform distribution for each standard deviation
//         priorParREInvVar[rest]             = arbitrary
//     if (priorREMean_InvVar[1] = _GammaIndep)
//         priorParREInvVar[0,...,nRandom-1]         = shape parameters of the gamma distributions
//         priorParREInvVar[nRandom,...,2*nRandom-1] = inverse scale (rate) parameters of the gamma distribution 
//
// G_lambda_a[]: initial values of lambda's and a coefficients (see argument lambda_a of the parametric Gspline2 constructor)
//
// G_dim[]: number of random effects and numbers of knots on each side of the zero knot for each margin of the G-spline
//           (see argument dimPar of the parametrix Gspline2 constructor)
//
// G_ipar[]:  (see argument penMix_iPar of the parametric Gspline2 constructor)
//
// G_dpar[]: basis standard deviations and knots
//            (see argument penMix_dPar of the parametric Gspline2 constructor)
//
// G_apar[]: (see argument aPar of the parametric Gspline2 constructor)
//
// G_aContrast[]: (see argument aContrast of the parametric Gspline2 constructor)
//                ignored and can be arbitrary in the models with BIVARIATE G-splines
//
// G_lambdaPrior[]: (see argument lambdaPar of the parametric Gspline2 constructor)
//
// allocRE[]: initial allocations to mixture components for random effects
//            on the scale -K,...,K
//
// iter[1]:
//
// nsimul[3]:                       vector with niter, nthin, nwrite where
//     nsimul[0] = niter:  number of simulated values (after discadring thinned ones)
//     nsimul[1] = nthin:  thinning interval
//     nsimul[2] = nwrite: how often is the information concerning the progress of the MCMC printed on the screen
//
// store[?]:                        mandatory storing of some simulated values
//     store[0]: store category probabilities?
//     store[1]: store individual random effects?
//     store[2]: store allocations (model with G-spline)?
//     store[3]: store all log-weights (model with bivariate G-spline)?
//
// mainSimul[1]:                    0/1, if 0 then all parameters are written to files only every nwrite iteration
//                                  (usually equal to 0 during the burn-in)
//
// precision[1]:    precision of the output
//
void
mcmc_cumlogit(char **dirP,
              const int *Y,                 const int *C,                   const int *n,                   const int *N,              const int *ni,
              const double *X,              const double *V,                const int *p,                   const int *q,
              double *Beta,                 const double *priorBeta,
              const double *XRE,            const double *VRE,              const int *pRE,                 const int *qRE,
              const int *pRE4C,             const int *qRE4C,               const int *bindxb4C,            const int *bindvb4C,
              double *RE,                   const int *REdist,              const int *hierarCenter,
              double *ParRE,                const int *priorREMean_InvVar,  const double *priorParREMean,   const double *priorParREInvVar,
              double *G_lambda_a,           const int *G_dim,               const int *G_ipar,              const double *G_dpar,      
              const double *G_lambdaPrior,  const int *G_apar,              const double *G_aContrast,       
              int *allocRE,
              int *iter,                    const int *nsimul,              const int *store,               const int *mainSimul,
              const int *precision,         int *err)
{
  try{
    const int _k_overrelax = 10;                              /** frequency of overrelaxation (model with G-spline) **/

    const char *fname = "mcmc_cumlogit.cpp";
    *err = 0;
    const int out_format[2] = {*precision, 1};                /** precision and width for output **/
    std::string dir = *dirP;

    int i;
    double *dP1;
    const double *cdP1;

    GetRNGstate();

    /*** Indicators for model components ***/
    const bool areFixed = (*p + (*q) > 0);
    const bool areRandom = (*pRE + (*qRE) > 0);
    const int nRandom = *hierarCenter ? (*C * *qRE + *pRE) : (*qRE + *pRE);      

    /*** Numbers of iterations etc. ***/
    const int niter = nsimul[0];
    const int nthin = nsimul[1];
    const int nwrite = nsimul[2];

    /** Mandatory storing of simulated values **/
    const int storeProb   = store[0];
    const int storeb      = store[1];
    const int storeAlloc  = store[2];
    const int storea      = store[3];
    const int storeAccept = 0;

    /*** Open files to store simulated values for writing ***/
    std::string iterpath      = dir + "/iteration.sim";
    std::string llpath        = dir + "/loglik.sim";
    std::string betaFpath     = dir + "/betaF.sim";
    std::string betaRpath     = dir + "/betaR.sim";
    std::string varRpath      = dir + "/varR.sim";
    std::string acceptpath    = dir + "/accept.sim";
    std::string probpath      = dir + "/probability.sim";
    std::string bpath         = dir + "/b.sim";
    std::string apath         = dir + "/logweight.sim";
    std::string wpath         = dir + "/weight.sim";
    std::string indpath       = dir + "/knotInd.sim";       
    std::string allocpath     = dir + "/alloc.sim";
    std::string lambdapath    = dir + "/lambda.sim";
    std::string gmomentpath   = dir + "/gmoment.sim";
    std::string betaRadjpath  = dir + "/betaRadj.sim";
    std::string varRadjpath   = dir + "/varRadj.sim";

    std::ofstream iterfile, llfile, betaFfile, betaRfile, varRfile, acceptfile, probfile, bfile, afile, wfile, indfile, allocfile, lambdafile, gmomentfile, betaRadjfile, varRadjfile;
    In_Output::openFile(iterfile, iterpath, 'a');                     /* append since there are already column labels present provided in R */
    In_Output::openFile(llfile, llpath, 'a');
    if (areFixed) In_Output::openFile(betaFfile, betaFpath, 'a');
    if (areRandom){
      if (*hierarCenter) In_Output::openFile(betaRfile, betaRpath, 'a');
      In_Output::openFile(varRfile, varRpath, 'a');
      In_Output::openFile(bfile, bpath, 'a');
    }
    if (storeAccept) In_Output::openFile(acceptfile, acceptpath, 'a');
    if (storeProb) In_Output::openFile(probfile, probpath, 'a');
    if (*REdist == mcmc_Random::_Gspline){
      In_Output::openFile(afile, apath, 'a');
      In_Output::openFile(allocfile, allocpath, 'a');
      In_Output::openFile(lambdafile, lambdapath, 'a');
      In_Output::openFile(gmomentfile, gmomentpath, 'a');

      In_Output::openFile(betaRadjfile, betaRadjpath, 'a');
      In_Output::openFile(varRadjfile, varRadjpath, 'a');

      if (nRandom == 2){
        In_Output::openFile(wfile, wpath, 'a');
        In_Output::openFile(indfile, indpath, 'a');
      }
    }

    /*** Files written only during the development ***/
    //std::string mixNpath = dir + "/mixtureN.sim";
    //std::ofstream mixNfile;
    //In_Output::openFile(mixNfile, mixNpath, 'o');

    /*** Response matrix (row vector) ***/
    const MatrixRect<int> y = MatrixRect<int>(1, *n, Y);

    /*** Category probabilities ***/
    MatrixRect<double> Prob = MatrixRect<double>(*C+1, *n);
    MatrixRect<double> PropProb = MatrixRect<double>(*C+1, *n);

    /*** Object to hold fixed effects ***/
    PredictorCL PredictFixed = PredictorCL(*p, *q, *n, *C, Beta, priorBeta, X, V);

    /*** Check distribution of random effects ***/
    switch (*REdist){
    case mcmc_Random::_None:
    case mcmc_Random::_Normal:
    case mcmc_Random::_Gspline:
      break;
    default:
      throw returnR(fname, "Unknown REdist argument", 1);
    }

    /*** Object to hold random effects ***/
    RandomCL RandomEff = RandomCL(*pRE, *qRE, *N, ni, *C, RE, ParRE, *REdist, 
                                  priorREMean_InvVar, priorParREMean, priorParREInvVar, XRE, VRE);

    /*** Object for G-spline ***/
    Gspline2 GsplineRE;
    BiGspline2 BiGsplineRE;
    MatrixRect<int> Alloc;                // allocation indicators
    MatrixRect<int> *MixtureN;            // numbers of observations in each mixture component
    int a_ipars[2];                       // needed by Gspline2::update_alla_lambda1() and BiGspline2::update_alla_lambdaBi() functions

    if (*REdist == mcmc_Random::_Gspline){

      switch (nRandom){
      case 1: 
        /*** G-spline object ***/
        GsplineRE = Gspline2(G_dim, G_ipar, G_dpar, G_apar, G_aContrast, G_lambdaPrior, G_lambda_a, RandomEff.REInvVarAconst()); 
        // Rprintf("\nInitial G-spline:\n");
        // GsplineRE.print();

        /*** Allocations ***/
        util_Gspline2::array2alloc(&Alloc, allocRE,  GsplineRE.KAconst(),  GsplineRE.dim(),  RandomEff.N());

        /*** Numbers of observations in each component ***/
        MixtureN = new MatrixRect<int>[GsplineRE.dim()];
        if (!MixtureN) throw returnR(fname, "Out of memory (object 'MixtureN')", 99);
        for (i = 0; i < GsplineRE.dim(); i++) MixtureN[i] = MatrixRect<int>(1, GsplineRE.length(i));
        util_Gspline2::alloc2mixtureN(MixtureN,  Alloc.aconst(),  GsplineRE.KAconst(),  GsplineRE.dim(),  RandomEff.N());

        a_ipars[0] = RandomEff.N();
        a_ipars[1] = MixtureN[0].aconst()[0];

        // Rprintf("\nInitial allocations:\n");
        // Alloc.printI(0);
        // for (i = 0; i < GsplineRE.dim(); i++){
        //   Rprintf("MixtureN[%d]:  ", i);
        //   MixtureN[i].printI(0);
        // }
        break;

      case 2:
        /*** G-spline object ***/
        BiGsplineRE = BiGspline2(G_dim, G_ipar, G_dpar, G_apar, G_lambdaPrior, G_lambda_a, RandomEff.REInvVarAconst()); 
        //Rprintf("\nInitial G-spline:\n");
        //BiGsplineRE.print();

        /*** Allocations ***/
        util_BiGspline2::array2allocBi(&Alloc, allocRE,  BiGsplineRE.Kconst(), RandomEff.N());

        /*** Numbers of observations in each component ***/
        MixtureN = new MatrixRect<int>[1];
        if (!MixtureN) throw returnR(fname, "Out of memory (object 'MixtureN')", 99);
        MixtureN[0] = MatrixRect<int>(BiGsplineRE.length(0), BiGsplineRE.length(1));
        util_BiGspline2::alloc2mixtureNBi(MixtureN->a(),  Alloc.aconst(),  BiGsplineRE.total_length(), RandomEff.N());

        a_ipars[0] = RandomEff.N();
        a_ipars[1] = MixtureN[0].aconst()[0];

        //Rprintf("\nInitial allocations:\n");
        //Alloc.printI(0);
        //Rprintf("MixtureN:  ");
        //MixtureN[0].printI(0);
        break;

      default:
        REprintf("nRandom = %d\n", nRandom);
        throw returnR(fname, "G-spline distribution of random effects only implemented uni- and bivariately.", 1);
      }
    }

    /*** Object to store acceptance indicators ***/
    int nAccept = 1;
    if (areRandom) nAccept += *N;
    MatrixRect<int> Accept = MatrixRect<int>(1, nAccept);

    /*** Log-likelihood, also check starting values ***/
    int anyZero;
    double LogLik;
    if (areFixed){
      PredictFixed.ll_cumlogit2(&LogLik, Prob.a(), &anyZero, RandomEff.etaAconst(), y.aconst(), 1, 2, ll_cumlogitFS2);
      if (!R_finite(LogLik)) throw returnR("STOP: Initial values lead to infinite likelihood (fixed|random).", 1);
    }
    if (areRandom){
      RandomEff.ll_cumlogit2(&LogLik, Prob.a(), &anyZero, PredictFixed.etaAconst(), y.aconst(), 1, 2, ll_cumlogitFS2);
      if (!R_finite(LogLik)) throw returnR("STOP: Initial values lead to infinite likelihood (random|fixed).", 1);
    }

    /** Main simulation: **/
      /* nullthIter ..... index for initial values (usually 0)                             */
      /* iter ........... counter of iterations                                            */
      /* iterTotal ...... total number of iterations done previously (thinnig included)    */
      /* iterTotalNow ... total number of iterations done with this function call          */
      /*                  -> this number would be used to compute acceptance rates         */  
      /* backs .......... how many times the carriage must be returned to print            */
      /*                  the next iteration number?                                       */
      /* writeAll ....... if equal to 1 then all values needed to restart the MCMC         */
      /*                  from the current point are written to files                      */
      /*                                                                                   */
    int nullthIter   = *iter;
    int lastIter     = nullthIter + niter;
    int iterTotal    = nullthIter*nthin;
    int iterTotalNow = 0;
    int backs        = 0;
    int writeAll     = 0;
    int overrelax    = 0;    /*** used by G-spline model only           ***/
    int writea;              /*** used by bivariate G-spline model only ***/
    int nlambda;             /*** used by bivariate G-spline model only ***/
    int witer;
    Rprintf("Iteration ");    

    double Egspline0, Egspline1, sdRandom0, sdRandom1;
    const double *betaP;
    double *betaRadj;
    double *varRadj;
    const int *indRandom;       /** to point to indeces of betaF corresponding to G-spline random effect when nRandom = 1 and !hierarCenter             **/
    const int *indRandom0;      /** to point to indeces of betaF corresponding to G-spline random effect (margin 0) when nRandom = 2 and !hierarCenter  **/
    const int *indRandom1;      /** to point to indeces of betaF corresponding to G-spline random effect (margin 1) when nRandom = 2 and !hierarCenter  **/
    int lindRandom, lindRandom0, lindRandom1, VXcase;

    switch (*REdist){   /*** REdist switch ***/

    /***** NO RANDOM EFFECTS *****/
    /***** ================= *****/
    case mcmc_Random::_None:
      for (*iter=nullthIter + 1; *iter <= lastIter; (*iter)++){
        for (witer = 1; witer <= nthin; witer++){  /**  thinning cycle  **/
          iterTotal++;                        /* = (iter-1)*nthin + witer                    */
          iterTotalNow++;                     /* = (iter-1)*nthin + witer - nullthIter*nthin */

          /** Update parameters **/
          PredictFixed.update1(Accept.a(), &LogLik, &Prob, &PropProb, RandomEff.etaconst(), &y);
        }    /** end of the thinning cycle  **/

        /**  Write to files  **/
	mcmc_cumlogitA::print_iter_infoC(writeAll, backs, *iter, nwrite, lastIter);
        if ((*mainSimul) || writeAll){
          writeToFile_1(iter, 1, iterfile, out_format[0], out_format[1]);
          writeToFile_1(&LogLik, 1, llfile, out_format[0], out_format[1]);
	  writeToFile_1(PredictFixed.ThetaAconst(), PredictFixed.nTheta(), betaFfile, out_format[0], out_format[1]);
          if (storeAccept) writeToFile_1(Accept.aconst(), Accept.length(), acceptfile, out_format[0], out_format[1]);
          if (storeProb) writeToFile_1(Prob.aconst(), Prob.length(), probfile, out_format[0], out_format[1]);    

          writeAll = 0;
        }
      }    /** end of the main loop over iter **/

      /*** Copy last sampled values to the output ***/
      PredictFixed.Thetaconst()->mat2array(Beta, 0);
      break;

    /***** NORMAL RANDOM EFFECTS *****/
    /***** ===================== *****/
    case mcmc_Random::_Normal:
      for (*iter=nullthIter + 1; *iter <= lastIter; (*iter)++){
        for (witer = 1; witer <= nthin; witer++){  /**  thinning cycle  **/
          iterTotal++;                        /* = (iter-1)*nthin + witer                    */
          iterTotalNow++;                     /* = (iter-1)*nthin + witer - nullthIter*nthin */

          /** Update parameters **/
          if (areFixed) PredictFixed.update1(Accept.a(), &LogLik, &Prob, &PropProb, RandomEff.etaconst(), &y);
          if (*hierarCenter) mcmc_Random::updateMeanRE1(RandomEff.REMean(), RandomEff.ThetaBar(),
				                        RandomEff.PropMean(), RandomEff.U(), RandomEff.I(),
                                                        RandomEff.Thetaconst(), RandomEff.REInvVarconst(), 
                                                        RandomEff.REMeanPriorWMeanconst(), RandomEff.REMeanPriorInvVarconst(), RandomEff.prior_for_REMean(),
                                                        RandomEff.nTheta(), RandomEff.N());
          RandomEff.updateRE1(Accept.a()+1, &LogLik, &Prob, &PropProb, PredictFixed.etaconst(), &y);
	  mcmc_Random::updateInvVarRE1(RandomEff.REInvVar(), RandomEff.REInvVarL(),
                                       RandomEff.Theta_REMean(), RandomEff.U(), RandomEff.I(), RandomEff.work_rwishart(),
                                       RandomEff.Thetaconst(), RandomEff.REMeanconst(), 
                                       RandomEff.REInvVarPriorDFconst(), RandomEff.REInvVarPriorInvScaleconst(), RandomEff.prior_for_REInvVar(),
                                       RandomEff.nTheta(), RandomEff.N());
        }    /** end of the thinning cycle  **/

        /**  Write to files  **/
	mcmc_cumlogitA::print_iter_infoC(writeAll, backs, *iter, nwrite, lastIter);
        if ((*mainSimul) || writeAll){
          writeToFile_1(iter, 1, iterfile, out_format[0], out_format[1]);
          writeToFile_1(&LogLik, 1, llfile, out_format[0], out_format[1]);
	  if (areFixed) writeToFile_1(PredictFixed.ThetaAconst(), PredictFixed.nTheta(), betaFfile, out_format[0], out_format[1]);
          if (*hierarCenter) writeToFile_1(RandomEff.REMeanAconst(), RandomEff.nTheta(), betaRfile, out_format[0], out_format[1]);

          RandomEff.invert_REInvVarFromL();    /** Compute inversion of InvVar **/
          writeToFile2_1(RandomEff.REVarAconst(), RandomEff.lInvVar(),
                         RandomEff.REInvVarAconst(), RandomEff.lInvVar(), varRfile, out_format[0], out_format[1]);

          if (storeb || writeAll) writeToFile_1(RandomEff.ThetaAconst(), RandomEff.N_nTheta(), bfile, out_format[0], out_format[1]);
          if (storeAccept) writeToFile_1(Accept.aconst(), Accept.length(), acceptfile, out_format[0], out_format[1]);
          if (storeProb) writeToFile_1(Prob.aconst(), Prob.length(), probfile, out_format[0], out_format[1]);    

          writeAll = 0;
        }
      }    /** end of the main loop over iter **/

      /*** Copy last sampled values to the output ***/
      if (areFixed) PredictFixed.Thetaconst()->mat2array(Beta, 0);
      RandomEff.Thetaconst()->mat2array(RE, 0);
      if (*hierarCenter) RandomEff.REMeanconst()->mat2array(ParRE, 0);
      RandomEff.REInvVarconst()->mat2array(ParRE+RandomEff.nTheta(), 0);     
      break;

    /***** G-SPLINE RANDOM EFFECTS *****/
    /***** ======================= *****/
    case mcmc_Random::_Gspline:

      switch (nRandom){    /*** Dimension switch ***/

      /***** UNIVARIATE G-SPLINE RANDOM EFFECTS *****/
      /***** ---------------------------------- *****/
      case 1:
        if (nRandom != GsplineRE.dim()){
          REprintf("nRandom=%d,  GsplineRE.dim()=%d\n", nRandom, GsplineRE.dim());
          throw returnR(fname, "Inconsistent nRandom and GsplineRE.dim(). Contact the author.", 1);
        }
        
        if (*hierarCenter){
          lindRandom = 1;
	}
        else{
          if (*qRE4C){
            lindRandom = *C;
            indRandom  = bindvb4C;            
          }
          else{
            if (*pRE4C){
              lindRandom = 1;
              indRandom  = bindxb4C;
            }
       	    else{
              REprintf("qRE4C=%d,  pRE4C=%d,  C=%d,  nRandom=%d\n", *qRE4C, *pRE4C, *C, nRandom);
              throw returnR(fname, "Inconsistent qRE4C/pRE4C and nRandom. Contact the author.", 1);
            }
          }
        }
        betaRadj  = (double*) calloc(lindRandom, sizeof(double));
        varRadj   = (double*) calloc(1, sizeof(double));
        if (!betaRadj || !varRadj) throw returnR(fname, "Not enough memory available.", 99);

        for (*iter=nullthIter + 1; *iter <= lastIter; (*iter)++){
          for (witer = 1; witer <= nthin; witer++){  /**  thinning cycle  **/
            iterTotal++;                        /* = (iter-1)*nthin + witer                    */
            iterTotalNow++;                     /* = (iter-1)*nthin + witer - nullthIter*nthin */
            overrelax = 1*(iterTotalNow/_k_overrelax != 0);

            /** Update parameters **/
            if (areFixed) PredictFixed.update1(Accept.a(), &LogLik, &Prob, &PropProb, RandomEff.etaconst(), &y);
	    util_Gspline2::updateAlloc1(Alloc.a(), MixtureN, &GsplineRE, RandomEff.ThetaAconst(), RandomEff.REMeanAconst(), RandomEff.N());  
            RandomEff.updateRE2(Accept.a()+1, &LogLik, &Prob, &PropProb, &GsplineRE, PredictFixed.etaconst(), &y, &Alloc);    /* at this moment only working for GsplineRE._dim=1 */
            if (*hierarCenter) mcmc_Random::updateMeanRE2(RandomEff.REMean(), RandomEff.ThetaBar(),
                                                          RandomEff.Thetaconst(),
                                                          RandomEff.REMeanPriorWMeanconst(), RandomEff.REMeanPriorInvVarconst(), RandomEff.prior_for_REMean(),
                                                          &GsplineRE, &Alloc,
                                                          RandomEff.nTheta(), RandomEff.N());                                                        
	    mcmc_Random::updateInvVarRE2(RandomEff.REInvVar(), RandomEff.REInvVarL(), RandomEff.REVar(),                           /* It updates also _REVar. */
                                         &GsplineRE, RandomEff.work_invVarSlice(),
                                         RandomEff.Thetaconst(), RandomEff.REMeanconst(),
                                         RandomEff.REInvVarPriorDFconst(), RandomEff.REInvVarPriorInvScaleconst(), RandomEff.prior_for_REInvVar(),
                                         &Alloc, 
                                         RandomEff.nTheta(), RandomEff.N(), overrelax);
            GsplineRE.update_alla_lambda1(a_ipars, MixtureN->aconst(), &overrelax);
          }    /** end of the thinning cycle  **/

          /**  Write to files  **/
	  mcmc_cumlogitA::print_iter_infoC(writeAll, backs, *iter, nwrite, lastIter);
          if ((*mainSimul) || writeAll){
            writeToFile_1(iter,    1, iterfile, out_format[0], out_format[1]);
            writeToFile_1(&LogLik, 1, llfile,   out_format[0], out_format[1]);
	    if (areFixed) writeToFile_1(PredictFixed.ThetaAconst(), PredictFixed.nTheta(), betaFfile, out_format[0], out_format[1]);

            Egspline0 = sqrt(RandomEff.REVarAconst()[0]) * GsplineRE.gMeanAconst()[0];
            if (*hierarCenter){
              writeToFile_1(RandomEff.REMeanAconst(), RandomEff.nTheta(), betaRfile, out_format[0], out_format[1]);

              betaRadj[0] = RandomEff.REMeanAconst()[0] + Egspline0;
            }
            else{
              betaP = PredictFixed.ThetaAconst();
              for (i = 0; i < lindRandom; i++) betaRadj[i] = betaP[indRandom[i]] + Egspline0;
            }
            writeToFile_1(betaRadj, lindRandom, betaRadjfile, out_format[0], out_format[1]);

            varRadj[0] = RandomEff.REVarAconst()[0] * GsplineRE.gVarAconst()[0];
            writeToFile_1(varRadj, 1, varRadjfile, out_format[0], out_format[1]);

            writeToFile2_1_diag(RandomEff.REVarAconst(), RandomEff.nTheta(),
                                RandomEff.REInvVarAconst(), RandomEff.nTheta(), varRfile, out_format[0], out_format[1]);
            GsplineRE.write_logweights2file(afile, out_format[0], out_format[1]);
            GsplineRE.gMoments();
            writeToFile2_1(GsplineRE.gMeanAconst(), GsplineRE.dim(), GsplineRE.gVarAconst(), GsplineRE.dim(), gmomentfile, out_format[0], out_format[1]);
            writeToFile_1(GsplineRE.lambdaAconst(), GsplineRE.dim(), lambdafile, out_format[0], out_format[1]);
            if (storeb || writeAll) writeToFile_1(RandomEff.ThetaAconst(), RandomEff.N_nTheta(), bfile, out_format[0], out_format[1]);
            if (storeAlloc || writeAll) writeToFile_1(Alloc.aconst(), RandomEff.N_nTheta(), allocfile, out_format[0], out_format[1]);
            if (storeAccept) writeToFile_1(Accept.aconst(), Accept.length(), acceptfile, out_format[0], out_format[1]);
            if (storeProb) writeToFile_1(Prob.aconst(), Prob.length(), probfile, out_format[0], out_format[1]);    

            writeAll = 0;
          }
        }    /** end of the main loop over iter **/

        /*** Copy last sampled values to the output ***/
        if (areFixed) PredictFixed.Thetaconst()->mat2array(Beta, 0);
        RandomEff.Thetaconst()->mat2array(RE, 0);
        if (*hierarCenter) RandomEff.REMeanconst()->mat2array(ParRE, 0);
        RandomEff.REInvVarconst()->mat2array(ParRE+RandomEff.nTheta(), 0);

        G_lambda_a[0] = GsplineRE.lambdaAconst()[0];
        dP1  = G_lambda_a + 1;
        cdP1 = GsplineRE.aconst()[0].aconst();
        for (i = 0; i < GsplineRE.length(0); i++){
          *dP1 = *cdP1;
          dP1++;
          cdP1++;
        }
        util_Gspline2::alloc2array(allocRE, &Alloc, GsplineRE.KAconst(),  GsplineRE.dim(),  RandomEff.N());                

        free(betaRadj);
        free(varRadj);

        break;

      /***** BIVARIATE G-SPLINE RANDOM EFFECTS *****/
      /***** --------------------------------- *****/
      case 2:
        if (nRandom != BiGspline2A::_dim){
          REprintf("nRandom=%d,  BiGspline2A::_dim=%d\n", nRandom, BiGspline2A::_dim);
          throw returnR(fname, "Inconsistent nRandom and BiGspline2::_dim. Contact the author.", 1);
        }
        nlambda = BiGsplineRE.equal_lambda() ? 1 : BiGspline2A::_dim;

        if (*hierarCenter){
          lindRandom0 = 1;
          lindRandom1 = 1;
        }
        else{
          if (*qRE4C && *pRE4C){
            VXcase = mcmc_cumlogitA::VX;
            lindRandom0 = *C;
            lindRandom1 = 1;
            indRandom0  = bindvb4C;
            indRandom1  = bindxb4C;
          }
          else{
            if (*qRE4C){
              VXcase = mcmc_cumlogitA::VV;
              lindRandom0 = *C;
              lindRandom1 = *C;
              indRandom0  = bindvb4C;
              indRandom1  = bindvb4C + 1;
            }
            else{
              if (*pRE4C){
                VXcase = mcmc_cumlogitA::XX;
                lindRandom0 = 1;
                lindRandom1 = 1;
                indRandom0  = bindxb4C;
                indRandom1  = bindxb4C + 1; 
              }
              else{
                REprintf("qRE4C=%d,  pRE4C=%d,  C=%d,  nRandom=%d\n", *qRE4C, *pRE4C, *C, nRandom);
                throw returnR(fname, "Inconsistent qRE4C/pRE4C and nRandom. Contact the author.", 1);
              }
            }
          }
        }
        lindRandom = lindRandom0 + lindRandom1;

        betaRadj = (double*) calloc(lindRandom, sizeof(double));
        varRadj  = (double*) calloc(3, sizeof(double));
        if (!betaRadj || !varRadj) throw returnR(fname, "Not enough memory available.", 99);

        for (*iter=nullthIter + 1; *iter <= lastIter; (*iter)++){
          for (witer = 1; witer <= nthin; witer++){  /**  thinning cycle  **/
            iterTotal++;                        /* = (iter-1)*nthin + witer                    */
            iterTotalNow++;                     /* = (iter-1)*nthin + witer - nullthIter*nthin */
            overrelax = 1*(iterTotalNow/_k_overrelax != 0);

            /** Update parameters **/
            if (areFixed) PredictFixed.update1(Accept.a(), &LogLik, &Prob, &PropProb, RandomEff.etaconst(), &y);
	    util_BiGspline2::updateAllocBi(Alloc.a(), MixtureN->a(), &BiGsplineRE, RandomEff.ThetaAconst(), RandomEff.REMeanAconst(), RandomEff.N());  
            RandomEff.updateRE2Bi(Accept.a()+1, &LogLik, &Prob, &PropProb, &BiGsplineRE, PredictFixed.etaconst(), &y, &Alloc);
            if (*hierarCenter) mcmc_Random::updateMeanRE2Bi(RandomEff.REMean(), RandomEff.ThetaBar(),
                                                            RandomEff.PropMean(), RandomEff.U(), RandomEff.I(),
                                                            RandomEff.Thetaconst(),
                                                            RandomEff.REMeanPriorWMeanconst(), RandomEff.REMeanPriorInvVarconst(), RandomEff.prior_for_REMean(),
                                                            &BiGsplineRE, &Alloc,
                                                            RandomEff.nTheta(), RandomEff.N());                                                        
	    mcmc_Random::updateInvVarRE2Bi(RandomEff.REInvVar(), RandomEff.REInvVarL(), RandomEff.REVar(),                           /* It updates also _REVar. */
                                           &BiGsplineRE, RandomEff.work_invVarSlice(),
                                           RandomEff.Thetaconst(), RandomEff.REMeanconst(),
                                           RandomEff.REInvVarPriorDFconst(), RandomEff.REInvVarPriorInvScaleconst(), RandomEff.prior_for_REInvVar(),
                                           &Alloc, 
                                           RandomEff.nTheta(), RandomEff.N(), overrelax);
            BiGsplineRE.update_alla_lambdaBi(a_ipars, MixtureN->aconst(), &overrelax);
          }    /** end of the thinning cycle  **/

          /**  Write to files  **/
	  mcmc_cumlogitA::print_iter_infoC(writeAll, backs, *iter, nwrite, lastIter);
          if ((*mainSimul) || writeAll){
            writeToFile_1(iter,    1, iterfile, out_format[0], out_format[1]);
            writeToFile_1(&LogLik, 1, llfile,   out_format[0], out_format[1]);
	    if (areFixed) writeToFile_1(PredictFixed.ThetaAconst(), PredictFixed.nTheta(), betaFfile, out_format[0], out_format[1]);

            sdRandom0 = sqrt(RandomEff.REVarAconst()[0]); 
            sdRandom1 = sqrt(RandomEff.REVarAconst()[2]); 
            Egspline0 = sdRandom0 * BiGsplineRE.gMeanconst()[0];
            Egspline1 = sdRandom1 * BiGsplineRE.gMeanconst()[1];
            if (*hierarCenter){
              writeToFile_1(RandomEff.REMeanAconst(), RandomEff.nTheta(), betaRfile, out_format[0], out_format[1]);

              betaRadj[0] = RandomEff.REMeanAconst()[0] + Egspline0;
              betaRadj[1] = RandomEff.REMeanAconst()[1] + Egspline1;
            }
            else{
              betaP = PredictFixed.ThetaAconst();
              switch (VXcase){
              case mcmc_cumlogitA::VX:
                for (i = 0; i < lindRandom0; i++) betaRadj[i] = betaP[indRandom0[i]] + Egspline0;
                betaRadj[lindRandom0] = betaP[*indRandom1] + Egspline1;
                break;

              case mcmc_cumlogitA::VV:
                for (i = 0; i < lindRandom0; i++) betaRadj[2*i]   = betaP[indRandom0[2*i]] + Egspline0;
                for (i = 0; i < lindRandom1; i++) betaRadj[2*i+1] = betaP[indRandom1[2*i]] + Egspline1;
                break;

              case mcmc_cumlogitA::XX:
                betaRadj[0] = betaP[*indRandom0] + Egspline0;
                betaRadj[1] = betaP[*indRandom1] + Egspline1;
                break;
              }
            }
            writeToFile_1(betaRadj, lindRandom, betaRadjfile, out_format[0], out_format[1]);

            varRadj[0] = RandomEff.REVarAconst()[0] * BiGsplineRE.gVarconst()[0];                                      /** variance in margin 0   **/
            varRadj[1] = sdRandom0 * sdRandom1 * BiGsplineRE.gVarconst()[1];                                           /** covariance             **/
            varRadj[2] = RandomEff.REVarAconst()[2] * BiGsplineRE.gVarconst()[2];                                      /** variance in margin 1   **/
            writeToFile_1(varRadj, 3, varRadjfile, out_format[0], out_format[1]);

            writeToFile2_1_diag(RandomEff.REVarAconst(), RandomEff.nTheta(),
                                RandomEff.REInvVarAconst(), RandomEff.nTheta(), varRfile, out_format[0], out_format[1]);
            writea = 1*(storea || writeAll);
            BiGsplineRE.write_logweights2file(afile, wfile, indfile, &writea, out_format[0], out_format[1]);
            BiGsplineRE.gMoments();
            writeToFile2_1(BiGsplineRE.gMeanconst(), BiGspline2A::_dim, BiGsplineRE.gVarconst(), BiGspline2A::_LTdim, gmomentfile, out_format[0], out_format[1]);
            writeToFile_1(BiGsplineRE.lambdaconst(), nlambda, lambdafile, out_format[0], out_format[1]);
            if (storeb || writeAll) writeToFile_1(RandomEff.ThetaAconst(), RandomEff.N_nTheta(), bfile, out_format[0], out_format[1]);
            if (storeAlloc || writeAll) writeToFile_1(Alloc.aconst(), RandomEff.N(), allocfile, out_format[0], out_format[1]);
            if (storeAccept) writeToFile_1(Accept.aconst(), Accept.length(), acceptfile, out_format[0], out_format[1]);
            if (storeProb) writeToFile_1(Prob.aconst(), Prob.length(), probfile, out_format[0], out_format[1]);    

            writeAll = 0;
          }

          //writeToFile_1(MixtureN->aconst(), MixtureN->length(), mixNfile, out_format[0], out_format[1]);
        }    /** end of the main loop over iter **/

        /*** Copy last sampled values to the output ***/
        if (areFixed) PredictFixed.Thetaconst()->mat2array(Beta, 0);
        RandomEff.Thetaconst()->mat2array(RE, 0);
        if (*hierarCenter) RandomEff.REMeanconst()->mat2array(ParRE, 0);
        RandomEff.REInvVarconst()->mat2array(ParRE+RandomEff.nTheta(), 0);

        dP1  = G_lambda_a;
        cdP1 = BiGsplineRE.lambdaconst();
        for (i = 0; i < BiGspline2A::_dim; i++){
          *dP1 = *cdP1;
          dP1++;
          cdP1++;       
        }
        cdP1 = BiGsplineRE.aAconst();
        for (i = 0; i < BiGsplineRE.total_length(); i++){
          *dP1 = *cdP1;
          dP1++;
          cdP1++;       
        }
        util_BiGspline2::alloc2arrayBi(allocRE, &Alloc, BiGsplineRE.Kconst(), RandomEff.N());

        free(betaRadj);
        free(varRadj);

        break;

      default:
        Rprintf("Dimension of the G-spline random effect is %d\n", GsplineRE.dim());
        throw returnR("Error in mcmc_cumlogit.cpp: mcmc_cumlogit(). Only implemented for UNI- and BIVARIATE G-spline random effects", 1);
      }   /*** End of dimension switch ***/

      break;
    }   /*** End of REdist switch ***/

    /**************************************************************************/

    /*** Close files ***/
    if (*REdist == mcmc_Random::_Gspline){
      gmomentfile.close();
      lambdafile.close();
      allocfile.close();
      afile.close();

      betaRadjfile.close();
      varRadjfile.close();

      if (nRandom == 2){
        wfile.close();
        indfile.close();
      }
    }
    if (storeProb) probfile.close();
    if (storeAccept) acceptfile.close();
    if (areRandom){
      bfile.close();
      varRfile.close();
      if (*hierarCenter) betaRfile.close();
    }
    if (areFixed) betaFfile.close();    
    llfile.close();
    iterfile.close();

    /*** Close files written only during the development ***/
    //mixNfile.close();

    /*** Cleaning ***/
    if (*REdist == mcmc_Random::_Gspline){
      delete [] MixtureN;
    }

    PutRNGstate();
    return;
  }  
  catch(returnR rr){
    *err = rr.errflag();
    PutRNGstate();
    return;
  }
}

void
print_iter_infoC(int &writeAll,  int &backs,  const int &iter,  const int &nwrite,  const int &lastIter)
{
  static int i;

  if (!(iter % nwrite) || iter == lastIter){
    writeAll = 1;
    for (i = 0; i < backs; i++) Rprintf("\b");
    Rprintf("%d", iter);
    backs = int(log10(double(iter))) + 1;
  }

  return;
}

}  /* end of extern "C" */

}  /*** end of the namespace mcmc_cumlogitA ***/

