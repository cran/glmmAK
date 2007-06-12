/*** mcmc_poisson.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  13/02/2007
//                      version allowing fixed effects only:                  14/02/2007
//                      version allowing normal random effects:               14/02/2007
//                      version allowing univariate G-spline random effects:  14/02/2007
//                      version allowing bivariate G-spline random effects:   10/04/2007
//
//    ADDITIONAL SMALLER FUNTIONS:
//              13/02/2007:  print_iter_infoP
//
//
// PURPOSE: Poisson regression with random effects
//          * main function to run MCMC
//
/* ********************************************************************************* */

#include "mcmc_poisson.h"

namespace mcmc_poissonA{

extern "C"{

// dirP:                            directory where to store simulated values
//
// Y[n]:                            response vector (observed counts)
// offset[n]:                       offset vector
// n[1]:                            number of observations
// N[1]:                            0, if there are no random effects in the model
//                                  equal to number of clusters if there are random effects in the model
// ni[N]:                           numbers of observations for each cluster
//
// X[p*n]:                          covariates for fixed effects
// p[1]:                            number of fixed effects
//
// Beta[p]:                         INPUT: initial value of the fixed effects                           
//                                  OUTPUT: last sampled value of the fixed effects
// priorBeta[2*p]:                  prior means and prior INVERSE variances for fixed effects (beta parameters)
//     priorBeta[0,...,p-1]         = prior means of fixed effects
//     priorBeta[p,...,2*p-1]       = prior inverse variances of fixed effects
//
// XRE[pRE*n]:                      covariates for random effects
// pRE[1]:                          number of random effects
// RE[nRandom * N]:                 INPUT: initial values of cluster specific random effects
//                                  OUTPUT: last sampled values of cluster specific random effects
// REdist[1]:                       distribution of random effects
//                                  see 'RandomCL.h' for the union of possible values
// hierarCenter[1]:                 1 if hierarchical centering of random effects is used
//                                   => nRandom = pRE
//                                  0 otherwise
//                                   => nRandom = pRE
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
//     store[0]: store expected counts?
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
mcmc_poisson(char **dirP,
             const int *Y,                 const double *offset,           const int *n,                   const int *N,              const int *ni,
             const double *X,              const int *p,
             double *Beta,                 const double *priorBeta,
             const double *XRE,            const int *pRE,
             double *RE,                   const int *REdist,              const int *hierarCenter,
             double *ParRE,                const int *priorREMean_InvVar,  const double *priorParREMean,   const double *priorParREInvVar,
             double *G_lambda_a,           const int *G_dim,               const int *G_ipar,              const double *G_dpar,      
             const double *G_lambdaPrior,  const int *G_apar,              const double *G_aContrast,       
             int *allocRE,
             int *iter,                    const int *nsimul,              const int *store,               const int *mainSimul,
             const int *precision,         int *err)
{
  try{
    const int _k_overrelax = 10;                     /** frequency of overrelaxation (model with G-spline) **/

    const char *fname = "mcmc_poisson.cpp";
    *err = 0;
    const int out_format[2] = {*precision, 1};                /** precision and width for output **/
    std::string dir = *dirP;

    int i;
    double *dP1;
    const int *iP1;
    const double *cdP1;

    GetRNGstate();

    /*** Indicators for model components ***/
    const bool areFixed = (*p > 0);
    const bool areRandom = (*pRE > 0);
    const int nRandom = *pRE;

    /*** Numbers of iterations etc. ***/
    const int niter = nsimul[0];
    const int nthin = nsimul[1];
    const int nwrite = nsimul[2];

    /** Mandatory storing of simulated values **/
    const int storeECount = store[0];
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
    std::string ecountpath    = dir + "/expectcount.sim";
    std::string bpath         = dir + "/b.sim";
    std::string apath         = dir + "/logweight.sim";
    std::string wpath         = dir + "/weight.sim";
    std::string indpath       = dir + "/knotInd.sim";       
    std::string allocpath     = dir + "/alloc.sim";
    std::string lambdapath    = dir + "/lambda.sim";
    std::string gmomentpath   = dir + "/gmoment.sim";
    std::string betaRadjpath  = dir + "/betaRadj.sim";
    std::string varRadjpath   = dir + "/varRadj.sim";

    std::ofstream iterfile, llfile, betaFfile, betaRfile, varRfile, acceptfile, ecountfile, bfile, afile, wfile, indfile, allocfile, lambdafile, gmomentfile, betaRadjfile, varRadjfile;
    In_Output::openFile(iterfile, iterpath, 'a');                     /* append since there are already column labels present provided in R */
    In_Output::openFile(llfile, llpath, 'a');
    if (areFixed) In_Output::openFile(betaFfile, betaFpath, 'a');
    if (areRandom){
      if (*hierarCenter) In_Output::openFile(betaRfile, betaRpath, 'a');
      In_Output::openFile(varRfile, varRpath, 'a');
      In_Output::openFile(bfile, bpath, 'a');
    }
    if (storeAccept) In_Output::openFile(acceptfile, acceptpath, 'a');
    if (storeECount) In_Output::openFile(ecountfile, ecountpath, 'a');
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

    /*** Response matrix (row vector) ***/
    //const MatrixRect<int> y = MatrixRect<int>(1, *n, Y);     /* commented on 14/02/2007 */

    /*** Log(y!)  ***/
    MatrixRect<double> log_y_factor = MatrixRect<double>(1, *n);
    dP1 = log_y_factor.a();
    iP1 = Y;
    for (i = 0; i < *n; i++){
      *dP1 = lgamma1p(double(*iP1));
      dP1++;
      iP1++;
    }

    /*** Offset (will be used to store offset + parts of the model we condition on) ***/
    MatrixRect<double> Offset = MatrixRect<double>(1, *n, offset);

    /*** Expected counts ***/
    MatrixRect<double> ECount     = MatrixRect<double>(1, *n);
    MatrixRect<double> PropECount = MatrixRect<double>(1, *n);

    /*** Object to hold fixed effects ***/
    PredictorPoiss PredictFixed = PredictorPoiss(*p, *n, Beta, priorBeta, X);

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
    RandomPoiss RandomEff = RandomPoiss(*pRE, *N, ni, RE, ParRE, *REdist, 
                                        priorREMean_InvVar, priorParREMean, priorParREInvVar, XRE);

    /*** Object for G-spline ***/
    Gspline2 GsplineRE;
    BiGspline2 BiGsplineRE;
    MatrixRect<int> Alloc;                // allocation indicators
    MatrixRect<int> *MixtureN;            // numbers of observations in each mixture component
    int a_ipars[2];                       // needed by Gspline2::update_alla_lambda1() function

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
        // Rprintf("\nInitial G-spline:\n");
        // BiGsplineRE.print();

        /*** Allocations ***/
        util_BiGspline2::array2allocBi(&Alloc, allocRE,  BiGsplineRE.Kconst(), RandomEff.N());

        /*** Numbers of observations in each component ***/
        MixtureN = new MatrixRect<int>[1];
        if (!MixtureN) throw returnR(fname, "Out of memory (object 'MixtureN')", 99);
        MixtureN[0] = MatrixRect<int>(BiGsplineRE.length(0), BiGsplineRE.length(1));
        util_BiGspline2::alloc2mixtureNBi(MixtureN->a(),  Alloc.aconst(),  BiGsplineRE.total_length(), RandomEff.N());

        a_ipars[0] = RandomEff.N();
        a_ipars[1] = MixtureN[0].aconst()[0];

        // Rprintf("\nInitial allocations:\n");
        // Alloc.printI(0);
        // Rprintf("MixtureN:  ");
        // MixtureN[0].printI(0);
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
      AK_BLAS_LAPACK::c_aPlusb(Offset.a(), offset, RandomEff.etaAconst(), *n);
      PredictFixed.ll_poisson(&LogLik, ECount.a(), Offset.aconst(), Y, log_y_factor.aconst(), 2);
      if (!R_finite(LogLik)) throw returnR("STOP: Initial values lead to infinite likelihood (fixed|random).", 1);
    }
    if (areRandom){
      AK_BLAS_LAPACK::c_aPlusb(Offset.a(), offset, PredictFixed.etaAconst(), *n);
      RandomEff.ll_poisson(&LogLik, ECount.a(), Offset.aconst(), Y, log_y_factor.aconst(),  2);
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
    int overrelax    = 0;    /*** used by G-spline model only ***/
    int writea;              /*** used by bivariate G-spline model only ***/
    int nlambda;             /*** used by bivariate G-spline model only ***/
    int witer;
    Rprintf("Iteration ");    

    double Egspline0, Egspline1, sdRandom0, sdRandom1;
    const double *betaP;
    double *betaRadj;
    double *varRadj;

    switch (*REdist){   /*** REdist switch ***/

    /***** NO RANDOM EFFECTS *****/
    /***** ================= *****/
    case mcmc_Random::_None:
      for (*iter=nullthIter + 1; *iter <= lastIter; (*iter)++){
        for (witer = 1; witer <= nthin; witer++){  /**  thinning cycle  **/
          iterTotal++;                        /* = (iter-1)*nthin + witer                    */
          iterTotalNow++;                     /* = (iter-1)*nthin + witer - nullthIter*nthin */

          /** Update parameters **/
          PredictFixed.update1(Accept.a(), &LogLik, ECount.a(), PropECount.a(), offset, Y, log_y_factor.aconst());
        }    /** end of the thinning cycle  **/

        /**  Write to files  **/
	mcmc_poissonA::print_iter_infoP(writeAll, backs, *iter, nwrite, lastIter);
        if ((*mainSimul) || writeAll){
          writeToFile_1(iter, 1, iterfile, out_format[0], out_format[1]);
          writeToFile_1(&LogLik, 1, llfile, out_format[0], out_format[1]);
	  writeToFile_1(PredictFixed.ThetaAconst(), PredictFixed.nTheta(), betaFfile, out_format[0], out_format[1]);
          if (storeAccept) writeToFile_1(Accept.aconst(), Accept.length(), acceptfile, out_format[0], out_format[1]);
          if (storeECount) writeToFile_1(ECount.aconst(), *n, ecountfile, out_format[0], out_format[1]);    

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
          if (areFixed){
            AK_BLAS_LAPACK::c_aPlusb(Offset.a(), offset, RandomEff.etaAconst(), *n);
            PredictFixed.update1(Accept.a(), &LogLik, ECount.a(), PropECount.a(), Offset.aconst(), Y, log_y_factor.aconst());
          }
          if (*hierarCenter) mcmc_Random::updateMeanRE1(RandomEff.REMean(), RandomEff.ThetaBar(),
				                        RandomEff.PropMean(), RandomEff.U(), RandomEff.I(),
                                                        RandomEff.Thetaconst(), RandomEff.REInvVarconst(), 
                                                        RandomEff.REMeanPriorWMeanconst(), RandomEff.REMeanPriorInvVarconst(), RandomEff.prior_for_REMean(),
                                                        RandomEff.nTheta(), RandomEff.N());
          AK_BLAS_LAPACK::c_aPlusb(Offset.a(), offset, PredictFixed.etaAconst(), *n);
          RandomEff.updateRE1(Accept.a()+1, &LogLik, ECount.a(), PropECount.a(), Offset.aconst(), Y, log_y_factor.aconst());
	  mcmc_Random::updateInvVarRE1(RandomEff.REInvVar(), RandomEff.REInvVarL(),
                                       RandomEff.Theta_REMean(), RandomEff.U(), RandomEff.I(), RandomEff.work_rwishart(),
                                       RandomEff.Thetaconst(), RandomEff.REMeanconst(), 
                                       RandomEff.REInvVarPriorDFconst(), RandomEff.REInvVarPriorInvScaleconst(), RandomEff.prior_for_REInvVar(),
                                       RandomEff.nTheta(), RandomEff.N());
        }    /** end of the thinning cycle  **/

        /**  Write to files  **/
	mcmc_poissonA::print_iter_infoP(writeAll, backs, *iter, nwrite, lastIter);
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
          if (storeECount) writeToFile_1(ECount.aconst(), *n, ecountfile, out_format[0], out_format[1]);    

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

        betaRadj  = (double*) calloc(1, sizeof(double));
        varRadj   = (double*) calloc(1, sizeof(double));
        if (!betaRadj || !varRadj) throw returnR(fname, "Not enough memory available.", 99);

        for (*iter=nullthIter + 1; *iter <= lastIter; (*iter)++){
          for (witer = 1; witer <= nthin; witer++){  /**  thinning cycle  **/
            iterTotal++;                        /* = (iter-1)*nthin + witer                    */
            iterTotalNow++;                     /* = (iter-1)*nthin + witer - nullthIter*nthin */
            overrelax = 1*(iterTotalNow/_k_overrelax != 0);

            /** Update parameters **/
            if (areFixed){
              AK_BLAS_LAPACK::c_aPlusb(Offset.a(), offset, RandomEff.etaAconst(), *n);
              PredictFixed.update1(Accept.a(), &LogLik, ECount.a(), PropECount.a(), Offset.aconst(), Y, log_y_factor.aconst());
            }
	    util_Gspline2::updateAlloc1(Alloc.a(), MixtureN, &GsplineRE, RandomEff.ThetaAconst(), RandomEff.REMeanAconst(), RandomEff.N());  
            AK_BLAS_LAPACK::c_aPlusb(Offset.a(), offset, PredictFixed.etaAconst(), *n);
            RandomEff.updateRE2(Accept.a()+1, &LogLik, ECount.a(), PropECount.a(), &GsplineRE, Offset.aconst(), Y, log_y_factor.aconst(), &Alloc);    /* at this moment only working for GsplineRE._dim=1 */
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
	  mcmc_poissonA::print_iter_infoP(writeAll, backs, *iter, nwrite, lastIter);
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
              betaP       = PredictFixed.ThetaAconst() + (*p - 1);
              betaRadj[0] = betaP[0] + Egspline0;
            }
            writeToFile_1(betaRadj, 1, betaRadjfile, out_format[0], out_format[1]);

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
            if (storeECount) writeToFile_1(ECount.aconst(), *n, ecountfile, out_format[0], out_format[1]);    

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

        betaRadj = (double*) calloc(2, sizeof(double));
        varRadj  = (double*) calloc(3, sizeof(double));
        if (!betaRadj || !varRadj) throw returnR(fname, "Not enough memory available.", 99);

        for (*iter=nullthIter + 1; *iter <= lastIter; (*iter)++){
          for (witer = 1; witer <= nthin; witer++){  /**  thinning cycle  **/
            iterTotal++;                        /* = (iter-1)*nthin + witer                    */
            iterTotalNow++;                     /* = (iter-1)*nthin + witer - nullthIter*nthin */
            overrelax = 1*(iterTotalNow/_k_overrelax != 0);

            /** Update parameters **/
            if (areFixed){
              AK_BLAS_LAPACK::c_aPlusb(Offset.a(), offset, RandomEff.etaAconst(), *n);
              PredictFixed.update1(Accept.a(), &LogLik, ECount.a(), PropECount.a(), Offset.aconst(), Y, log_y_factor.aconst());
            }
	    util_BiGspline2::updateAllocBi(Alloc.a(), MixtureN->a(), &BiGsplineRE, RandomEff.ThetaAconst(), RandomEff.REMeanAconst(), RandomEff.N());  
            AK_BLAS_LAPACK::c_aPlusb(Offset.a(), offset, PredictFixed.etaAconst(), *n);
            RandomEff.updateRE2Bi(Accept.a()+1, &LogLik, ECount.a(), PropECount.a(), &BiGsplineRE, Offset.aconst(), Y, log_y_factor.aconst(), &Alloc);
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
	  mcmc_poissonA::print_iter_infoP(writeAll, backs, *iter, nwrite, lastIter);
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
              betaP       = PredictFixed.ThetaAconst() + (*p - 2);
              betaRadj[0] = betaP[0] + Egspline0;
              betaRadj[1] = betaP[1] + Egspline1;
            }
            writeToFile_1(betaRadj, 2, betaRadjfile, out_format[0], out_format[1]);

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
            if (storeECount) writeToFile_1(ECount.aconst(), *n, ecountfile, out_format[0], out_format[1]);    

            writeAll = 0;
          }
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
        Rprintf("Dimension of the G-spline random effect is %d\n", nRandom);
        throw returnR("Error in mcmc_poisson.cpp: mcmc_poisson(). Only implemented for UNI- and BIVARIATE G-spline random effects", 1);
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
    if (storeECount) ecountfile.close();
    if (storeAccept) acceptfile.close();
    if (areRandom){
      bfile.close();
      varRfile.close();
      if (*hierarCenter) betaRfile.close();
    }
    if (areFixed) betaFfile.close();    
    llfile.close();
    iterfile.close();


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
print_iter_infoP(int &writeAll,  int &backs,  const int &iter,  const int &nwrite,  const int &lastIter)
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

}  /*** end of the namespace mcmc_poissonA ***/

}  /* end of extern "C" */


