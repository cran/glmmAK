#*** cumlogitRE.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              akom@email.cz
##
##           CREATED:  11/09/2006
##   WORKING VERSION:  03/10/2006  (normal random effects)
##                     16/01/2007  (univariate G-spline random effects)
##                     10/04/2007  (bivariate G-spline random effects)
## MAJOR MOIFICATION:  11/10/2006
##                     (internally, regression coefficients are stored in the order 'v' covariates first, 'x' covariates then)
##
## PURPOSE: Cumulative logit model with random effects
##          fitted using Bayesian specification and MCMC
##
## FUNCTIONS: cumlogitRE
## 
#* ********************************************************************************* */

cumlogitRE <- function(y, v, x, vb, xb, cluster,                       
      intcpt.random=FALSE,
      hierar.center=FALSE,                       
      drandom=c("normal", "gspline"),
      C=1,
      logit.order=c("decreasing", "increasing"),
      prior.fixed,
      prior.random,
      prior.gspline,
      init.fixed,
      init.random,
      init.gspline,                 
      nsimul = list(niter=10, nthin=1, nburn=0, nwrite=10),
      store = list(prob=FALSE, b=FALSE, alloc=FALSE, acoef=FALSE),
      dir=getwd(),
      precision=8)
{
  thispackage <- "glmmAK"
  #thispackage <- NULL

  nsimul <- checknsimul.glmmAK(nsimul)
  
  ## Design matrices and check of some input arguments
  DES <- design.cumlogit(y=y, v=v, x=x, vb=vb, xb=xb, cluster=cluster, intcpt.random=intcpt.random, hierar.center=hierar.center, C=C)
  if (hierar.center){
    nFixed <- C*DES$nv + DES$nx
    nRandom <- C*DES$nvb + DES$nxb
  }
  else{
    nFixed <-  C*DES$nv + DES$nx      ## nxb is included in nx, nvb is included in nv
    nRandom <- DES$nvb + DES$nxb
  }  
  
  ## Order of logits
  logit.order <- match.arg(logit.order)  
  if (logit.order == "increasing"){
    stop("logit.order = increasing not yet implemented")
#    yFIT <- C - DES$y
  }
  else{
    yFIT <- DES$y
  }      
  
  ## Prior for fixed effects
  prior.fixed <- prior.fixed.glmmAK(prior.fixed=prior.fixed, nFixed=nFixed)
  PriorFixed <- attr(prior.fixed, "PriorFixed")

  
  ## Distribution of random effects
  if (nRandom){
    drandom <- match.arg(drandom)
    if (drandom == "normal") drandomI <- 1
    else{
      if (drandom == "gspline") drandomI <- 2
    }
  }  
  else{
    drandom <- "none"
    drandomI <- 0
  }  
  
  ## Prior for random effects
  prior.random <- prior.random.glmmAK(prior.random=prior.random, nRandom=nRandom, drandom=drandom, hierar.center=hierar.center)
  PriorMRandom <- attr(prior.random, "PriorMRandom")
  PriorDRandom <- attr(prior.random, "PriorDRandom")
  PriorMDRandom <- attr(prior.random, "PriorMDRandom")  

  
  ## Prior for G-spline
  prior.gspline <- prior.gspline.glmmAK(prior.gspline=prior.gspline, nRandom=nRandom, drandom=drandom)
  Gdim <- attr(prior.gspline, "Dim")
  GiPar <- attr(prior.gspline, "iPar")
  GdPar <- attr(prior.gspline, "dPar")
  GlambdaPrior <- attr(prior.gspline, "lambdaPrior")  
  GaPar <- attr(prior.gspline, "aPar")
  GaContrast <- attr(prior.gspline, "aContrast")

  
  ## Model with fixed effects only fitted using ML
  if (hierar.center){    
    fit.init <- fit.cumlogit(y=DES$y, v=cbind(DES$v, DES$vb), x=cbind(DES$x, DES$xb), C=C, logit.order=logit.order)
    prob.init <- prob00.cumlogit(coef=fit.init$coefficients, v=cbind(DES$v, DES$vb), x=cbind(DES$x, DES$xb), C=C, logit.order=logit.order)
    init.fixedML <- fit.init$coefficients[c(DES$bindv, DES$bindx)]
    init.MrandomML <- fit.init$coefficients[c(DES$bindvb, DES$bindxb)]
    init.DrandomML <- fit.init$expect.vcov[c(DES$bindvb, DES$bindxb), c(DES$bindvb, DES$bindxb)]
  }
  else{
    fit.init <- fit.cumlogit(y=DES$y, v=DES$v, x=DES$x, C=C, logit.order=logit.order)
    prob.init <- prob00.cumlogit(coef=fit.init$coefficients, v=DES$v, x=DES$x, C=C, logit.order=logit.order)
    init.fixedML <- fit.init$coefficients
    init.MrandomML <- rep(0, nRandom)
    if (nRandom){
      if (DES$nxb & DES$nvb) NAAM <- c(colnames(DES$vb), colnames(DES$xb))
      else{
        if (DES$nxb) NAAM <- colnames(DES$xb)
        else         NAAM <- colnames(DES$vb)
      }  
      names(init.MrandomML) <- NAAM
    }  
    init.DrandomML <- fit.init$expect.vcov[DES$ind.logit1, DES$ind.logit1]  ## take only block corresponding to the first logit
  }  
  
  ## Initial values for fixed effects
  init.fixed <- init.fixed.glmmAK(init.fixed=init.fixed, init.fixedML=init.fixedML, nFixed=nFixed)
  
  ## Initial values for random effects and their parameters
  init.random <- init.random.glmmAK(init.random=init.random, init.MrandomML=init.MrandomML, init.DrandomML=init.DrandomML,
                                     nRandom=nRandom, drandom=drandom, N=DES$N, hierar.center=hierar.center)
  InitParRandom <- attr(init.random, "InitParRandom")

  ## Initial values for G-spline parameters
  init.gspline <- init.gspline.glmmAK(init.gspline=init.gspline, prior.gspline=prior.gspline,
                                       init.random=init.random, nRandom=nRandom, drandom=drandom)
  InitLambda.a <- attr(init.gspline, "lambda.a")
    
  ## Write headers to files
  headers.glmmAK(dir=dir, store=store, DES=DES, prob.init=prob.init, init.fixedML=init.fixedML, init.MrandomML=init.MrandomML, hierar.center=hierar.center,
                  drandom=drandom, model="cumlogit", prior.gspline=prior.gspline)

  ## Conversion of store to integer
  storeI <- rep(0, 4)
  names(storeI) <- c("prob", "b", "alloc", "acoef")
  if (!is.null(store$prob))  if (store$prob) storeI["prob"] <- 1
  if (!is.null(store$b))     if (store$b & nRandom) storeI["b"] <- 1
  if (!is.null(store$alloc)) if (store$alloc & nRandom & drandom=="gspline") storeI["alloc"] <- 1
  if (!is.null(store$acoef)) if (store$acoef & nRandom==2 & drandom=="gspline") storeI["acoef"] <- 1
  
  ## Compute quantities to determine the space needed to be allocated
  ##   and numbers of iterations in different phases
  if (nsimul$nburn >= nsimul$niter) nsimul$nburn <- nsimul$niter - 1
  if (nsimul$nburn < 0) nsimul$nburn <- 0
 
  if (nsimul$nburn == 0) nruns <- 1
  else                   nruns <- 2

  nrun <- numeric(2)
  nrun[2] <- nsimul$niter - nsimul$nburn
  nrun[1] <- nsimul$nburn

  nwrite.run <- nrun
  nwrite.run[nsimul$nwrite <= nrun] <- nsimul$nwrite   
  max.nwrite <- max(nwrite.run)

  ## Combine similar parameters into one vector  
  #dims <- c(nobs, as.numeric(doubly))
  nsimul.run1 <- c(nrun[1], nsimul$nthin, nwrite.run[1])
  nsimul.run2 <- c(nrun[2], nsimul$nthin, nwrite.run[2])
  
  cat("Simulation started on                       ", date(), "\n", sep = "")
  if (nrun[1] > 0){
    fit <- .C("mcmc_cumlogit",
                   dir=as.character(dir),
                   Y=as.integer(yFIT),
                   C=as.integer(C),
                   n=as.integer(DES$ny),
                   N=as.integer(DES$N),
                   ni=as.integer(DES$ni),
                   X=as.double(t(DES$x)),
                   V=as.double(t(DES$v)),
                   p=as.integer(DES$nx),
                   q=as.integer(DES$nv),
                   Beta=as.double(init.fixed),
                   priorBeta=as.double(PriorFixed),
                   XRE=as.double(t(DES$xb)),
                   VRE=as.double(t(DES$vb)),
                   pRE=as.integer(DES$nxb),
                   qRE=as.integer(DES$nvb),
                   pRE4C=as.integer(DES$nxb4C),
                   qRE4C=as.integer(DES$nvb4C),
                   bindxb4C=as.integer(DES$bindxb4C - 1),
                   bindvb4C=as.integer(DES$bindvb4C - 1),
                   RE=as.double(t(init.random$b)),
                   REdist=as.integer(drandomI),
                   hierarCenter=as.integer(hierar.center),
                   ParRE=as.double(InitParRandom),
                   priorREMean_InvVar=as.integer(PriorMDRandom),              
                   priorParREMean=as.double(PriorMRandom),
                   priorParREInvVar=as.double(PriorDRandom),
                   G.lambda.a=as.double(InitLambda.a),
                   G.dim=as.integer(Gdim),
                   G.ipar=as.integer(GiPar),
                   G.dpar=as.double(GdPar),
                   G.lambdaPrior=as.double(GlambdaPrior),
                   G.apar=as.integer(GaPar),
                   G.aContrast=as.double(GaContrast),
                   allocRE=as.integer(t(init.gspline$alloc)),
                   iter=as.integer(0),
                   nsimul = as.integer(nsimul.run1),
                   store=as.integer(storeI),
                   mainSimul=as.integer(0),
                   precision=as.integer(precision),              
                   err=as.integer(0),
              PACKAGE=thispackage
         )
    
    cat("\n")
    if (fit$err != 0) stop ("Something went wrong during the simulation.")
    fit$iter <- fit$iter - 1
    cat("Burn-up finished on                         ", date(), "   (iteration ", fit$iter, ")", "\n", sep = "")

    ## Rewrite sampled values by new files
    headers.glmmAK(dir=dir, store=store, DES=DES, prob.init=prob.init, init.fixedML=init.fixedML, init.MrandomML=init.MrandomML, hierar.center=hierar.center,
                    drandom=drandom, model="cumlogit", prior.gspline=prior.gspline)

    
  }
  else{
    fit <- list(Beta=init.fixed, RE=t(init.random$b), ParRE=InitParRandom,
                G.lambda.a=InitLambda.a, allocRE=t(init.gspline$alloc), iter=0)
  }  
  
  ## Main simulation
  fit <- .C("mcmc_cumlogit",
                 dir=as.character(dir),
                 Y=as.integer(yFIT),
                 C=as.integer(C),
                 n=as.integer(DES$ny),
                 N=as.integer(DES$N),
                 ni=as.integer(DES$ni),
                 X=as.double(t(DES$x)),
                 V=as.double(t(DES$v)),
                 p=as.integer(DES$nx),
                 q=as.integer(DES$nv),
                 Beta=as.double(fit$Beta),
                 priorBeta=as.double(PriorFixed),
                 XRE=as.double(t(DES$xb)),
                 VRE=as.double(t(DES$vb)),
                 pRE=as.integer(DES$nxb),
                 qRE=as.integer(DES$nvb),
                 pRE4C=as.integer(DES$nxb4C),
                 qRE4C=as.integer(DES$nvb4C),
                 bindxb4C=as.integer(DES$bindxb4C - 1),
                 bindvb4C=as.integer(DES$bindvb4C - 1),
                 RE=as.double(fit$RE),
                 REdist=as.integer(drandomI),
                 hierarCenter=as.integer(hierar.center),
                 ParRE=as.double(fit$ParRE),
                 priorREMean_InvVar=as.integer(PriorMDRandom),
                 priorParREMean=as.double(PriorMRandom),
                 priorParREInvVar=as.double(PriorDRandom),
                 G.lambda.a=as.double(fit$G.lambda.a),
                 G.dim=as.integer(Gdim),
                 G.ipar=as.integer(GiPar),
                 G.dpar=as.double(GdPar),
                 G.lambdaPrior=as.double(GlambdaPrior),
                 G.apar=as.integer(GaPar),
                 G.aContrast=as.double(GaContrast),
                 allocRE=as.integer(fit$allocRE),            
                 iter=as.integer(fit$iter),
                 nsimul = as.integer(nsimul.run2),
                 store=as.integer(storeI),
                 mainSimul=as.integer(1),
                 precision=as.integer(precision),            
                 err=as.integer(0),
              PACKAGE=thispackage            
       )

  cat("\n")  
  if (fit$err != 0) stop ("Something went wrong during the simulation.")
  fit$iter <- fit$iter - 1  
  cat("Simulation finished on                      ", date(), "   (iteration ", fit$iter, ")", "\n", sep = "")  

  RET <- list(prior.fixed=prior.fixed, prior.random=prior.random, prior.gspline=prior.gspline,
              init.fixed=init.fixed, init.random=init.random, init.gspline=init.gspline,
              design=DES)
  return(RET)  
}

