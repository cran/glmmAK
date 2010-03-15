#*** glmmAK.help.R ***/
##
##     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
##              akom@email.cz
##
##         CREATED:  26/09/2006
##
## PURPOSE: GALMM models
##          * low level functions
##          * not to be called by the users
##
## FUNCTIONS: checknsimul.glmmAK   (06/02/2007)
##
##            prior.fixed.glmmAK
##            init.fixed.glmmAK
##            prior.random.glmmAK
##            init.random.glmmAK
##            prior.gspline.glmmAK
##            init.gspline.glmmAK
##
##            headers.glmmAK
## 
#* ********************************************************************************* */

### checknsimul.glmmAK
### ----------------------
checknsimul.glmmAK <- function(nsimul)
{
  if(length(nsimul) == 0) innsimul <- "arnost"
  else                    innsimul <- names(nsimul)

  tmp <- match("niter", innsimul, nomatch=NA)
  if(is.na(tmp)) stop("nsimul$niter must be given")
  if (nsimul$niter <= 0) stop("nsimul$niter must be positive")
  
  tmp <- match("nburn", innsimul, nomatch=NA)
  if(is.na(tmp)) stop("nsimul$nburn must be given")
  if (nsimul$nburn < 0) stop("nsimul$nburn must be non-negative")
  if (nsimul$nburn > nsimul$niter) stop("nsimul$nburn must not be higher than nsimul$niter")

  tmp <- match("nwrite", innsimul, nomatch=NA)
  if(is.na(tmp)) nsimul$nwrite <- nsimul$niter
  
  tmp <- match("nthin", innsimul, nomatch=NA)
  if(is.na(tmp)) nsimul$nthin <- 1
  
  return(nsimul)  
}

#* ********************************************************************************* */
prior.fixed.glmmAK <- function(prior.fixed, nFixed)
{
  if (nFixed){
    if (missing(prior.fixed)) stop("prior.fixed must be given")
    
    if (is.null(prior.fixed$mean)) stop("Prior mean of the fixed effects not given")
    if (length(prior.fixed$mean) == 1) prior.fixed$mean <- rep(prior.fixed$mean, nFixed)
    if (length(prior.fixed$mean) != nFixed) stop("Incorrect prior.fixed$mean specified")
    
    if (is.null(prior.fixed$var)) stop("Prior variance of the fixed effects not given")    
    if (length(prior.fixed$var) == 1) prior.fixed$var <- rep(prior.fixed$var, nFixed)    
    if (length(prior.fixed$var) != nFixed) stop("Incorrect prior.fixed$var specified")
    if (any(prior.fixed$var <= 0)) stop("prior.fixed$var must be positive")   
    fixed.ivar <- 1/prior.fixed$var
    
    PriorFixed <- c(prior.fixed$mean, fixed.ivar)
    names(PriorFixed) <- paste(rep(c("mean", "invVar"), each=nFixed), rep(1:nFixed, 2), sep="")
  }
  else{
    prior.fixed <- list(mean=0, var=0)
    PriorFixed <- 0
  }

  attr(prior.fixed, "PriorFixed") <- PriorFixed
  return(prior.fixed)  
}

#* ********************************************************************************* */
init.fixed.glmmAK <- function(init.fixed, init.fixedML, nFixed)
{
  if (!nFixed){
    init.fixed <- 0
    return(init.fixed)
  }  
  
  if (missing(init.fixed)){
    init.fixed <- as.numeric(init.fixedML)
  }
  else{
    if (length(init.fixed) != nFixed) stop("length(init.fixed) not consistent with the number of fixed effects")
    init.fixed <- as.numeric(init.fixed)
  }  

  return(init.fixed)  
}  

#* ********************************************************************************* */
prior.random.glmmAK <- function(prior.random, nRandom, drandom=c("normal", "gspline"), hierar.center)
{
  if (nRandom){
    if (missing(prior.random)) stop("prior.random must be given")

    drandom <- match.arg(drandom)

    #####  Prior for means of random effects  #####
    #####  =================================  #####
    if (hierar.center){
      if (is.null(prior.random$Mdistrib)) prior.random$Mdistrib <- "normal"
      random.Mdistrib <- pmatch(prior.random$Mdistrib, c("fixed", "normal"), nomatch=0) - 1
      if (random.Mdistrib < 0){
        stop("Incorrect prior.random$Mdistrib specified")
      }  
      if (random.Mdistrib == 0){
        prior.random$Mdistrib <- "fixed"
        prior.random$Mmean <- rep(0, nRandom)
        prior.random$Mvar <- rep(0, nRandom)
        PriorMRandom <- rep(0, 2*nRandom)
        names(PriorMRandom) <- paste(rep(c("mean", "invVar"), each=nRandom), rep(1:nRandom, 2), sep="")
      }  
      else{
        prior.random$Mdistrib <- "normal"
        if (is.null(prior.random$Mmean)) stop("Prior mean of the means of random effects not given")
        if (length(prior.random$Mmean) == 1) prior.random$Mmean <- rep(prior.random$Mmean, nRandom)
        if (length(prior.random$Mmean) != nRandom) stop("Incorrect prior.random$Mmean specified")
    
        if (is.null(prior.random$Mvar)) stop("Prior variance of the means of random effects not given")
        if (length(prior.random$Mvar) == 1) prior.random$Mvar <- rep(prior.random$Mvar, nRandom)
        if (length(prior.random$Mvar) != nRandom) stop("Incorrect prior.random$Mvar specified")
        if (any(prior.random$Mvar <= 0)) stop("prior.random$Mvar must be positive")
        random.iMvar <- 1/prior.random$Mvar
    
        PriorMRandom <- c(prior.random$Mmean, random.iMvar)
        names(PriorMRandom) <- paste(rep(c("mean", "invVar"), each=nRandom), rep(1:nRandom, 2), sep="")        
      }  
    }
    else{
      prior.random$Mdistrib <- "fixed"
      random.Mdistrib <- 0
      prior.random$Mmean <- rep(0, nRandom)
      prior.random$Mvar <- rep(0, nRandom)
      PriorMRandom <- rep(0, 2*nRandom)
    }

    
    #####  Prior for covariance matrix of NORMALLY distributed random effects  #####
    #####  ==================================================================  #####       
    if (drandom == "normal"){      

      if (is.null(prior.random$Ddistrib)) prior.random$Ddistrib <- "wishart"
      random.Ddistrib <- pmatch(prior.random$Ddistrib, c("fixed", "wishart", "sduniform", "gamma"), nomatch=0) - 1
      if (random.Ddistrib < 0) stop("Incorrect prior.random$Ddistrib specified")

      ### Fixed covariance matrix
      if (random.Ddistrib == 0){
        prior.random$Ddistrib <- "fixed"
        prior.random$Ddf <- 0
        prior.random$Dshape <- 0
        prior.random$DinvScale <- matrix(0, nrow=nRandom, ncol=nRandom)
        prior.random$Dupper <- 0
        random.Discale <- prior.random$DinvScale[lower.tri(prior.random$DinvScale, diag=TRUE)]
        PriorDRandom <- c(prior.random$Ddf, random.Discale)
        names(PriorDRandom) <- c("df", paste("invScale", 1:length(random.Discale), sep=""))        
      }

      ### Wishart covariance matrix
      if (random.Ddistrib == 1){
        prior.random$Ddistrib <- "wishart"
        prior.random$Dupper <- 0
        prior.random$Dshape <- 0        
        
        if (is.null(prior.random$Ddf)) stop("Prior degrees of freedom for the Wishart prior of var(b) not given")
        if (length(prior.random$Ddf) != 1) stop("Incorrect prior.random$Ddf specified")
        if (prior.random$Ddf <= nRandom - 1) stop(paste("prior.random$Ddf must be > ", nRandom - 1, sep=""))

        if (is.null(prior.random$DinvScale)) stop("Prior inverse scale for the Wishart prior of var(b) not given")
        if (nRandom == 1){
          if (length(prior.random$DinvScale) != 1) stop("Incorrect prior.random$DinvScale specified")
          if (prior.random$DinvScale <= 0) stop("prior.random$DinvScale must be positive")
          random.Discale <- prior.random$DinvScale        
        }
        else{            
          if (length(prior.random$DinvScale) == 1){
            if (prior.random$DinvScale <= 0) stop("prior.random$DinvScale must be positive")
            prior.random$DinvScale = diag(prior.random$DinvScale, nrow=nRandom, ncol=nRandom)
          }
          if (!is.matrix(prior.random$DinvScale)) stop("prior.random$DinvScale must be a matrix")
          if (nrow(prior.random$DinvScale) != nRandom | ncol(prior.random$DinvScale) != nRandom)
            stop(paste("prior.random$DinvScale must have ", nRandom, " rows and columns", sep=""))
          tmp <- t(prior.random$DinvScale)
          prior.random$DinvScale[upper.tri(prior.random$DinvScale)] <- tmp[upper.tri(tmp)]
          CHOL <- try(chol(prior.random$DinvScale), silent=TRUE)
          if (class(CHOL) == "try-error") stop("prior.random$DinvScale not positive definite")
          random.Discale <- prior.random$DinvScale[lower.tri(prior.random$DinvScale, diag=TRUE)]
        }

        PriorDRandom <- c(prior.random$Ddf, random.Discale)
        names(PriorDRandom) <- c("df", paste("invScale", 1:length(random.Discale), sep=""))
      }

      ### SD uniform covariance matrix
      if (random.Ddistrib == 2){        
        prior.random$Ddistrib <- "sduniform"
        prior.random$Ddf <- 0
        prior.random$Dshape <- 0        
        prior.random$DinvScale <- 0

        if (nRandom > 1) stop(paste("prior.random$Ddistrib may not be equal to ", dQuote("sduniform"), " when there are multivariate normal random effects in the model", sep=""))
        if (is.null(prior.random$Dupper)) stop("Prior upper limit of the uniform prior of sd(b) not given")
        if (length(prior.random$Dupper) == 1) prior.random$Dupper <- rep(prior.random$Dupper, nRandom)
        if (length(prior.random$Dupper) != nRandom) stop("Incorrect prior.random$Dupper specified")
        if (any(prior.random$Dupper <= 0)) stop("prior.random$Dupper must be > 0")
  
        PriorDRandom <- prior.random$Dupper
        names(PriorDRandom) <- paste("upper", 1:nRandom)        
      }

      ### Gamma indep. covariance matrix
      if (random.Ddistrib == 3){
        if (nRandom == 1){
          prior.random$Ddistrib <- "gamma"
          prior.random$Ddf <- 0
          prior.random$Dupper <- 0

          if (is.null(prior.random$Dshape)) stop("Shape parameter of the gamma prior for ivar(b) not given")
          if (length(prior.random$Dshape) != nRandom) stop("Incorrect prior.random$Dshape specified")
          if (any(prior.random$Dshape <= 0)) stop("prior.random$Dshape must be > 0")

          if (is.null(prior.random$DinvScale)) stop("Rate (inverse scale) parameters of the gamma prior for ivar(b) not given")
          if (length(prior.random$DinvScale) != nRandom) stop("Incorrect prior.random$DinvScale specified")
          if (any(prior.random$DinvScale <= 0)) stop("prior.random$DinvScale must be > 0")        

          PriorDRandom <- c(2*prior.random$Dshape, 2*prior.random$DinvScale)     ## Wishart parametrization
          names(PriorDRandom) <- c("df", "invScale1")
          random.Ddistrib <- 1                    
        } else{   
          stop(paste("prior.random$Ddistrib = ", dQuote("gamma"), " not allowed when drandom = ", dQuote("normal"), "and random effects are not univariate", sep=""))
        }   
      }        
    }    ## end of if (drandom == "normal")

    
    #####  Prior for covariance matrix of G-SPLINE distributed random effects  #####
    #####  ==================================================================  #####           
    if (drandom == "gspline"){

      if (is.null(prior.random$Ddistrib)) prior.random$Ddistrib <- "gamma"
      random.Ddistrib <- pmatch(prior.random$Ddistrib, c("fixed", "wishart", "sduniform", "gamma"), nomatch=0) - 1
      if (random.Ddistrib < 0) stop("Incorrect prior.random$Ddistrib specified")

      ### Fixed covariance matrix
      if (random.Ddistrib == 0){
        prior.random$Ddistrib <- "fixed"
        prior.random$Ddf <- 0
        prior.random$Dshape <- 0        
        prior.random$DinvScale <- matrix(0, nrow=nRandom, ncol=nRandom)
        prior.random$Dupper <- 0
        random.Discale <- prior.random$DinvScale[lower.tri(prior.random$DinvScale, diag=TRUE)]
        PriorDRandom <- c(prior.random$Ddf, random.Discale)
        names(PriorDRandom) <- c("df", paste("invScale", 1:length(random.Discale), sep=""))        
      }

      ### Wishart covariance matrix
      if (random.Ddistrib == 1){
        stop(paste("prior.random$Ddistrib = ", dQuote("wishart"), " not allowed when drandom = ", dQuote("gspline"), sep=""))
      }

      ### SD uniform covariance matrix
      if (random.Ddistrib == 2){        
        prior.random$Ddistrib <- "sduniform"
        prior.random$Ddf <- 0
        prior.random$Dshape <- 0        
        prior.random$DinvScale <- 0

        if (is.null(prior.random$Dupper)) stop("Prior upper limit of the uniform prior of sd(b) not given")
        if (length(prior.random$Dupper) == 1) prior.random$Dupper <- rep(prior.random$Dupper, nRandom)
        if (length(prior.random$Dupper) != nRandom) stop("Incorrect prior.random$Dupper specified")
        if (any(prior.random$Dupper <= 0)) stop("prior.random$Dupper must be > 0")
  
        PriorDRandom <- prior.random$Dupper
        names(PriorDRandom) <- paste("upper", 1:nRandom)
      }  

      ### Gamma indep. covariance matrix
      if (random.Ddistrib == 3){
        prior.random$Ddistrib <- "gamma"
        prior.random$Ddf <- 0
        prior.random$Dupper <- 0        

        if (is.null(prior.random$Dshape)) stop("Shape parameters of gamma priors for var(b) not given")
        if (length(prior.random$Dshape) == 1) prior.random$Dshape <- rep(prior.random$Dshape, nRandom)
        if (length(prior.random$Dshape) != nRandom) stop("Incorrect prior.random$Dshape specified")
        if (any(prior.random$Dshape <= 0)) stop("prior.random$Dshape must be > 0")
        
        if (is.null(prior.random$DinvScale)) stop("Rate (inverse scale) parameters of gamma priors for var(b) not given")
        if (length(prior.random$DinvScale) == 1) prior.random$DinvScale <- rep(prior.random$DinvScale, nRandom)
        if (length(prior.random$DinvScale) != nRandom) stop("Incorrect prior.random$DinvScale specified")
        if (any(prior.random$DinvScale <= 0)) stop("prior.random$DinvScale must be > 0")        

        PriorDRandom <- c(prior.random$Dshape, prior.random$DinvScale)
        names(PriorDRandom) <- paste(rep(c("shape", "invScale"), each=nRandom), rep(1:nRandom, 2), sep="")
      }  
    }    ## end of if (drandom == "gspline")
  
  }  ## end of if (nRandom)
  else{  
    prior.random <- list(Mmean=0, Mvar=0, Mdistrib="fixed", Ddf=0, Dshape=0, DinvScale=0, Dupper=0, Ddistrib="fixed")
    PriorMRandom <- PriorDRandom <- random.Mdistrib <- random.Ddistrib <- 0
  }

  attr(prior.random, "PriorMRandom") <- PriorMRandom
  attr(prior.random, "PriorDRandom") <- PriorDRandom
  attr(prior.random, "PriorMDRandom") <- c(random.Mdistrib, random.Ddistrib)
  return(prior.random)
}

#* ********************************************************************************* */
init.random.glmmAK <- function(init.random, init.MrandomML, init.DrandomML, nRandom, drandom=c("normal", "gspline"), N, hierar.center)
{
  if (nRandom){
    if (hierar.center) Initb <- matrix(rep(init.MrandomML, N), nrow=N, ncol=nRandom, byrow=TRUE)
    else               Initb <- matrix(rep(0, N*nRandom), nrow=N, ncol=nRandom)

    drandom <- match.arg(drandom)
    
    if (missing(init.random)) init.random <- list()
    
    if (length(init.random)) ininit <- names(init.random)
    else{
      ininit <- "arnost"
      init.random <- list()
    }  

    tmp <- match("b", ininit, nomatch=NA)
    if (is.na(tmp)) init.random$b <- Initb

    tmp <- match("mean", ininit, nomatch=NA)
    if (is.na(tmp)){
      if (hierar.center) init.random$mean <- init.MrandomML
      else               init.random$mean <- rep(0, nRandom)
    }  

    tmp <- match("var", ininit, nomatch=NA)
    if (is.na(tmp)){
      if (drandom == "normal") init.random$var <- init.DrandomML
      if (drandom == "gspline"){
        if (nRandom == 1) init.random$var <- init.DrandomML
        else              init.random$var <- diag(init.DrandomML)
      }  
    }  
    
    
    #####  Initial means of random effects  #####
    #####  ===============================  #####
    if (length(init.random$mean) != nRandom) stop(paste("init.random$mean must have the length ", nRandom, sep=""))
    if (!hierar.center & any(init.random$mean != 0)) warning("Do you really want to have the mean of random effects different from zero?")


    #####  Initial variances of random effects  #####
    #####  ===================================  #####
    if (nRandom == 1){
      if (length(init.random$var) != 1) stop("init.random$var must have length 1")
      if (init.random$var <= 0) stop("init.random$var must be positive")
      InitInvDRandom <- 1/init.random$var
    }
    else{
      if (drandom == "normal"){      
        if (!is.matrix(init.random$var)) stop("init.random$var must be matrix")
        if (nrow(init.random$var) != nRandom | ncol(init.random$var) != nRandom) stop(paste("init.random$var must have ", nRandom, " rows and columns", sep=""))
        tmp <- t(init.random$var)
        init.random$var[upper.tri(init.random$var)] <- tmp[upper.tri(tmp)]
        CHOL <- try(chol(init.random$var), silent=TRUE)
        if (class(CHOL) == "try-error") stop("init.random$var not positive definite")
        InitInvDRandom <- chol2inv(CHOL)
        InitInvDRandom <- InitInvDRandom[lower.tri(InitInvDRandom, diag=TRUE)]
      }
      if (drandom == "gspline"){
        if (length(init.random$var) != nRandom) stop(paste("init.random$var must be a vector of length ", nRandom, sep=""))
        if (any(init.random$var <= 0)) stop("Negative values in init.random$var supplied")
        InitInvDRandom <- diag(nRandom)      ## nRandom is > 1 here
        diag(InitInvDRandom) <- 1/init.random$var
        InitInvDRandom <- InitInvDRandom[lower.tri(InitInvDRandom, diag=TRUE)]        
      }  
    }

    InitParRandom <- c(init.random$mean, InitInvDRandom)

    #####  Initial values of random effects  #####
    #####  ================================  #####    
    if (nRandom == 1){
      if (length(init.random$b) != N) stop("Incorrect init.random$b supplied")
    }
    else{
      if (is.data.frame(init.random$b)) init.random$b <- as.matrix(init.random$b)
      if (!is.matrix(init.random$b)) stop("init.random$b must be matrix")
      if (nrow(init.random$b) != N | ncol(init.random$b) != nRandom) stop(paste("init.random$b must have ", N, " rows and ", nRandom, " columns", sep=""))
    }        
  }  ## end of if (nRandom)
  
  else{
    init.random <- list(b=0, mean=0, var=0)
    InitParRandom <- 0
  }  
  
  attr(init.random, "InitParRandom") <- InitParRandom
  return(init.random)  
}  

#* ********************************************************************************* */
## \item{simplified}{if TRUE, most of the components of prior.gspline are not checked. The simplified version
##    is, e.g., used by predictive functions when we do not need all the prior information}
##
prior.gspline.glmmAK <- function(prior.gspline, nRandom, drandom=c("normal", "gspline"), simplified=FALSE)
{
  if (drandom != "none") drandom <- match.arg(drandom)
  if (nRandom <= 0 | drandom == "normal"){
    prior.gspline <- list(K=0, delta=0, sigma=0, CARorder=0, neighbor.system="uniCAR",
                          Ldistrib="fixed", Lequal=TRUE, Lshape=0, Linvscale=0, Lupper=0,
                          Aident="reference", Areference=0, Acontrast=0, AtypeUpdate="slice")
    Dim <- 0
    names(Dim) <- "nRandom"
    dPar <- 0
    names(dPar) <- "sigma"
    lambdaPrior <- rep(0, 2)
    names(lambdaPrior) <- c("Lshape", "LinvScale")    
    iPar <- c(0, 1, 0, 0, 0)
    names(iPar) <- c("AtypeUpdate", "Lequal", "Ldistrib", "CARorder", "neighbor.system")
    aPar <- rep(0, 2)
    names(aPar) <- c("Aident", "Areference")
    Acontrast <- 0
    names(Acontrast) <- "aContrast"
  }  
  else{    
    if (length(prior.gspline)) inprior <- names(prior.gspline)
    else{
      inprior <- "arnost"
      prior.gspline <- list()
    }  

    ### K ###  
    tmp <- match("K", inprior, nomatch=NA)
    if (is.na(tmp)) prior.gspline$K <- rep(15, nRandom)
    else{
      if (length(prior.gspline$K) == 1) prior.gspline$K <- rep(prior.gspline$K, nRandom)
      if (length(prior.gspline$K) != nRandom) stop(paste("prior.gspline$K must be a vector of length ", nRandom, sep=""))
      if (any(prior.gspline$K < 0)) stop("All components of prior.gspline$K must be non-negative")
    }  
    Dim <- c(nRandom, prior.gspline$K)
    names(Dim) <- c("nRandom", paste("K", 1:nRandom, sep=""))

    ### knots, sigma ###
    tmp <- match("delta", inprior, nomatch=NA)
    if (is.na(tmp)) prior.gspline$delta <- rep(0.3, nRandom)
    else{
      if (length(prior.gspline$delta) == 1) prior.gspline$delta <- rep(prior.gspline$delta, nRandom)
      if (length(prior.gspline$delta) != nRandom) stop(paste("prior.gspline$delta must be a vector of length ", nRandom, sep=""))
      if (any(prior.gspline$delta <= 0)) stop("All components of prior.gspline$delta must be positive")
    }  
    knots <- numeric()
    for (i in 1:nRandom){
      Knots <- ((-prior.gspline$K[i]):prior.gspline$K[i])*prior.gspline$delta[i]
      names(Knots) <- paste("mu", i, ".", (-prior.gspline$K[i]):prior.gspline$K[i], sep="")
      knots <- c(knots, Knots)
    }  
    
    tmp <- match("sigma", inprior, nomatch=NA)
    if (is.na(tmp)) prior.gspline$sigma <- (2/3)*prior.gspline$delta
    else{
      if (length(prior.gspline$sigma) == 1) prior.gspline$sigma <- rep(prior.gspline$sigma, nRandom)
      if (length(prior.gspline$sigma) != nRandom) stop(paste("prior.gspline$sigma must be a vector of length ", nRandom, sep=""))
      if (any(prior.gspline$sigma <= 0)) stop("All components of prior.gspline$sigma must be positive")
    }
    Sigma <- prior.gspline$sigma
    names(Sigma) <- paste("sigma", 1:nRandom, sep="")
    dPar <- c(Sigma, knots)

    ### Lequal, Ldistrib, Lshape, LinvScale, Lupper ###    
    if (simplified){
      prior.gspline$Lequal    <- TRUE
      prior.gspline$Ldistrib  <- "fixed"
      prior.gspline$Lshape    <- rep(0, nRandom)
      prior.gspline$LinvScale <- rep(0, nRandom)      
      prior.gspline$Lupper    <- rep(0, nRandom)      
    }
    else{
      tmp <- match("Lequal", inprior, nomatch=NA)
      if (is.na(tmp)) prior.gspline$Lequal <- FALSE
      else{
        if (length(prior.gspline$Lequal) != 1) stop("prior.gspline$Lequal must be of length 1")
        if (!is.logical(prior.gspline$Lequal)) stop("prior.gspline$Lequal must be logical")
      }
      
      tmp <- match("Ldistrib", inprior, nomatch=NA)
      if (is.na(tmp)) prior.gspline$Ldistrib <- "gamma"    
      gspline.Ldistrib <- pmatch(prior.gspline$Ldistrib, c("fixed", "gamma", "sduniform"), nomatch=0) - 1
      if (gspline.Ldistrib < 0) stop("Incorrect prior.gspline$Ldistrib specified")
      if (gspline.Ldistrib == 0){
        prior.gspline$Ldistrib <- "fixed"
        prior.gspline$Lshape <- rep(0, nRandom)
        prior.gspline$LinvScale <- rep(0, nRandom)      
        prior.gspline$Lupper <- rep(0, nRandom)
        lambdaPrior <- rep(0, 2*nRandom)
      }
      else{
        if (gspline.Ldistrib == 1){
          prior.gspline$Ldistrib <- "gamma"
          tmp <- match("Lshape", inprior, nomatch=NA)
          if (is.na(tmp)) stop("prior.gspline$Lshape must be given")
          if (length(prior.gspline$Lshape) == 1) prior.gspline$Lshape <- rep(prior.gspline$Lshape, nRandom)
          if (prior.gspline$Lequal) prior.gspline$Lshape <- rep(prior.gspline$Lshape[1], nRandom)        
          if (length(prior.gspline$Lshape) != nRandom) stop(paste("prior.gspline$Lshape must be a vector of length ", nRandom, sep=""))
          if (any(prior.gspline$Lshape <= 0)) stop("All components of prior.gspline$Lshape must be positive")
        
          tmp <- match("LinvScale", inprior, nomatch=NA)        
          if (is.na(tmp)) stop("prior.gspline$LinvScale must be given")
          if (length(prior.gspline$LinvScale) == 1) prior.gspline$LinvScale <- rep(prior.gspline$LinvScale, nRandom)
          if (prior.gspline$Lequal) prior.gspline$LinvScale <- rep(prior.gspline$LinvScale[1], nRandom)        
          if (length(prior.gspline$LinvScale) != nRandom) stop(paste("prior.gspline$LinvScale must be a vector of length ", nRandom, sep=""))
          if (any(prior.gspline$LinvScale <= 0)) stop("All components of prior.gspline$LinvScale must be positive")        

          prior.gspline$Lupper <- rep(0, nRandom)
          lambdaPrior <- c(prior.gspline$Lshape, prior.gspline$LinvScale)
          names(lambdaPrior) <- c(paste("Lshape", 1:nRandom, sep=""), paste("LinvScale", 1:nRandom, sep=""))        
        }
        else{
          prior.gspline$Ldistrib <- "sduniform"
          tmp <- match("Lupper", inprior, nomatch=NA)
          if (is.na(tmp)) stop("prior.gspline$Lupper must be given")
          if (length(prior.gspline$Lupper) == 1) prior.gspline$Lupper <- rep(prior.gspline$Lupper, nRandom)
          if (prior.gspline$Lequal) prior.gspline$Lupper <- rep(prior.gspline$Lupper[1], nRandom)        
          if (length(prior.gspline$Lupper) != nRandom) stop(paste("prior.gspline$Lupper must be a vector of length ", nRandom, sep=""))
          if (any(prior.gspline$Lupper <= 0)) stop("All components of prior.gspline$Lupper must be positive")

          prior.gspline$Lshape <- rep(0, nRandom)
          prior.gspline$LinvScale <- rep(0, nRandom)
          lambdaPrior <- c(prior.gspline$Lupper, rep(0, nRandom))
          names(lambdaPrior) <- c(paste("Lupper", 1:nRandom, sep=""), paste("null", 1:nRandom, sep=""))        
        }  
      }
    }  

    ### CARorder ###
    if (simplified){
      prior.gspline$CARorder <- rep(3, nRandom)
    }
    else{    
      tmp <- match("CARorder", inprior, nomatch=NA)
      if (is.na(tmp)) prior.gspline$CARorder <- rep(3, nRandom)
      if (length(prior.gspline$CARorder) == 1) prior.gspline$CARorder <- rep(prior.gspline$CARorder, nRandom)
      if (length(prior.gspline$CARorder) != nRandom) stop(paste("prior.gspline$CARorder must be a vector of length ", nRandom, sep=""))
      if (any(prior.gspline$CARorder < 0)) stop("prior.gspline$CARorder must be non-negative")
    }  

    ### neighbor.system
    if (simplified){
      prior.gspline$neighbor.system <- "uniCAR"
    }
    else{
      tmp <- match("neighbor.system", inprior, nomatch=NA)
      if (is.na(tmp)) prior.gspline$neighbor.system <- "uniCAR"
      gspline.neighbor.system <- pmatch(prior.gspline$neighbor.system, c("uniCAR", "eight.neighbors", "twelve.neighbors"), nomatch=0) - 1
      if (gspline.neighbor.system < 0) stop("Incorrect prior.gspline$neighbor.system specified")
      if (gspline.neighbor.system != 0){
        prior.gspline$Lequal    <- TRUE
      }
    }     
    
    ### AtypeUpdate ###
    if (simplified){
      prior.gspline$AtypeUpdate <- "slice"
    }
    else{   
      tmp <- match("AtypeUpdate", inprior, nomatch=NA)
      if (is.na(tmp)) prior.gspline$AtypeUpdate <- "slice"
      gspline.AtypeUpdate <- pmatch(prior.gspline$AtypeUpdate, c("slice", "ars.quantile", "ars.mode", "block"), nomatch=0) - 1
      if (gspline.AtypeUpdate < 0) stop("Incorrect prior.gspline$AtypeUpdate specified")    
    
      iPar <- c(gspline.AtypeUpdate, 1*prior.gspline$Lequal, gspline.Ldistrib, prior.gspline$CARorder, gspline.neighbor.system)
      names(iPar) <- c("AtypeUpdate", "Lequal", "Ldistrib", paste("CARorder", 1:nRandom, sep=""), "neighbor.system")
    }
      
    ### Aident, Areference ###
    if (simplified){
      prior.gspline$Aident <- "reference"
      prior.gspline$Areference <- rep(0, nRandom)
    }  
    else{    
      tmp <- match("Aident", inprior, nomatch=NA)
      if (is.na(tmp)) prior.gspline$Aident <- "reference"

      prior.gspline$Acontrast <- list()
      Acontrast <- numeric()
      Aident <- pmatch(prior.gspline$Aident, c("mean", "reference"), nomatch=0) - 1
      names(Aident) <- "Aident"
      if (Aident < 0) stop("Incorrect prior.gspline$Aident specified")
      if (Aident == 0){
        if (gspline.AtypeUpdate == 3) warning("I am not sure how the program works when mean contrast and block update of a coefficients are used together")
      
        prior.gspline$Aident <- "mean"
        prior.gspline$Areference <- rep(0, nRandom)
        names(prior.gspline$Areference) <- paste("aref", 1:nRandom, sep="")
        Areference <- prior.gspline$Areference
        for (i in 1:nRandom){
          prior.gspline$Acontrast[[i]] <- rep(1, 2*prior.gspline$K[i] + 1)/(2*prior.gspline$K[i] + 1)
          names(prior.gspline$Acontrast[[i]]) <- paste("a", i, ".", (-prior.gspline$K[i]):prior.gspline$K[i], sep="")
          Acontrast <- c(Acontrast, prior.gspline$Acontrast[[i]])
        }
      }     
      else{
        prior.gspline$Aident <- "reference"
        tmp <- match("Areference", inprior, nomatch=NA)
        if (is.na(tmp)) prior.gspline$Areference <- rep(0, nRandom)
        if (length(prior.gspline$Areference) == 1) prior.gspline$Areference <- rep(prior.gspline$Areference, nRandom)
        if (length(prior.gspline$Areference) != nRandom) stop(paste("prior.gspline$Areference must be a vector of length ", nRandom, sep=""))
        if (any(abs(prior.gspline$Areference) > prior.gspline$K)) stop("prior.gspline$Areference outside the range of knots indeces")
        Areference <- prior.gspline$Areference
        names(Areference) <- paste("aref", 1:nRandom, sep="")      
        for (i in 1:nRandom){      
          prior.gspline$Acontrast[[i]] <- rep(0, 2*prior.gspline$K[i] + 1)
          prior.gspline$Acontrast[[i]][prior.gspline$Areference[i] + prior.gspline$K[i] + 1] <- 1
          names(prior.gspline$Acontrast[[i]]) <- paste("a", i, ".", (-prior.gspline$K[i]):prior.gspline$K[i], sep="")
          Acontrast <- c(Acontrast, prior.gspline$Acontrast[[i]])        
        }
      }
      aPar <- c(Aident, Areference)
    }
  }  

  attr(prior.gspline, "Dim") <- Dim
  attr(prior.gspline, "dPar") <- dPar  
  if (!simplified){
    attr(prior.gspline, "lambdaPrior") <- lambdaPrior
    attr(prior.gspline, "iPar") <- iPar
    attr(prior.gspline, "aPar") <- aPar
    attr(prior.gspline, "aContrast") <- Acontrast
  }  
  return(prior.gspline)    
}  


#* ********************************************************************************* */
init.gspline.glmmAK <- function(init.gspline, prior.gspline, init.random, nRandom, drandom=c("normal", "gspline"))
{
  require(smoothSurv)
  
  if (drandom != "none") drandom <- match.arg(drandom)
  if (nRandom <= 0 | drandom == "normal"){
    init.gspline <- list(lambda=0, weights=list(0), acoef=list(0), alloc=0)
    lambda.a <- c(0, 0)
    names(lambda.a) <- c("lambda", "a")
  }
  else{
    if (missing(init.gspline)) init.gspline <- list()

    if (length(init.gspline)) ininit <- names(init.gspline)
    else{
      ininit <- "arnost"
      init.gspline <- list()
    }  

    ### lambda ###
    tmp <- match("lambda", ininit, nomatch=NA)
    if (is.na(tmp)){
      if (prior.gspline$Ldistrib == "fixed") stop("init.gspline$lambda must be given")
      else{
        if (prior.gspline$Ldistrib == "gamma"){
          if (prior.gspline$Lequal) init.gspline$lambda <- rep(rgamma(1, shape=prior.gspline$Lshape[1], rate=prior.gspline$LinvScale[1]), nRandom) + 0.001
          else                      init.gspline$lambda <- rgamma(nRandom, shape=prior.gspline$Lshape, rate=prior.gspline$LinvScale) + 0.001
        }
        else{
          if (prior.gspline$Ldistrib == "sduniform"){
            if (prior.gspline$Lequal) init.gspline$lambda <- rep(1/(runif(1, min=0, max=prior.gspline$Lupper[1]))^2, nRandom)
            else                      init.gspline$lambda <- 1/(runif(nRandom, min=0, max=prior.gspline$Lupper))^2
          }
          else{
            stop("Incorrect prior.gspline$Ldistrib specified")
          }  
        }  
      }
    }
    if (length(init.gspline$lambda) == 1) init.gspline$lambda <- rep(init.gspline$lambda, nRandom)
    if (prior.gspline$Lequal) init.gspline$lambda <- rep(init.gspline$lambda[1], nRandom)        
    if (length(init.gspline$lambda) != nRandom) stop(paste("init.gspline$lambda must be a vector of length ", nRandom, sep=""))
    if (any(init.gspline$lambda <= 0)) stop("All components of init.gspline$lambda must be positive")
    Lambda <- init.gspline$lambda
    names(Lambda) <- paste("lambda", 1:nRandom, sep="")

    
    ### weights and a coefficients ###
    aaa <- list()
    www <- list()
    
    tmp <- match("weights", ininit, nomatch=NA)
    if (is.na(tmp)){                        
      for (i in 1:nRandom){   ## Give more or less normal distribution as the initial G-spline in each margin
                              ## by minimizing the 3rd order penalty
        Sdspline <- prior.gspline$sigma[i]
        if (Sdspline >= 0.95) Sdspline <- 0.95
        Knots <- ((-prior.gspline$K[i]):prior.gspline$K[i])*prior.gspline$delta[i]
        minp <- minPenalty(knots=Knots, sdspline=Sdspline, difforder=3, info=FALSE)            ### from package smoothSurv
        if (minp$fail) stop(paste("Unable to guess initial 'a' coefficients for margin ", i, ", give your own", sep=""))
        aaa[[i]] <- minp$spline[, "a coef."]
        aaa[[i]] <- aaa[[i]] - t(prior.gspline$Acontrast[[i]]) %*% aaa[[i]]
        names(aaa[[i]]) <- paste("a", i, ".", (-prior.gspline$K[i]):prior.gspline$K[i], sep="")
        exp.a <- exp(aaa[[i]])
        www[[i]] <- exp.a/sum(exp.a)
        names(www[[i]]) <- paste("w", i, ".", (-prior.gspline$K[i]):prior.gspline$K[i], sep="")          
      }
      if (nRandom == 1){
        init.gspline$acoef <- aaa[[1]]
        init.gspline$weights <- www[[1]]
      }
      else{
        if (nRandom == 2){       ## For bivariate models, take a product of marginal weights as joint weights (uncorrelated margins)
          if (prior.gspline$Aident != "reference") stop("prior.gspline$Aident must be reference")      
          init.gspline$weights <- www[[1]] %o% www[[2]]
          iref <- prior.gspline$Areference + prior.gspline$K + 1
          wref <- init.gspline$weights[iref[1], iref[2]]
          init.gspline$acoef <- log(init.gspline$weights/wref)
          rownames(init.gspline$acoef) <- paste("a", 1, ".", (-prior.gspline$K[1]):prior.gspline$K[1], sep="")
          colnames(init.gspline$acoef) <- paste("a", 2, ".", (-prior.gspline$K[2]):prior.gspline$K[2], sep="")

          ### TEMPORAR
          ##expacoef <- exp(init.gspline$acoef)
          ##print(sumexpa1 <- apply(expacoef, 1, sum))
          ##print(sumexpa2 <- apply(expacoef, 2, sum))
        }
        else{
          stop("only implemented for uni- and bivariate random effects")
        }  
      }  
    }
    else{
      if (nRandom == 1){
        if (!is.numeric(init.gspline$weights)) stop("init.gspline$weights must be a numeric vector")
        if (length(init.gspline$weights) != 2*prior.gspline$K[1] + 1) stop("init.gspline$weights has incorrect length")
        init.gspline$weights <- init.gspline$weights/sum(init.gspline$weights)
        wreftemp <- max(init.gspline$weights)
        init.gspline$acoef <- log(init.gspline$weights/wreftemp)
        init.gspline$acoef <- init.gspline$acoef - t(prior.gspline$Acontrast[[1]]) %*% init.gspline$acoef
        names(init.gspline$acoef) <- paste("a", 1, ".", (-prior.gspline$K[1]):prior.gspline$K[1], sep="")        
      }
      else{
        if (nRandom == 2){
          if (prior.gspline$Aident != "reference") stop("prior.gspline$Aident must be reference")          
          if (!is.matrix(init.gspline$weights)) stop("init.gspline$weights must be a matrix")
          if (nrow(init.gspline$weights) != 2*prior.gspline$K[1] + 1) stop("init.gspline$weights has incorrect number of rows")
          if (ncol(init.gspline$weights) != 2*prior.gspline$K[2] + 1) stop("init.gspline$weights has incorrect number of columns")
          iref <- prior.gspline$Areference + prior.gspline$K + 1
          wref <- init.gspline$weights[iref[1], iref[2]]
          init.gspline$acoef <- log(init.gspline$weights/wref)
          rownames(init.gspline$acoef) <- paste("a", 1, ".", (-prior.gspline$K[1]):prior.gspline$K[1], sep="")
          colnames(init.gspline$acoef) <- paste("a", 2, ".", (-prior.gspline$K[2]):prior.gspline$K[2], sep="")          
        }
        else{
          stop("only implemented for uni- and bivariate random effects")
        }  
      }  
    }
    lambda.a <- c(Lambda, init.gspline$acoef)
    
    
    ### allocations ###
    tmp <- match("alloc", ininit, nomatch=NA)
    if (is.na(tmp)){
      if (nRandom == 1) Std.Dev <- sqrt(init.random$var)
      else{
        if (length(init.random$var) != nRandom) stop(paste("init.random$var should be a vector of length ", nRandom, " at this stage", sep=""))
        Std.Dev <- sqrt(init.random$var)
      }
      init.gspline$alloc <- maxPosterProb(data=init.random$b, intercept=init.random$mean, std.dev=Std.Dev,
                                          K=prior.gspline$K, delta=prior.gspline$delta, sigma=prior.gspline$sigma)
    }
    else{
      if (nRandom == 1){
        nb <- length(init.random$b)
        if (length(init.gspline$alloc) != nb) stop("init.gspline$alloc has incorrect length")
        if (any(abs(init.gspline$alloc) > prior.gspline$K[1])) stop("init.gspline$alloc out of range")
      }
      else{
        nb <- nrow(init.random$b)
        if (nrow(init.random$alloc) != nb) stop("init.gspline$alloc has incorrect number of rows")
        if (ncol(init.random$alloc) != nRandom) stop("init.gspline$alloc has incorrect number of columns")
        Ktemp <- matrix(rep(prior.gspline$K, nb), ncol=nRandom, nrow=nb, byrow=TRUE)
        if (any(abs(init.gspline$alloc) > Ktemp)) stop("init.gspline$alloc out of range")        
      }        
    }  
  }  

  attr(init.gspline, "lambda.a") <- lambda.a
  return(init.gspline)  
}  


#* ********************************************************************************* */
headers.glmmAK <- function(dir, store, DES, prob.init, init.fixedML, init.MrandomML, hierar.center,
                            drandom=c("normal", "gspline"), model=c("cumlogit", "logpoisson"), prior.gspline)
{
  if (drandom != "none") drandom <- match.arg(drandom)
  model <- match.arg(model)
  if (missing(store)) instore <- "arnost"
  else                instore <- names(store)

  if (model == "cumlogit"){
    nFixed <- DES$nx + DES$C*DES$nv
    nRandom <- DES$nxb + DES$C*DES$nvb
  }else{
    if (model == "logpoisson"){
      nFixed <- DES$nx
      nRandom <- DES$nxb
    }else{
      stop("Unknown model")
    }  
  }  
  name.cluster <- names(table(DES$cluster))

  files.in.dir <- dir(dir)
  
  ### Index of the iteration (iteration.sim) ###
  sink(paste(dir, "/iteration.sim", sep = ""), append = FALSE)
  cat("iteration", "\n")
  sink()

  ### Regression parameters: fixed effects (betaF.sim) ###
  if (nFixed){
    sink(paste(dir, "/betaF.sim", sep = ""), append = FALSE)
    cat(names(init.fixedML), "\n", sep="   ")
    sink()
  }
  else{
    if ("betaF.sim" %in% files.in.dir) file.remove(paste(dir, "/betaF.sim", sep = ""))
  }  
  
  ### Log-likelihood (loglik.sim) ###
  sink(paste(dir, "/loglik.sim", sep = ""), append = FALSE)
  cat("value", "\n")
  sink()

  ### Acceptance indicators (accept.sim) ###
  #sink(paste(dir, "/accept.sim", sep = ""), append = FALSE)
  #if (nFixed & nRandom) cat("betaF", paste("b:", name.cluster, sep=""), "\n", sep="  ")
  #else{
  #  if (nFixed) cat("betaF", "\n", sep="  ")
  #  if (nRandom) cat(paste("b:", name.cluster, sep=""), "\n", sep="  ")
  #}  
  #sink()

  if (model == "cumlogit"){
    ### Category probabilities (probability.sim) ###
    tmp <- match("prob", instore, nomatch=NA)
    if (is.na(tmp)){
      if ("probability.sim" %in% files.in.dir) file.remove(paste(dir, "/probability.sim", sep = ""))
    }
    else{
      if (store$prob){
        sink(paste(dir, "/probability.sim", sep = ""), append = FALSE)
        cat(paste(colnames(prob.init), ":", rep(1:DES$ny, each=DES$C+1), sep=""), "\n", sep="   ")
        sink()
      }else{
        if ("probability.sim" %in% files.in.dir) file.remove(paste(dir, "/probability.sim", sep = ""))
      }  
    }      
  }else{
    if (model == "logpoisson"){
      ### Expected counts (probability.sim) ###
      tmp <- match("ecount", instore, nomatch=NA)
      if (is.na(tmp)){
        if ("expectcount.sim" %in% files.in.dir) file.remove(paste(dir, "/expectcount.sim", sep = ""))
      }
      else{
        if (store$ecount){
          sink(paste(dir, "/expectcount.sim", sep = ""), append = FALSE)
          cat(paste("ecount", 1:DES$ny, sep=""), "\n", sep="   ")
          sink()
        }else{
          if ("expectcount.sim" %in% files.in.dir) file.remove(paste(dir, "/expectcount.sim", sep = ""))
        }  
      }        
    }else{
      stop("Unknown model")
    }  
  }  

  if (nRandom){

    ### Means of random effects (betaR.sim). Only if there is hierarchical centering. ###
    if (hierar.center){    
      sink(paste(dir, "/betaR.sim", sep = ""), append = FALSE)
      cat(names(init.MrandomML), "\n", sep="   ")
      sink()
    }
    else{
      if ("betaR.sim" %in% files.in.dir) file.remove(paste(dir, "/betaR.sim", sep = ""))
    }  

    ### Covariance matrix of random effects and its inverse (varR.sim) ###
    ### -> differently for normal and G-spline random effects          ###
    sink(paste(dir, "/varR.sim", sep = ""), append = FALSE)    
    if (drandom == "normal"){
      tmp <- diag(nRandom)    
      rindex <- row(tmp)[lower.tri(row(tmp), diag=TRUE)]
      cindex <- col(tmp)[lower.tri(col(tmp), diag=TRUE)]
      cat(paste("varR.", rindex, ".", cindex, sep=""), paste("ivarR.", rindex, ".", cindex, sep=""), "\n", sep="   ")
    }
    if (drandom == "gspline"){
      rindex <- 1:nRandom
      cat(paste("varR.", rindex, sep=""), paste("ivarR.", rindex, sep=""), "\n", sep="   ")
    }
    sink()            
    
    ### Values of random effects ###
    tmp <- match("b", instore, nomatch=NA)
    if (is.na(tmp)){
#     if ("b.sim" %in% files.in.dir) file.remove(paste(dir, "/b.sim", sep = ""))    ### We will store b's anyway but not at each iteration
      sink(paste(dir, "/b.sim", sep = ""), append = FALSE)
      cat(paste(rep(names(init.MrandomML), DES$N), ":", rep(name.cluster, each=nRandom), sep=""), "\n", sep="   ")            
      sink()            
    }
    else{
      sink(paste(dir, "/b.sim", sep = ""), append = FALSE)
      cat(paste(rep(names(init.MrandomML), DES$N), ":", rep(name.cluster, each=nRandom), sep=""), "\n", sep="   ")      
      sink()
    }  
  }
  else{
    if ("betaR.sim" %in% files.in.dir) file.remove(paste(dir, "/betaR.sim", sep = ""))
    if ("varR.sim" %in% files.in.dir)  file.remove(paste(dir, "/varR.sim", sep = ""))
    if ("b.sim" %in% files.in.dir)     file.remove(paste(dir, "/b.sim", sep = ""))    
  }  

  if (drandom == "gspline"){
    
    ### Basis information on the G-spline (gspline.sim) ###
    if (missing(prior.gspline)) stop("prior.gspline must be given")
    sink(paste(dir, "/gspline.sim", sep = ""), append = FALSE)
    Dim <- attr(prior.gspline, "Dim")
    if (is.null(Dim)) stop("prior.gspline should have an attribute Dim at this stage")
    dPar <- attr(prior.gspline, "dPar")
    if (is.null(dPar)) stop("prior.gspline should have an attribute dPar at this stage")
    cat(names(Dim), names(dPar), "\n", sep="  ")
    cat(Dim, dPar, "\n", sep="  ")    
    sink()

    ### Log-weights (logweight.sim) ###
    if (nRandom == 1){
      nlogweight <- paste("a", ".", (-prior.gspline$K[1]):prior.gspline$K[1], sep="")

      if ("weight.sim" %in% files.in.dir)  file.remove(paste(dir, "/weight.sim", sep = ""))
      if ("knotInd.sim" %in% files.in.dir) file.remove(paste(dir, "/knotInd.sim", sep = ""))              
    }
    else{
      if (nRandom == 2){
        ll1 <- 2*prior.gspline$K[1] + 1
        ll2 <- 2*prior.gspline$K[2] + 1                
        rind <- rep((-prior.gspline$K[1]):prior.gspline$K[1], ll2)
        cind <- rep((-prior.gspline$K[2]):prior.gspline$K[2], each=ll1)
        nlogweight <- paste("a", ".", rind, ".", cind, sep="")

        sink(paste(dir, "/weight.sim", sep = ""), append = FALSE)
        cat(paste("w.", 1:(ll1*ll2), sep=""), "\n", sep="   ")
        sink()

        sink(paste(dir, "/knotInd.sim", sep = ""), append = FALSE)        
        cat(c("ncomponents", paste("ind.", 1:(ll1*ll2), sep="")), "\n", sep="   ")
        sink()        
      }
      else{
        stop("nRandom must be <= 2 when distribution of random effects is G-spline")
      }  
    }  
    sink(paste(dir, "/logweight.sim", sep = ""), append = FALSE)
    cat(nlogweight, "\n", sep="   ")
    sink()

    ### Moments of the G-spline
    sink(paste(dir, "/gmoment.sim", sep = ""), append = FALSE)
    if (nRandom == 1){
      cat("gmean   gvar\n")
    }
    else{
      if (nRandom == 2){
        cat("gmean1   gmean2   gvar1   gcovar21   gvar2\n")
      }
      else{
        stop("nRandom must be <= 2 when distribution of random effects is G-spline")
      }  
    }  
    sink()    

    ### Adjusted betaR (effect of covariates that enter into random effects)
    sink(paste(dir, "/betaRadj.sim", sep = ""), append = FALSE)
    cat(names(init.MrandomML), "\n", sep="   ")
    sink()

    ### Adjusted varR (covariance matrix of the random effect taking into account G-spline variance)
    sink(paste(dir, "/varRadj.sim", sep = ""), append = FALSE)
    if (nRandom == 1){
      cat("varR.1\n")
    }else{
      if (nRandom == 2){
        cat("varR.1.1   varR.2.1   varR.2.2\n")
      }else{
        stop("nRandom must be <= 2 when distribution of random effects is G-spline")
      }  
    }      
    sink()
    
    ### Allocations (alloc.sim) ###
    tmp <- match("alloc", instore, nomatch=NA)
    if (is.na(tmp)){
#     if ("alloc.sim" %in% files.in.dir) file.remove(paste(dir, "/alloc.sim", sep = ""))    ### We will store alloc's anyway but not at each iteration
      if (nRandom == 1){
        sink(paste(dir, "/alloc.sim", sep = ""), append = FALSE)
        #cat(paste(rep(names(init.MrandomML), DES$N), ":", rep(name.cluster, each=nRandom), sep=""), "\n", sep="   ")
        cat(paste("b:", name.cluster, sep=""), "\n", sep="   ")
        sink()
      }else{
        if (nRandom == 2){
          sink(paste(dir, "/alloc.sim", sep = ""), append = FALSE)
          cat(paste("b:", name.cluster, sep=""), "\n", sep="   ")
          sink()
        }else{
          stop("nRandom must be <= 2 when distribution of random effects is G-spline")
        }  
      }  
    }
    else{
      if (nRandom == 1){
        sink(paste(dir, "/alloc.sim", sep = ""), append = FALSE)
        #cat(paste(rep(names(init.MrandomML), DES$N), ":", rep(name.cluster, each=nRandom), sep=""), "\n", sep="   ")
        cat(paste("b:", name.cluster, sep=""), "\n", sep="   ")
        sink()
      }else{
        if (nRandom == 2){
          sink(paste(dir, "/alloc.sim", sep = ""), append = FALSE)
          cat(paste("b:", name.cluster, sep=""), "\n", sep="   ")
          sink()
        }else{
          stop("nRandom must be <= 2 when distribution of random effects is G-spline")
        }  
      }  
    }

    ### MRF precisions (lambda.sim) ###
    sink(paste(dir, "/lambda.sim", sep = ""), append = FALSE)
    if (prior.gspline$Lequal) cat("lambda\n")
    else                      cat(paste("lambda.", 1:nRandom, sep=""), "\n", sep="   ")
    sink()    
  }
  else{
    if ("gspline.sim" %in% files.in.dir)   file.remove(paste(dir, "/gspline.sim", sep = ""))
    if ("logweight.sim" %in% files.in.dir) file.remove(paste(dir, "/logweight.sim", sep = ""))
    if ("alloc.sim" %in% files.in.dir)     file.remove(paste(dir, "/alloc.sim", sep = ""))
    if ("gmoment.sim" %in% files.in.dir)   file.remove(paste(dir, "/gmoment.sim", sep = ""))

    if ("weight.sim" %in% files.in.dir)  file.remove(paste(dir, "/weight.sim", sep = ""))
    if ("knotInd.sim" %in% files.in.dir) file.remove(paste(dir, "/knotInd.sim", sep = ""))

    if ("betaRadj.sim" %in% files.in.dir) file.remove(paste(dir, "/betaRadj.sim", sep = ""))
    if ("varRadj.sim" %in% files.in.dir)  file.remove(paste(dir, "/varRadj.sim", sep = ""))        
  }    
  
  return(invisible(dir))    
}  
  
