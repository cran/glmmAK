###################################################
### chunk number 1: options
###################################################
options(width=85)


###################################################
### chunk number 2: invoke-man-pages-in-html eval=FALSE
###################################################
help(logpoissonRE, package=glmmAK, htmlhelp=TRUE)
#help(logpoissonRE.predict, package=glmmAK, htmlhelp=TRUE)
help(summaryGspline2, package=glmmAK, htmlhelp=TRUE)


###################################################
### chunk number 3: start
###################################################
library(glmmAK)
root <- "/home/komarek/Rlib/glmmAK/Doc/"
setwd(root)
data(epileptic)
data(epilepticBC)


###################################################
### chunk number 4: data.summary
###################################################
summary(epileptic)
summary(epilepticBC)


###################################################
### chunk number 5: priorBeta-fixedEffects
###################################################
prior.fixed <- list(mean=0, var=10000)


###################################################
### chunk number 6: priorPGM
###################################################
prior.gspline <- list(K=15, delta=0.3, sigma=0.2, CARorder=3,
                      Ldistrib="gamma", Lequal=FALSE, Lshape=1, LinvScale=0.005,
                      AtypeUpdate="slice")


###################################################
### chunk number 7: priorPGM-alternative
###################################################
prior.gspline.Alternative <- list(K=c(15, 10), delta=c(0.3, 0.6), sigma=c(0.2, 0.4), CARorder=3,
                      Ldistrib="gamma", Lequal=FALSE, Lshape=c(1, 0.001), LinvScale=c(0.005, 0.001),
                      AtypeUpdate="slice")


###################################################
### chunk number 8: priorb-randomEffects-inPGM-GLMM-nhc
###################################################
prior.random.gspl.nhc <- list(Ddistrib="gamma", Dshape=1, DinvScale=0.005)


###################################################
### chunk number 9: priorb-randomEffects-inPGM-GLMM-hc
###################################################
prior.random.gspl.hc  <- list(Mdistrib="normal", Mmean=0, Mvar=10000, Ddistrib="gamma", Dshape=1, DinvScale=0.005)


###################################################
### chunk number 10: priorb-randomEffects-inPGM-GLMM-nhc-uniform-prior-on-tau
###################################################
prior.random.gspl.nhc.Unif <- list(Ddistrib="sduniform", Dupper=c(100, 200))


###################################################
### chunk number 11: priorb-randomEffects-inNormal-GLMM-nhc
###################################################
prior.random.norm.nhc  <- list(Ddistrib="wishart", Ddf=2, DinvScale=0.005)


###################################################
### chunk number 12: priorb-randomEffects-inNormal-GLMM-hc
###################################################
prior.random.norm.hc  <- list(Mdistrib="normal", Mmean=0, Mvar=10000, Ddistrib="wishart", Ddf=2, DinvScale=0.005)


###################################################
### chunk number 13: create.dirs.to.store.MCMC
###################################################
if (!("chEpileptic" %in% dir(root))) dir.create(paste(root, "chEpileptic", sep=""))
dirNames <- c("PGM_nhc", "PGM_hc", "Normal_nhc", "Normal_hc")
dirPaths <- paste(root, "chEpileptic/", dirNames, "/", sep="")
for (i in 1:length(dirPaths)){
  if(!(dirNames[i] %in% dir(paste(root, "chEpileptic", sep="")))) dir.create(dirPaths[i])
}
names(dirPaths) <- c("PGM_nhc", "PGM_hc", "Normal_nhc", "Normal_hc")


###################################################
### chunk number 14: print.dirs.to.store.MCMC
###################################################
print(dirPaths)


###################################################
### chunk number 15: create.matrix.of.covariates
###################################################
X2mat <- epilepticBC[,c("Base", "Trt", "Base.Trt", "Age")]
Xb2mat <- data.frame(Visit=epilepticBC[,"Visit"])


###################################################
### chunk number 16: print.first.few.rows.of.covariate.matrices
###################################################
print(X2mat[1:6,])
print(Xb2mat[1:6,])


###################################################
### chunk number 17: length.of.MCMC
###################################################
nsimul <- list(niter=2000, nthin=10, nburn=1000, nwrite=100)


###################################################
### chunk number 18: MCMC.for.PGM-GLMM-nhc eval=FALSE
###################################################
fit.PGM.nhc <- logpoissonRE(y=epilepticBC$Seizure, x=X2mat, xb=Xb2mat, cluster=epilepticBC$id,
                            intcpt.random=TRUE, hierar.center=FALSE, drandom="gspline",
                            prior.fixed=prior.fixed, prior.random=prior.random.gspl.nhc, prior.gspline=prior.gspline,
                            nsimul=nsimul, store=list(ecount=FALSE, b=TRUE), dir=dirPaths["PGM_nhc"])


###################################################
### chunk number 19: MCMC.for.PGM-GLMM-hc eval=FALSE
###################################################
fit.PGM.hc <- logpoissonRE(y=epilepticBC$Seizure, x=X2mat, xb=Xb2mat, cluster=epilepticBC$id,
                           intcpt.random=TRUE, hierar.center=TRUE, drandom="gspline",
                           prior.fixed=prior.fixed, prior.random=prior.random.gspl.hc, prior.gspline=prior.gspline,
                           nsimul=nsimul, store=list(ecount=FALSE, b=TRUE), dir=dirPaths["PGM_hc"])


###################################################
### chunk number 20: MCMC.for.Normal-GLMM-nhc eval=FALSE
###################################################
fit.Normal.nhc <- logpoissonRE(y=epilepticBC$Seizure, x=X2mat, xb=Xb2mat, cluster=epilepticBC$id,
                               intcpt.random=TRUE, hierar.center=FALSE, drandom="normal",
                               prior.fixed=prior.fixed, prior.random=prior.random.norm.nhc,
                               nsimul=nsimul, store=list(ecount=FALSE, b=TRUE), dir=dirPaths["Normal_nhc"])


###################################################
### chunk number 21: MCMC.for.Normal-GLMM-hc eval=FALSE
###################################################
fit.Normal.hc <- logpoissonRE(y=epilepticBC$Seizure, x=X2mat, xb=Xb2mat, cluster=epilepticBC$id,
                              intcpt.random=TRUE, hierar.center=TRUE, drandom="normal",
                              prior.fixed=prior.fixed, prior.random=prior.random.norm.hc,
                              nsimul=nsimul, store=list(ecount=FALSE, b=TRUE), dir=dirPaths["Normal_hc"])


###################################################
### chunk number 22: length.of.MCMC.in.the.paper
###################################################
nsimul <- list(niter=50000, nthin=130, nburn=25000, nwrite=1000)


###################################################
### chunk number 23: read.chains eval=FALSE
###################################################
chPGM.nhc <- glmmAK.files2coda(dir=dirPaths["PGM_nhc"], drandom="gspline", skip=0)
chPGM.hc  <- glmmAK.files2coda(dir=dirPaths["PGM_hc"], drandom="gspline", skip=0)
chNormal.nhc <- glmmAK.files2coda(dir=dirPaths["Normal_nhc"], drandom="normal", skip=0)
chNormal.hc  <- glmmAK.files2coda(dir=dirPaths["Normal_hc"], drandom="normal", skip=0)


###################################################
### chunk number 24: read.chains.PGM-GLMM-nhc eval=FALSE
###################################################
iters    <- scanFH(paste(dirPaths["PGM_nhc"], "iteration.sim", sep=""))
betaF    <- scanFH(paste(dirPaths["PGM_nhc"], "betaF.sim", sep=""))
betaRadj <- scanFH(paste(dirPaths["PGM_nhc"], "betaRadj.sim", sep=""))
varRadj  <- scanFH(paste(dirPaths["PGM_nhc"], "varRadj.sim", sep=""))

chPGM.nhc <- mcmc(data.frame(Base=betaF[,"Base"], Trt=betaF[,"Trt"], Base.Trt=betaF[,"Base.Trt"], Age=betaF[,"Age"],
                             Intcpt=betaRadj[,"(Intercept)"], Visit=betaRadj[,"Visit"],
                             SDIntcpt=sqrt(varRadj[,"varR.1.1"]), SDVisit=sqrt(varRadj[,"varR.2.2"]), 
                             Corr=varRadj[,"varR.2.1"]/sqrt(varRadj[,"varR.1.1"]*varRadj[,"varR.2.2"])),
                  start=iters[1,1])
rm(list=c("iters", "betaF", "betaRadj", "varRadj"))


###################################################
### chunk number 25: read.chains.PGM-GLMM-hc eval=FALSE
###################################################
iters    <- scanFH(paste(dirPaths["PGM_hc"], "iteration.sim", sep=""))
betaF    <- scanFH(paste(dirPaths["PGM_hc"], "betaF.sim", sep=""))
betaRadj <- scanFH(paste(dirPaths["PGM_hc"], "betaRadj.sim", sep=""))
varRadj  <- scanFH(paste(dirPaths["PGM_hc"], "varRadj.sim", sep=""))

chPGM.hc <- mcmc(data.frame(Base=betaF[,"Base"], Trt=betaF[,"Trt"], Base.Trt=betaF[,"Base.Trt"], Age=betaF[,"Age"],
                            Intcpt=betaRadj[,"(Intercept)"], Visit=betaRadj[,"Visit"],
                            SDIntcpt=sqrt(varRadj[,"varR.1.1"]), SDVisit=sqrt(varRadj[,"varR.2.2"]), 
                            Corr=varRadj[,"varR.2.1"]/sqrt(varRadj[,"varR.1.1"]*varRadj[,"varR.2.2"])),
                 start=iters[1,1])
rm(list=c("iters", "betaF", "betaRadj", "varRadj"))


###################################################
### chunk number 26: read.chains.Normal-GLMM-nhc eval=FALSE
###################################################
iters <- scanFH(paste(dirPaths["Normal_nhc"], "iteration.sim", sep=""))
betaF <- scanFH(paste(dirPaths["Normal_nhc"], "betaF.sim", sep=""))
varR  <- scanFH(paste(dirPaths["Normal_nhc"], "varR.sim", sep=""))

chNormal.nhc <- mcmc(data.frame(Base=betaF[,"Base"], Trt=betaF[,"Trt"], Base.Trt=betaF[,"Base.Trt"], Age=betaF[,"Age"],
                                Intcpt=betaF[,"(Intercept)"], Visit=betaF[,"Visit"],
                                SDIntcpt=sqrt(varR[,"varR.1.1"]), SDVisit=sqrt(varR[,"varR.2.2"]), 
                                Corr=varR[,"varR.2.1"]/sqrt(varR[,"varR.1.1"]*varR[,"varR.2.2"])),
                  start=iters[1,1])
rm(list=c("iters", "betaF", "varR"))


###################################################
### chunk number 27: read.chains.Normal-GLMM-hc eval=FALSE
###################################################
iters <- scanFH(paste(dirPaths["Normal_hc"], "iteration.sim", sep=""))
betaF <- scanFH(paste(dirPaths["Normal_hc"], "betaF.sim", sep=""))
betaR <- scanFH(paste(dirPaths["Normal_hc"], "betaR.sim", sep=""))
varR  <- scanFH(paste(dirPaths["Normal_hc"], "varR.sim", sep=""))

chNormal.hc <- mcmc(data.frame(Base=betaF[,"Base"], Trt=betaF[,"Trt"], Base.Trt=betaF[,"Base.Trt"], Age=betaF[,"Age"],
                               Intcpt=betaR[,"(Intercept)"], Visit=betaR[,"Visit"],
                               SDIntcpt=sqrt(varR[,"varR.1.1"]), SDVisit=sqrt(varR[,"varR.2.2"]), 
                               Corr=varR[,"varR.2.1"]/sqrt(varR[,"varR.1.1"]*varR[,"varR.2.2"])),
                  start=iters[1,1])
rm(list=c("iters", "betaF", "betaR", "varR"))


###################################################
### chunk number 28: summary.coda eval=FALSE
###################################################
summary(chPGM.nhc)
summary(chPGM.hc)
summary(chNormal.nhc)
summary(chNormal.hc)


###################################################
### chunk number 29: Bayesian.P-values eval=FALSE
###################################################
params <- c("Base", "Trt", "Base.Trt", "Age", "Visit")
BPvalue(chPGM.nhc[,params])
BPvalue(chPGM.hc[,params])
BPvalue(chNormal.nhc[,params])
BPvalue(chNormal.hc[,params])


###################################################
### chunk number 30: HPD.intervals eval=FALSE
###################################################
HPDinterval(chPGM.nhc, prob=0.95)
HPDinterval(chPGM.hc, prob=0.95)
HPDinterval(chNormal.nhc, prob=0.95)
HPDinterval(chNormal.hc, prob=0.95)


###################################################
### chunk number 31: standardized.random.effect.density eval=FALSE
###################################################
knots1 <- seq(-4.5, 4.5, by=0.3)
knots2 <- seq(-4.5, 4.5, by=0.3)
sigma1 <- 0.2
sigma2 <- 0.2
grid1 <- seq(-3, 3, length=20)
grid2 <- seq(-3, 3, length=20) 

### PGM GLMM(nhc)
stPGM.nhc <- summaryGspline2(x1=grid1, x2=grid2, mu1=knots1, mu2=knots2, sigma1=sigma1, sigma2=sigma2, standard=TRUE, 
                             probs=c(0.025, 0.25, 0.5, 0.75, 0.975), values=FALSE, dir=dirPaths["PGM_nhc"])

### PGM GLMM(hc)
stPGM.hc <- summaryGspline2(x1=grid1, x2=grid2, mu1=knots1, mu2=knots2, sigma1=sigma1, sigma2=sigma2, standard=TRUE, 
                            probs=c(0.025, 0.25, 0.5, 0.75, 0.975), values=FALSE, dir=dirPaths["PGM_hc"])


###################################################
### chunk number 32: figure-stdensity-PGMnhc eval=FALSE
###################################################
obj  <- stPGM.nhc$summary
obj1 <- stPGM.nhc$summary1
obj2 <- stPGM.nhc$summary2

par(mfrow=c(2, 2), bty="n", mar=c(4, 4, 1, 0)+0.1)

### Joint density (posterior mean only)
contour(obj$x1, obj$x2, obj$Mean, col="red", xlab="b1[st]", ylab="b2[st]")
persp(obj$x1, obj$x2, obj$Mean, col="seagreen3", theta=30, phi=60, xlab="b1[st]", ylab="b2[st]", zlab="g(b1[st],b2[st])")

### Marginal random intercept density (posterior mean, 2.5% and 97.5% quantiles)
plot(obj1$x, obj1[,"97.5%"], type="l", lty=1, col="red", xlab="b1[st]", ylab="g(b1[st])", main="Random intercept")
lines(obj1$x, obj1[,"2.5%"], lty=2, col="red")
lines(obj1$x, obj1$Mean, lty=1, col="blue")

### Marginal random Visit effect density (posterior mean, 2.5% and 97.5% quantiles)
plot(obj2$x, obj2[,"97.5%"], type="l", lty=1, col="red", xlab="b2[st]", ylab="g(b2[st])", main="Random Visit effect")
lines(obj2$x, obj2[,"2.5%"], lty=2, col="red")
lines(obj2$x, obj2$Mean, lty=1, col="blue")


###################################################
### chunk number 33: extract-id-from-data
###################################################
IDNR <- unique(epilepticBC$id)
IDNR0 <- unique(subset(epilepticBC, Trt==0)$id)
IDNR1 <- unique(subset(epilepticBC, Trt==1)$id)
index.tr0 <- (1:length(IDNR))[IDNR %in% IDNR0]
index.tr1 <- (1:length(IDNR))[IDNR %in% IDNR1]


###################################################
### chunk number 34: read-b-in-PGM-GLMM-nhc eval=FALSE
###################################################
betab.PGMnhc <- scanFH(paste(dirPaths["PGM_nhc"], "/betaF.sim", sep=""))[,c("(Intercept)", "Visit")]
b.PGMnhc <- scanFH(paste(dirPaths["PGM_nhc"], "b.sim", sep = "")) + as.matrix(betab.PGMnhc)
colnames(b.PGMnhc) <- paste(c("Intcpt", "Visit"), rep(IDNR, each=2), sep="")


###################################################
### chunk number 35: mean-median-bIntcpt-in-PGM-GLMM-nhc eval=FALSE
###################################################
indIntcpt <- seq(1, ncol(b.PGMnhc)-1, by=2)
bIntcptMean.PGMnhc <- apply(b.PGMnhc[,indIntcpt], 2, mean)
bIntcptMedian.PGMnhc <- apply(b.PGMnhc[,indIntcpt], 2, median)


###################################################
### chunk number 36: mean-median-bVisit-in-PGM-GLMM-nhc eval=FALSE
###################################################
indVisit <- seq(2, ncol(b.PGMnhc), by=2)
bVisitMean.PGMnhc <- apply(b.PGMnhc[,indVisit], 2, mean)
bVisitMedian.PGMnhc <- apply(b.PGMnhc[,indVisit], 2, median)


###################################################
### chunk number 37: figure-posterior-mean-and-median-of-individual-random-effects-in-PGM-GLMM-nhc eval=FALSE
###################################################
showid <- c(112, 135, 225, 227, 232)
index.show <- IDNR %in% showid

par(mfrow=c(2, 1), bty="n", mar=c(4, 4, 4, 1)+0.1)

### Posterior means
plot(bIntcptMean.PGMnhc[index.tr0], bVisitMean.PGMnhc[index.tr0], pch=1, col="red", xlab="beta1+b1", ylab="beta2+b2", 
     xlim=range(bIntcptMean.PGMnhc), ylim=range(bVisitMean.PGMnhc), main="Posterior means")
points(bIntcptMean.PGMnhc[index.tr1], bVisitMean.PGMnhc[index.tr1], pch=7, col="darkgreen")
text(bIntcptMean.PGMnhc[index.show]+0.005, bVisitMean.PGMnhc[index.show], labels=IDNR[index.show], pos=4)

### Posterior medians
plot(bIntcptMedian.PGMnhc[index.tr0], bVisitMedian.PGMnhc[index.tr0], pch=1, col="red", xlab="beta1+b1", ylab="beta2+b2", 
     xlim=range(bIntcptMedian.PGMnhc), ylim=range(bVisitMedian.PGMnhc), main="Posterior medians")
points(bIntcptMedian.PGMnhc[index.tr1], bVisitMedian.PGMnhc[index.tr1], pch=7, col="darkgreen")
text(bIntcptMedian.PGMnhc[index.show]+0.005, bVisitMedian.PGMnhc[index.show], labels=IDNR[index.show], pos=4)


