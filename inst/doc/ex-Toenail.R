###################################################
### chunk number 1: options
###################################################
options(width=85)


###################################################
### chunk number 2: invoke-man-pages-in-html eval=FALSE
###################################################
help(cumlogitRE, package=glmmAK, htmlhelp=TRUE)
#help(cumlogitRE.predict, package=glmmAK, htmlhelp=TRUE)
help(summaryGspline1, package=glmmAK, htmlhelp=TRUE)


###################################################
### chunk number 3: start
###################################################
library(glmmAK)
root <- "/home/komarek/Rlib/glmmAK/Doc/"
setwd(root)
data(toenail)


###################################################
### chunk number 4: data.summary
###################################################
summary(toenail)


###################################################
### chunk number 5: priorBeta-fixedEffects
###################################################
prior.fixed <- list(mean=0, var=10000)


###################################################
### chunk number 6: priorPGM-slice-sampler
###################################################
prior.gspline.Slice <- list(K=15, delta=0.3, sigma=0.2, CARorder=3,
                            Ldistrib="gamma", Lequal=FALSE, Lshape=1, LinvScale=0.005,
                            AtypeUpdate="slice")


###################################################
### chunk number 7: priorPGM-alternative-block-update
###################################################
prior.gspline.Block <- list(K=15, delta=0.3, sigma=0.2, CARorder=3,
                            Ldistrib="gamma", Lequal=FALSE, Lshape=1, LinvScale=0.005,
                            AtypeUpdate="block")


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
prior.random.gspl.nhc.Unif <- list(Ddistrib="sduniform", Dupper=100)


###################################################
### chunk number 11: priorb-randomEffects-inNormal-GLMM-nhc
###################################################
prior.random.norm.nhc  <- list(Ddistrib="gamma", Dshape=1, DinvScale=0.005)


###################################################
### chunk number 12: priorb-randomEffects-inNormal-GLMM-hc
###################################################
prior.random.norm.hc  <- list(Mdistrib="normal", Mmean=0, Mvar=10000, Ddistrib="gamma", Dshape=1, DinvScale=0.005)


###################################################
### chunk number 13: priorb-randomEffects-inPGM-GLMM-nhc-uniform-prior-on-tau
###################################################
prior.random.norm.nhc.Unif <- list(Ddistrib="sduniform", Dupper=100)


###################################################
### chunk number 14: create.dirs.to.store.MCMC
###################################################
if (!("chToenail" %in% dir(root))) dir.create(paste(root, "chToenail", sep=""))
dirNames <- c("PGM_nhc", "PGM_hc", "Normal_nhc", "Normal_hc")
dirPaths <- paste(root, "chToenail/", dirNames, "/", sep="")
for (i in 1:length(dirPaths)){
  if(!(dirNames[i] %in% dir(paste(root, "chToenail", sep="")))) dir.create(dirPaths[i])
}
names(dirPaths) <- c("PGM_nhc", "PGM_hc", "Normal_nhc", "Normal_hc")


###################################################
### chunk number 15: print.dirs.to.store.MCMC
###################################################
print(dirPaths)


###################################################
### chunk number 16: create.matrix.of.covariates
###################################################
iXmat <- data.frame(trt=toenail$trt, time=toenail$time, trt.time=toenail$trt*toenail$time)


###################################################
### chunk number 17: length.of.MCMC
###################################################
nsimul <- list(niter=4000, nthin=100, nburn=2000, nwrite=100)


###################################################
### chunk number 18: MCMC.for.PGM-GLMM-nhc eval=FALSE
###################################################
fit.PGM.nhc <- cumlogitRE(y=toenail$infect, x=iXmat, cluster=toenail$idnr,
                     intcpt.random=TRUE, hierar.center=FALSE, C=1, drandom="gspline",
                     prior.fixed=prior.fixed, prior.random=prior.random.gspl.nhc, 
                     prior.gspline=prior.gspline.Slice,
                     nsimul=nsimul, store=list(prob=FALSE, b=TRUE), dir=dirPaths["PGM_nhc"])


###################################################
### chunk number 19: MCMC.for.PGM-GLMM-hc eval=FALSE
###################################################
fit.PGM.hc <- cumlogitRE(y=toenail$infect, x=iXmat, cluster=toenail$idnr,
                     intcpt.random=TRUE, hierar.center=TRUE, C=1, drandom="gspline",
                     prior.fixed=prior.fixed, prior.random=prior.random.gspl.hc, 
                     prior.gspline=prior.gspline.Slice,
                     nsimul=nsimul, store=list(prob=FALSE, b=TRUE), dir=dirPaths["PGM_hc"])


###################################################
### chunk number 20: MCMC.for.Normal-GLMM-nhc eval=FALSE
###################################################
fit.Normal.nhc <- cumlogitRE(y=toenail$infect, x=iXmat, cluster=toenail$idnr,
                     intcpt.random=TRUE, hierar.center=FALSE, C=1, drandom="normal",
                     prior.fixed=prior.fixed, prior.random=prior.random.norm.nhc, 
                     nsimul=nsimul, store=list(prob=FALSE, b=TRUE), dir=dirPaths["Normal_nhc"])


###################################################
### chunk number 21: MCMC.for.Normal-GLMM-hc eval=FALSE
###################################################
fit.Normal.hc <- cumlogitRE(y=toenail$infect, x=iXmat, cluster=toenail$idnr,
                     intcpt.random=TRUE, hierar.center=TRUE, C=1, drandom="normal",
                     prior.fixed=prior.fixed, prior.random=prior.random.norm.hc, 
                     nsimul=nsimul, store=list(prob=FALSE, b=TRUE), dir=dirPaths["Normal_hc"])


###################################################
### chunk number 22: read.chains eval=FALSE
###################################################
chPGM.nhc <- glmmAK.files2coda(dir=dirPaths["PGM_nhc"], drandom="gspline", skip=0)
chPGM.hc  <- glmmAK.files2coda(dir=dirPaths["PGM_hc"], drandom="gspline", skip=0)
chNormal.nhc <- glmmAK.files2coda(dir=dirPaths["Normal_nhc"], drandom="normal", skip=0)
chNormal.hc  <- glmmAK.files2coda(dir=dirPaths["Normal_hc"], drandom="normal", skip=0)


###################################################
### chunk number 23: read.chains.PGM-GLMM-nhc eval=FALSE
###################################################
iters    <- scanFH(paste(dirPaths["PGM_nhc"], "iteration.sim", sep=""))
betaF    <- scanFH(paste(dirPaths["PGM_nhc"], "betaF.sim", sep=""))
betaRadj <- scanFH(paste(dirPaths["PGM_nhc"], "betaRadj.sim", sep=""))
varRadj  <- scanFH(paste(dirPaths["PGM_nhc"], "varRadj.sim", sep=""))

chPGM.nhc <- mcmc(data.frame(Trt=betaF[,"trt"], Time=betaF[,"time"], Trt.Time=betaF[,"trt.time"],
                             Meanb=betaRadj[,"(Intercept)"], SDb=sqrt(varRadj[,"varR.1"])),
                  start=iters[1,1])
rm(list=c("iters", "betaF", "betaRadj", "varRadj"))


###################################################
### chunk number 24: read.chains.PGM-GLMM-hc eval=FALSE
###################################################
iters    <- scanFH(paste(dirPaths["PGM_hc"], "iteration.sim", sep=""))
betaF    <- scanFH(paste(dirPaths["PGM_hc"], "betaF.sim", sep=""))
betaRadj <- scanFH(paste(dirPaths["PGM_hc"], "betaRadj.sim", sep=""))
varRadj  <- scanFH(paste(dirPaths["PGM_hc"], "varRadj.sim", sep=""))

chPGM.hc <- mcmc(data.frame(Trt=betaF[,"trt"], Time=betaF[,"time"], Trt.Time=betaF[,"trt.time"],
                            Meanb=betaRadj[,"(Intercept)"], SDb=sqrt(varRadj[,"varR.1"])),
                  start=iters[1,1])
rm(list=c("iters", "betaF", "betaRadj", "varRadj"))


###################################################
### chunk number 25: read.chains.Normal-GLMM-nhc eval=FALSE
###################################################
iters <- scanFH(paste(dirPaths["Normal_nhc"], "iteration.sim", sep=""))
betaF <- scanFH(paste(dirPaths["Normal_nhc"], "betaF.sim", sep=""))
varR  <- scanFH(paste(dirPaths["Normal_nhc"], "varR.sim", sep=""))

chNormal.nhc <- mcmc(data.frame(Trt=betaF[,"trt"], Time=betaF[,"time"], Trt.Time=betaF[,"trt.time"],
                                Meanb=betaF[,"(Intercept)"], SDb=sqrt(varR[,"varR.1.1"])),
                     start=iters[1,1])
rm(list=c("iters", "betaF", "varR"))


###################################################
### chunk number 26: read.chains.Normal-GLMM-hc eval=FALSE
###################################################
iters <- scanFH(paste(dirPaths["Normal_hc"], "iteration.sim", sep=""))
betaF <- scanFH(paste(dirPaths["Normal_hc"], "betaF.sim", sep=""))
betaR <- scanFH(paste(dirPaths["Normal_hc"], "betaR.sim", sep=""))
varR  <- scanFH(paste(dirPaths["Normal_hc"], "varR.sim", sep=""))

chNormal.hc <- mcmc(data.frame(Trt=betaF[,"trt"], Time=betaF[,"time"], Trt.Time=betaF[,"trt.time"],
                               Meanb=betaR[,"(Intercept)"], SDb=sqrt(varR[,"varR.1.1"])),
                  start=iters[1,1])
rm(list=c("iters", "betaF", "betaR", "varR"))


###################################################
### chunk number 27: summary.coda eval=FALSE
###################################################
summary(chPGM.nhc)
summary(chPGM.hc)
summary(chNormal.nhc)
summary(chNormal.hc)


###################################################
### chunk number 28: Bayesian.P-values eval=FALSE
###################################################
BPvalue(chPGM.nhc)
BPvalue(chPGM.hc)
BPvalue(chNormal.nhc)
BPvalue(chNormal.hc)


###################################################
### chunk number 29: HPD.intervals eval=FALSE
###################################################
HPDinterval(chPGM.nhc, prob=0.95)
HPDinterval(chPGM.hc, prob=0.95)
HPDinterval(chNormal.nhc, prob=0.95)
HPDinterval(chNormal.hc, prob=0.95)


###################################################
### chunk number 30: standardized.random.intercept.density eval=FALSE
###################################################
knots <- seq(-4.5, 4.5, by=0.3)
sigma <- 0.2
grid <- seq(-2, 4.5, length=300)

### PGM GLMM(nhc)
stPGM.nhc <- summaryGspline1(x=grid, mu=knots, sigma=sigma, standard=TRUE, 
                             probs=c(0.025, 0.25, 0.5, 0.75, 0.975), values=TRUE, dir=dirPaths["PGM_nhc"])

### PGM GLMM(hc)
stPGM.hc <- summaryGspline1(x=grid, mu=knots, sigma=sigma, standard=TRUE, 
                            probs=c(0.025, 0.25, 0.5, 0.75, 0.975), values=TRUE, dir=dirPaths["PGM_hc"])


###################################################
### chunk number 31: figure-stdensity-PGMnhc eval=FALSE
###################################################
par(bty="n", mar=c(4, 4, 1, 1)+0.1)
layout(matrix(c(1,1,2,2, 3,4,5,6, 7,8,9,10, 11,12,13,14), ncol=4, byrow=TRUE))

### Posterior mean
plot(stPGM.nhc$summary$x, stPGM.nhc$summary$Mean, type="l", xlab="b[st]", ylab="g(b[st])", col="blue")

### Posterior median and 2.5%, 97.5% quantiles
plot(stPGM.nhc$summary$x, stPGM.nhc$summary[,"97.5%"], type="l", lty=2, xlab="b[st]", ylab="g(b[st])", col="red")
lines(stPGM.nhc$summary$x, stPGM.nhc$summary[,"2.5%"], lty=2, col="red")
lines(stPGM.nhc$summary$x, stPGM.nhc$summary[,"50%"], lty=1, col="blue")

### Sampled densities at selected iterations
ylim <- c(0, 2.8)
#for (iters in c(1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 950, 1000)){
for (iters in c(1, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 1900, 2000)){
  plot(stPGM.nhc$summary$x, stPGM.nhc$values[iters,], type="l", xlab="", ylab="", col="darkblue", ylim=ylim)
  title(main=paste("Iteration ", nsimul$nburn+iters, sep=""))
}


###################################################
### chunk number 32: unstandardized.random.intercept.density-knots-etc
###################################################
knots <- seq(-4.5, 4.5, by=0.3)
sigma <- 0.2
grid <- seq(-10, 13, length=300)


###################################################
### chunk number 33: unstandardized.random.intercept.density-PGM-nhc eval=FALSE
###################################################
shift.nhc <- scanFH(paste(dirPaths["PGM_nhc"], "betaF.sim", sep = ""))[, "(Intercept)"]
scale.nhc <- sqrt(scanFH(paste(dirPaths["PGM_nhc"], "varR.sim", sep = ""))[, "varR.1"])
PGM.nhc <- summaryGspline1(x=grid, mu=knots, sigma=sigma, standard=FALSE, 
                           intcpt=shift.nhc, scale=scale.nhc,
                           probs=c(0.025, 0.25, 0.5, 0.75, 0.975), values=TRUE, dir=dirPaths["PGM_nhc"])


###################################################
### chunk number 34: unstandardized.random.intercept.density-PGM-hc eval=FALSE
###################################################
shift.hc <- scanFH(paste(dirPaths["PGM_hc"], "betaR.sim", sep = ""))[, "(Intercept)"]
scale.hc <- sqrt(scanFH(paste(dirPaths["PGM_hc"], "varR.sim", sep = ""))[, "varR.1"])
PGM.hc <- summaryGspline1(x=grid, mu=knots, sigma=sigma, standard=FALSE, 
                          intcpt=shift.hc, scale=scale.hc,
                          probs=c(0.025, 0.25, 0.5, 0.75, 0.975), values=TRUE, dir=dirPaths["PGM_hc"])


###################################################
### chunk number 35: figure-density-PGMnhc eval=FALSE
###################################################
par(bty="n", mar=c(4, 4, 1, 1)+0.1)
layout(matrix(c(1,1,2,2, 3,4,5,6, 7,8,9,10, 11,12,13,14), ncol=4, byrow=TRUE))

### Posterior mean
plot(PGM.nhc$summary$x, PGM.nhc$summary$Mean, type="l", xlab="b", ylab="g(b)", col="blue")

### Posterior median and 2.5%, 97.5% quantiles
plot(PGM.nhc$summary$x, PGM.nhc$summary[,"97.5%"], type="l", lty=2, xlab="b", ylab="g(b)", col="red")
lines(PGM.nhc$summary$x, PGM.nhc$summary[,"2.5%"], lty=2, col="red")
lines(PGM.nhc$summary$x, PGM.nhc$summary[,"50%"], lty=1, col="blue")

### Sampled densities at selected iterations
ylim <- c(0, 1)
#for (iters in c(1, 100, 200, 300, 400, 500, 600, 700, 800, 900, 950, 1000)){
for (iters in c(1, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 1900, 2000)){
  plot(PGM.nhc$summary$x, PGM.nhc$values[iters,], type="l", xlab="", ylab="", col="darkblue", ylim=ylim)
  title(main=paste("Iteration ", nsimul$nburn+iters, sep=""))
}


###################################################
### chunk number 36: extract-id-from-data
###################################################
IDNR <- unique(toenail$idnr)
IDNR0 <- unique(subset(toenail, trt==0)$idnr)
IDNR1 <- unique(subset(toenail, trt==1)$idnr)
index.tr0 <- (1:length(IDNR))[IDNR %in% IDNR0]
index.tr1 <- (1:length(IDNR))[IDNR %in% IDNR1]


###################################################
### chunk number 37: read-b-in-PGM-GLMM-nhc eval=FALSE
###################################################
beta1.PGMnhc <- scanFH(paste(dirPaths["PGM_nhc"], "betaF.sim", sep = ""))[,"(Intercept)"]
b.PGMnhc <- scanFH(paste(dirPaths["PGM_nhc"], "b.sim", sep = "")) + beta1.PGMnhc
colnames(b.PGMnhc) <- IDNR


###################################################
### chunk number 38: mean-median-b-in-PGM-GLMM-nhc eval=FALSE
###################################################
bMean.PGMnhc <- apply(b.PGMnhc, 2, mean)
bMedian.PGMnhc <- apply(b.PGMnhc, 2, median)


###################################################
### chunk number 39: figure-posterior-mean-and-median-of-individual-random-effects-in-PGM-GLMM-nhc eval=FALSE
###################################################
xlim <- c(-7, 11)
ylim <- c(0, 0.6)

layout(matrix(c(0,1,1,0, 2,2,3,3, 0,4,4,0, 5,5,6,6), ncol=4, byrow=TRUE))
par(bty="n", mar=c(4, 4, 4, 1)+0.1)

### Histograms of posterior means
hist(bMean.PGMnhc, prob=TRUE, xlab="beta1+b", ylab="Density", col="seagreen3", xlim=xlim, ylim=ylim, main="Posterior mean (all patients)")
hist(bMean.PGMnhc[index.tr0], prob=TRUE, xlab="beta1+b", ylab="Density", col="skyblue4", xlim=xlim, ylim=ylim, main="Control")
hist(bMean.PGMnhc[index.tr1], prob=TRUE, xlab="beta1+b", ylab="Density", col="skyblue4", xlim=xlim, ylim=ylim, main="Treatment")

### Histograms of posterior medians
hist(bMedian.PGMnhc, prob=TRUE, xlab="beta1+b", ylab="Density", col="seagreen3", xlim=xlim, ylim=ylim, main="Posterior median (all patients)")
hist(bMedian.PGMnhc[index.tr0], prob=TRUE, xlab="beta1+b", ylab="Density", col="skyblue4", xlim=xlim, ylim=ylim, main="Control")
hist(bMedian.PGMnhc[index.tr1], prob=TRUE, xlab="beta1+b", ylab="Density", col="skyblue4", xlim=xlim, ylim=ylim, main="Treatment")


