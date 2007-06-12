epilepticBC <- read.table("epileptic.dat", header=TRUE)

epilepticBC <- epilepticBC[order(epilepticBC$visit),]
epilepticBC <- epilepticBC[order(epilepticBC$id),]

epilepticBC$seizure0 <- rep(epilepticBC[seq(1, 5*59-4, by=5), "seizure"], each=5)
epilepticBC <- epilepticBC[-seq(1, 5*59-4, by=5),]

## Covariates as in the papers Breslow and Clayton (1993), Kleinman and Ibrahim (1998):
epilepticBC$Seizure  <- epilepticBC$seizure
epilepticBC$Base     <- log(epilepticBC$seizure0/4)
epilepticBC$Trt      <- epilepticBC$trt
epilepticBC$Base.Trt <- epilepticBC$Base*epilepticBC$Trt
epilepticBC$Age      <- log(epilepticBC$age)
epilepticBC$Visit    <- (2*epilepticBC$visit - 5)/10

epilepticBC <- subset(epilepticBC, select=c("id", "visit", "seizure0", "age", "Seizure", "Base", "Trt", "Base.Trt", "Age", "Visit"))

