library(amen)
library(coda)
load("result.RData")
thetaX.thin<-thetaX.iter[,1:1000*20]
matplot(t(thetaX.thin),type="l")

thetaX.mean<-apply(thetaX.thin,1,mean)

thetaX.sd<-apply(thetaX.thin,1,sd)
thetaX.eff<-apply(thetaX.thin,1,effectiveSize)
thetaX.se<-thetaX.sd/sqrt(thetaX.eff)

thetaX.confint<-apply(thetaX.thin,1,quantile,c(0.025,0.975))

##
thetaX1.thin<-thetaX1.iter[,1:1000*20,drop=F]
matplot(t(thetaX1.thin),type="l")

thetaX1.mean<-apply(thetaX1.thin,1,mean)

thetaX1.sd<-apply(thetaX1.thin,1,sd)
thetaX1.eff<-apply(thetaX1.thin,1,effectiveSize)
thetaX1.se<-thetaX1.sd/sqrt(thetaX1.eff)

thetaX1.confint<-apply(thetaX1.thin,1,quantile,c(0.025,0.975))

##
thetaZ.thin<-thetaZ.iter[,1:1000*20,drop=F]
matplot(t(thetaZ.thin),type="l")

thetaZ.mean<-apply(thetaZ.thin,1,mean)

thetaZ.sd<-apply(thetaZ.thin,1,sd)
thetaZ.eff<-apply(thetaZ.thin,1,effectiveSize)
thetaZ.se<-thetaZ.sd/sqrt(thetaZ.eff)

thetaZ.confint<-apply(thetaZ.thin,1,quantile,c(0.025,0.975))

##
thetaZ1.thin<-thetaZ1.iter[,1:1000*20,drop=F]
matplot(t(thetaZ1.thin),type="l")

thetaZ1.mean<-apply(thetaZ1.thin,1,mean)

thetaZ1.sd<-apply(thetaZ1.thin,1,sd)
thetaZ1.eff<-apply(thetaZ1.thin,1,effectiveSize)
thetaZ1.se<-thetaZ1.sd/sqrt(thetaZ1.eff)

thetaZ1.confint<-apply(thetaZ1.thin,1,quantile,c(0.025,0.975))
