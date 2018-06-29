###
load("result.RData")  


thetaZ<-apply(thetaZ.iter[-c(6:7),],1,quantile,probs=c(.025,.5,.975)) 

thetaX<-apply(thetaX.iter,1,quantile,probs=c(.025,.5,.975))


colnames(thetaZ)<-c("alpha1","alpha2","beta1","beta2","beta3","h") 

colnames(thetaX)<-c("b","a","c1","c2") 


thetaZ<-thetaZ[,c(3:5,1,2,6)] 


apply(thetaZ1.iter,1,quantile,probs=c(.025,.5,.975))  
apply(thetaX1.iter,1,quantile,probs=c(.025,.5,.975))  


acf(thetaZ.iter[8,] )$acf[11]  


pdf("hfmix.pdf",height=3,width=8,family="Times")
par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot( t(thetaZ.iter)[,8] ,type="l",xlab="iteration",ylab=expression(h)) 
abline(h=median(t(thetaZ.iter)[,8] ),lty=2) 
dev.off()

