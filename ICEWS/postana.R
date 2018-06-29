load("mcmc_output")
load("YICEWS") 

m<-dim(MPS)[1] 
n<-dim(XHX)[3] 
r<-ncol(YPAR)-2   

cnames<-rownames(MPS) 

## -- effective sample sizes
nburn<-2500 ; odens<-10
burned<-(1:nrow(YPAR))[-(1:(nburn/odens))] 
nsscan<-length( burned )

range( coda::effectiveSize(  cbind(XPAR[burned,], YPAR[burned,]) ) )
mean(  coda::effectiveSize(  cbind(XPAR[burned,], YPAR[burned,]) ))


## -- some plots 
r<-ncol(YPAR)-2 
par(mfrow=c(3,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(YPAR[,1],type="l",ylab="Y autocor")
matplot(YPAR[,1+1:r],type="l",ylab="homophily")
matplot(XPAR[,1:r^2],type="l",ylab="X autocor")
matplot(XPAR[,r^2+1:r^2],type="l",ylab="contagion")
plot(YPAR[,2+r],type="l",ylab=expression(sigma^2))
matplot(XSS,type="l")


## -- Y model posterior quantiles 
apply(YPAR[burned,],2,quantile,prob=c(.025,.5,.975)) 


## -- relative SS
apply( sweep( SSDY,1,apply(SSDY,1,sum),"/")[burned,]  ,2,mean) 


## -- X model posterior quantities
round( apply(XPAR[burned,],2,quantile,prob=c(.025,.5,.975)) ,3) 

## -- relative SS
apply( sweep( SSDX,1,apply(SSDX,1,sum),"/")[burned,]  ,2,mean) 





## -- country plots 
XPM<-array(0,dim=c(m,r,n)) 
HP<-NULL
for(i in 1:n)
{
  eX<-eigen(XHX[,,i]/nsscan)
  XPM[,,i]<-eX$vec[,1:r]
  XPM[,,i]<-sweep(XPM[,,i],2,sign(XPM[match("United States",cnames) ,,i]),"*")
  HP<-rbind(HP,eX$val[1:r]) 
  XPM[,,i]<-XPM[,,i]%*%diag( sqrt(eX$val[1:r]) )
}


pdf("tsattributes.pdf",height=6,width=8,family="Times") 

par(mfrow=c(2,1),mar=c(1,3,1,1),mgp=c(1.75,.75,0), oma=c(4,1,1,1)) 

cs<-isonames<-list()
cs[[1]]<-match( c("United States","Iran","United Kingdom"), cnames)  
cs[[2]]<-match( c("Russian Federation","Ukraine","Germany"), cnames)
isonames[[1]]<-c("USA","IRN","UKG")  
isonames[[2]]<-c("RUS","UKR","GER") 

cols<-list( c("red","green","blue"), c("red","lightblue","orange")) 

for(k in 1:2)
{
matplot(t(XPM[cs[[k]],k,-c(1,n)] ) ,type="l",lty=1,lwd=3,xlim=c(0,n+1),
    ylab="",xlab="",xaxt="n",col=cols[[k]] )
if(k==2) 
axis(1,at=seq(1,n,by=15),labels=dimnames(Y)[[3]][-c(1,n)][seq(1,n-2,by=15)],
     las=3 )
text(rep(n+2,length(cs[[k]])),XPM[cs[[k]],k,n-1],labels=isonames[[k]])
text(rep(-1,length(cs[[k]])),XPM[cs[[k]],k,2],labels=isonames[[k]])
abline(h=0)

}
dev.off() 



ACFX<-apply(XPAR[burned,],2,function(x){ acf(x)$acf[2]  } ) 
ACFY<-apply(YPAR[burned,],2,function(x){ acf(x)$acf[2]  } )

range(ACFX)
mean(ACFX)

range(ACFY)
mean(ACFY)

range(c(ACFX,ACFY))
mean(c(ACFX,ACFY))


#### ---- 
pdf("hmix.pdf",height=3,width=8,family="Times")
par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
matplot(YPAR[burned,2:3],type="l",lty=1,xlab="iteration",ylab="")
abline(h=apply(YPAR[burned,2:3],2,median) ,col=1:2,lty=2)  
legend(1750,.0675,legend=c(expression(h[11]),expression(h[22])),
       lty=c(1,1),col=c("black","red"),bty="n") 
dev.off()


