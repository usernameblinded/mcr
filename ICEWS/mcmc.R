#### ---- data and functions
source("functions.R")
#m<-50 ; source("irdata.R")
load("YICEWS") 


#### ---- starting values and priors
r<-2 
m<-dim(Y)[1]
n<-dim(Y)[3]
Y[is.na(Y)]<-0

## -- regression model for Y
M<-apply(Y[,,-1],c(1,2),mean,na.rm=TRUE)  
fit<-lm(  c( sweep(Y[,,-1],c(1,2),M,"-") ) ~ -1 + c(Y[,,-n])) 
alpha<-fit$coef 
s2<-mean(fit$res^2)

## -- X's represent residual variation
X<-array(dim=c(m,r,n))

for(i in 1:(n-1))
{
  E<- Y[,,i+1] - M - alpha*Y[,,i] ; E[is.na(E)]<-0

  sE<-eigen(E)
  X[,,i]<-sE$vec[,1:r]
  if(i>1){ X[,,i]<-X[,,i]%*%
     diag(sign(apply(X[,,i,drop=F]*X[,,i-1,drop=F],2,mean)),nrow=r) }
}
X[,,n]<-X[,,n-1]

fit<-lm(  c(X[,,-1]) ~ -1 + c(X[,,-n]))
X<-X/sqrt(mean(fit$res^2)) 



#### ---- MCMC  
YPAR<-XPAR<-XSS<-NULL
XHX<-array(0,dim=dim(Y)) ; dimnames(XHX)[[1]]<-dimnames(Y)[[1]] 
MPS<-Y[,,1]*0 ; TPS<-X[,,1]*0

seed<-1 ; set.seed(seed)
plot<-TRUE
nscan<-25000 ; nburn<-2500 ; odens<-10
for(s in 1:(nscan+nburn))
{

  ## -- Gibbs updates 
  betaX<-rbetaX_fc(Y,X) ; T<-betaX$T ; A<-betaX$A ; C<-betaX$C  

  betaY<-rbetaY_fc(Y,X,s2) ; M<-betaY$M ; alpha<-betaY$alpha ; H<-betaY$H

  s2<-rs2_fc(Y,M,alpha,H,X)

  X<-rX_fc(Y,M,alpha,H,s2,X,T,A,C)

  ## -- save and plot 
  if(s%%odens==0) 
  {  

    cat(s," ",s2,"\n")

    if(s>nburn)
    {
      for(i in 1:n){ XHX[,,i]<-XHX[,,i] + X[,,i]%*%H%*%t(X[,,i]) }
      MPS<-MPS+M
      TPS<-TPS+T 
    } 

    YPAR<-rbind(YPAR,c(alpha,diag(H),s2))
    XPAR<-rbind(XPAR,c(A,C))
    XSS<-rbind(XSS,apply(X^2,2,mean))

    if(plot)
    {
      par(mfrow=c(3,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))  
      plot(YPAR[,1],type="l",ylab="Y autocor")
      matplot(YPAR[,1+1:r],type="l",ylab="homophily")
      matplot(XPAR[,1:r^2],type="l",ylab="X autocor")  
      matplot(XPAR[,r^2+1:r^2],type="l",ylab="contagion")  
      plot(YPAR[,2+r],type="l",ylab=expression(sigma^2))  
      matplot(XSS,type="l")
    } 

    if(s%%1000==0){ save(YPAR,XPAR,XSS,XHX,MPS,TPS,file="mcmc_output") }
  }

} 


 
