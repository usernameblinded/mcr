# n for time and M for number of nodes
library(MCMCpack)
library(amen) # need to use mvrnorm and zscores
library(truncnorm)

rXIinv_fc<-function(X,lambda){
  ed<-eigen(X)
  d<-1/(ed$values+lambda)
  ev<-ed$vectors
  inv<-tcrossprod(tcrossprod(ev,diag(d)),ev)
  return(inv)
}

norm2<-function(x){
  return(sum((x)^2))
}

rZcoef_fc<-function(Z,Xdyad,X,s2,transpose=T,rc=T,mu0=c(1,0,rep(0,dim(Xdyad)[3]+dim(X)[2]*2+dim(X)[2]^2)),var0=100){
  M<-dim(Z)[1]
  N<-dim(Z)[3]
  p<-dim(X)[2]
  npara<-length(mu0)
  ndyad<-dim(Xdyad)[3]
  y<-c(Z[,,2:N])
  xy<-c(Z[,,2:N-1])
  xyt<-c(aperm(Z[,,2:N-1],c(2,1,3)))
  xdyad<-rep(1,N-1)%x%matrix(c(Xdyad),ncol=ndyad)
  xr<-xc<-xm<-NULL
  for (n in 2:N-1){
    xr<-rbind(xr,matrix(rep(1,M)%x%X[,,n],ncol=p))
    xc<-rbind(xc,matrix(X[,,n]%x%rep(1,M),ncol=p))
    xm<-rbind(xm,matrix(X[,,n]%x%X[,,n],ncol=p^2))
  }
  Design<-cbind(xy,xyt,xdyad,xr,xc,xm)
  columns<-rep(TRUE,npara)
  if(!transpose){
    columns[2]<-F
  }
  if (!rc){
    columns[1:(2*p)+ndyad+2]<-F
  }
  Design<-Design[,columns]
  mu0<-mu0[columns]
  Design<-Design[!is.na(y),]
  y<-matrix((y[!is.na(y)]),ncol=1)
  Design.cross<-crossprod(Design)
  Q<-Design.cross/s2+1/var0*diag(ncol(Design))
  L<-crossprod(Design,y)/s2+1/var0*mu0
  Sigma<-s2*rXIinv_fc(Design.cross,s2/var0)
  Mu<-crossprod(Sigma,L)
  theta<-rmvnorm(1,Mu,Sigma)
  theta2<-rep(0,npara)
  theta2[columns]<-theta
  return(theta2)
}

rZ1coef_fc<-function(Z,Xdyad,s1,mu0=rep(0,dim(Xdyad)[3]),var0=100){
  M<-dim(Z)[1]
  N<-dim(Z)[3]
  ndyad<-dim(Xdyad)[3]
  y<-c(Z[,,1])
  Design<-matrix(c(Xdyad),ncol=ndyad)
  Design<-Design[!is.na(y),,drop=F]
  y<-matrix((y[!is.na(y)]),ncol=1)
  Design.cross<-crossprod(Design)
  Q<-Design.cross/s1+1/var0*diag(1,ncol=ncol(Design),nrow=ncol(Design))
  L<-crossprod(Design,y)/s1+1/var0*mu0
  Sigma<-solve(Q)
  Mu<-crossprod(Sigma,L)
  theta<-rmvnorm(1,Mu,Sigma)
  return(theta)
}

rsigma1_fc<-function(Z,Xdyad,thetaZ1,gamma.prior=c(3,1/3)){
  M<-dim(Z)[1]
  N<-dim(Z)[3]
  Z10<-Z[,,1]-apply(sweep(Xdyad,3,thetaZ1,"*"),1:2,sum); diag(Z10)<-NA
  W<-sum(Z10^2,na.rm=T)
  sigma2<-rgamma(1,M*(M-1)/2+gamma.prior[1],prod(gamma.prior)+W/2)
  sigma2<-1/sigma2
  return(sigma2)
}

rXcoef_fc<-function(X,Z,Xdyad,t2,transpose=T,mu0=NULL,var0=100){
  M<-dim(X)[1]
  p<-dim(X)[2]
  N<-dim(X)[3]
  npara<-3*p*p+p
  sex<-Xdyad[,1,1]
  if (is.null(mu0)){
    mu0<-c(rep(0,p),diag(1,ncol=p,nrow=p),rep(0,2*p*p))
  }
  y<-matrix(c(X[,,2:N]))
  x<-zx<-ztx<-NULL
  for (n in 2:N-1){
    Zn0<-Z[,,n];
    diag(Zn0)<-0
    x<-rbind(x,cbind(diag(p)%x%matrix(sex,ncol=1),diag(p)%x%X[,,n]))
    zx<-rbind(zx,diag(p)%x%(Zn0%*%X[,,n]))
    ztx<-rbind(ztx,diag(p)%x%crossprod(Zn0,X[,,n]))
  }
  Design<-cbind(x,zx,ztx)
  columns<-rep(T,npara)
  if (!transpose){
    columns[1:(p^2)+2*p*p+p]<-F
  }
  Design<-Design[,columns]
  mu0<-mu0[columns]
  Design.cross<-crossprod(Design)
  Q<-Design.cross/t2+1/var0*diag(1,ncol=ncol(Design),nrow=ncol(Design))
  L<-crossprod(Design,y)/t2+1/var0*mu0
  Sigma<-t2*rXIinv_fc(Design.cross,t2/var0)
  Mu<-crossprod(Sigma,L)
  theta<-rmvnorm(1,Mu,Sigma)
  theta2<-rep(0,npara)
  theta2[columns]<-theta
  return(theta2)
}

rX1coef_fc<-function(X,Xdyad,t1,mu0=0,var0=100){
  p<-dim(X)[2]
  y<-c(X[,,1])
  Design<-diag(p)%x%Xdyad[,1,1]
  Design<-Design[!is.na(y),,drop=F]
  y<-matrix((y[!is.na(y)]),ncol=1)
  Design.cross<-crossprod(Design)
  Q<-Design.cross/t1+1/var0*diag(1,ncol=ncol(Design),nrow=ncol(Design))
  L<-crossprod(Design,y)/t1+1/var0*mu0
  Sigma<-solve(Q)
  Mu<-crossprod(Sigma,L)
  theta<-rmvnorm(1,Mu,Sigma)
  return(theta)
}

rtau1_fc<-function(X,Xdyad,thetaX1,gamma.prior=c(3,0.2)){
  M<-dim(X)[1]
  p<-dim(X)[2]
  tau2<-NULL
  for (j in 1:p){
    W<-sum((X[,j,1]-Xdyad[,1,1]*thetaX1[j])^2)
    tau2[j]<-rgamma(1,M/2+gamma.prior[1],prod(gamma.prior)+W/2)
    tau2[j]<-1/tau2[j]
  }
  
  return(tau2)
}

rZ_ord_fc<-function(Z,Xdyad,X,thetaZ,thetaZ1,thetaX,thetaX1,Y,thresZ,s1,t2,s2,mean.prior=0,var.prior=100){
  label<-is.na(Y)
  M<-dim(Z)[1]
  N<-dim(Z)[3]
  ranksY<-seq(min(Y,na.rm=T),max(Y,na.rm=T))
  Nd<-dim(Xdyad)[3]
  p<-dim(X)[2]
  Xdyad.prod1<-apply(sweep(Xdyad,3,thetaZ1,"*"),1:2,sum)
  Xdyad.prod2<-apply(sweep(Xdyad,3,thetaZ[1:Nd+2],"*"),1:2,sum)
  Zm<-matrix(thetaZ[1:(p^2)+2+Nd+2*p],ncol=p)
  Zr<-thetaZ[1:p+2+Nd]
  Zc<-thetaZ[1:p+2+Nd+p]
  sex<-matrix(Xdyad[,1,1],ncol=1)
  
  G<-A<-B<-array(dim=dim(X))
  G[,,1]<-X[,,1]-sweep(matrix(rep(sex,p),ncol=p),2,thetaX1,"*")
  A[,,1]<-X[,,1]%*%matrix(thetaX[1:(p^2)+p^2+p],ncol=p)
  B[,,1]<-X[,,1]%*%matrix(thetaX[1:(p^2)+2*p^2+p],ncol=p)
  for (n in 2:N){
    G[,,n]<-X[,,n]-sweep(matrix(rep(sex,p),ncol=p),2,thetaX[1:p],"*")-
            X[,,n-1]%*%matrix(thetaX[1:(p^2)+p],ncol=p)
    A[,,n]<-tcrossprod(X[,,n],matrix(thetaX[1:(p^2)+p^2+p],ncol=p))
    B[,,n]<-tcrossprod(X[,,n],matrix(thetaX[1:(p^2)+2*p^2+p],ncol=p))
  }
  
  for (n in 1:N){
    diag(Z[,,n])<-0
  }
  
  
  for (n in 1:N){
    for (i in 1:M){
      for (j in 1:M){
        if (i!=j){
          if (n==1){
            sig<-1/(1/var.prior+1/s1+(thetaZ[1]^2+thetaZ[2]^2)/s2+sum(A[j,,n]^2)/t2+sum(B[i,,n]^2)/t2)
            L<-mean.prior/var.prior+Xdyad.prod1[i,j]/s1+
              thetaZ[1]*(Z[i,j,n+1]-thetaZ[2]*Z[j,i,n]-Xdyad.prod2[i,j]-sum(X[i,,n]*Zr)-sum(X[j,,n]*Zc)-sum((X[i,,n]%*%Zm)*X[j,,n]))/s2+
              thetaZ[2]*(Z[j,i,n+1]-thetaZ[1]*Z[i,j,n]-Xdyad.prod2[j,i]-sum(X[j,,n]*Zr)-sum(X[i,,n]*Zc)-sum((X[j,,n]%*%Zm)*X[i,,n]))/s2+
              sum(A[j,,n]*(G[i,,n+1]-Z[i,,n]%*%A[,,n]-t(Z[,i,n])%*%B[,,n]+Z[i,j,n]*A[j,,n]+Z[i,i,n]*A[i,,n]+Z[i,i,n]*B[i,,n]))/t2+
              sum(B[i,,n]*(G[j,,n+1]-Z[j,,n]%*%A[,,n]-t(Z[,j,n])%*%B[,,n]+Z[i,j,n]*B[i,,n]+Z[j,j,n]*A[j,,n]+Z[j,j,n]*B[j,,n]))/t2
            ez<-sig*L
            sdz<-sqrt(sig)
            if (label[i,j,n]){
              Z[i,j,n]<-rnorm(1,ez,sdz)
            }
            if (!label[i,j,n]){
              rankid<-which(ranksY==Y[i,j,n])
              Z[i,j,n]<-rtruncnorm(1,thresZ[rankid],thresZ[rankid+1],ez,sdz)
            }
          }
          
          if (n>1 & n<N){
            sig<-1/(1/var.prior+1/s2+(thetaZ[1]^2+thetaZ[2]^2)/s2+sum(A[j,,n]^2)/t2+sum(B[i,,n]^2)/t2)
            L<-mean.prior/var.prior+(thetaZ[1]*Z[i,j,n-1]+thetaZ[2]*Z[j,i,n-1]+Xdyad.prod2[i,j]+sum(X[i,,n-1]*Zr)+sum(X[j,,n-1]*Zc)+sum((X[i,,n-1]%*%Zm)*X[j,,n-1]))/s2+
              thetaZ[1]*(Z[i,j,n+1]-thetaZ[2]*Z[j,i,n]-Xdyad.prod2[i,j]-sum(X[i,,n]*Zr)-sum(X[j,,n]*Zc)-sum((X[i,,n]%*%Zm)*X[j,,n]))/s2+
              thetaZ[2]*(Z[j,i,n+1]-thetaZ[1]*Z[i,j,n]-Xdyad.prod2[j,i]-sum(X[j,,n]*Zr)-sum(X[i,,n]*Zc)-sum((X[j,,n]%*%Zm)*X[i,,n]))/s2+
              sum(A[j,,n]*(G[i,,n+1]-Z[i,,n]%*%A[,,n]-t(Z[,i,n])%*%B[,,n]+Z[i,j,n]*A[j,,n]+Z[i,i,n]*A[i,,n]+Z[i,i,n]*B[i,,n]))/t2+
              sum(B[i,,n]*(G[j,,n+1]-Z[j,,n]%*%A[,,n]-t(Z[,j,n])%*%B[,,n]+Z[i,j,n]*B[i,,n]+Z[j,j,n]*A[j,,n]+Z[j,j,n]*B[j,,n]))/t2
            ez<-sig*L
            sdz<-sqrt(sig)
            if (label[i,j,n]){
              Z[i,j,n]<-rnorm(1,ez,sdz)
            }
            if (!label[i,j,n]){
              rankid<-which(ranksY==Y[i,j,n])
              Z[i,j,n]<-rtruncnorm(1,thresZ[rankid],thresZ[rankid+1],ez,sdz)
            }
          }
          
          if (n==N){
            sig<-1/(1/var.prior+1/s2)
            L<-mean.prior/var.prior+(thetaZ[1]*Z[i,j,n-1]+thetaZ[2]*Z[j,i,n-1]+Xdyad.prod2[i,j]+sum(X[i,,n-1]*Zr)+sum(X[j,,n-1]*Zc)+sum((X[i,,n-1]%*%Zm)*X[j,,n-1]))/s2
            ez<-sig*L
            sdz<-sqrt(sig)
            if (label[i,j,n]){
              Z[i,j,n]<-rnorm(1,ez,sdz)
            }
            if (!label[i,j,n]){
              rankid<-which(ranksY==Y[i,j,n])
              Z[i,j,n]<-rtruncnorm(1,thresZ[rankid],thresZ[rankid+1],ez,sdz)
            }
          } 
        }
      }
    }
  }
  
  for (n in 1:N){
    diag(Z[,,n])<-NA
  }

  return (Z)
}

rX_ord_fc<-function(Z,Xdyad,X,thetaZ,thetaX,thetaX1,W,thresX,t1,t2,s2,mu0=0,var0=100){
  label<-is.na(W)
  ranksW<-list()
  
  M<-dim(Z)[1]
  N<-dim(Z)[3]
  Nd<-dim(Xdyad)[3]
  p<-dim(X)[2]
  sex<-matrix(Xdyad[,1,1],ncol=1)
  
  for (n in 1:N){
    diag(Z[,,n])<-NA
  }
  
  Zm<-matrix(thetaZ[1:(p^2)+2+Nd+2*p],ncol=p)
  
  for (j in 1:p){
    ranksW[[j]]<-seq(min(W[,j,],na.rm=T),max(W[,j,],na.rm=T))
  }

  Ax<-t(matrix(thetaX[1:(p^2)+p],ncol=p))
  Bx<-t(matrix(thetaX[1:(p^2)+p^2+p],ncol=p))
  Cx<-t(matrix(thetaX[1:(p^2)+2*p^2+p],ncol=p))
  
  Sz<-Z[,,2:N]-thetaZ[1]*Z[,,2:N-1]-thetaZ[2]*aperm(Z[,,2:N-1],c(2,1,3))-
    array(rep(apply(sweep(Xdyad,3,thetaZ[1:Nd+2],"*"),1:2,sum),N-1),dim=c(M,M,N-1))
  
  n<-1
  Zn0<-Z[,,n]; diag(Zn0)<-0
  for (i in 1:M){
    for (j in 1:p){
      Sr<-Sc<-Sm<-matrix(0,nrow=M,ncol=M)
      Sr[i,]<-thetaZ[Nd+2+j]
      Sc[,i]<-thetaZ[Nd+2+p+j]
      Sm[i,]<-rowSums(sweep(matrix(X[,,n],ncol=p),2,c(Zm[j,]),"*"))
      Sm[,i]<-rowSums(sweep(matrix(X[,,n],ncol=p),2,c(Zm[,j]),"*"))
      Srcm<-Sr+Sc+Sm
      diag(Srcm)<-NA
      Sij<-Sz[,,n]-tcrossprod(X[,,n],outer(rep(1,M),thetaZ[(Nd+3):(Nd+2+p)]))-
        tcrossprod(outer(rep(1,M),thetaZ[(Nd+3+p):(Nd+2+2*p)]),X[,,n])-
        tcrossprod(X[,,n],tcrossprod(X[,,n],matrix(thetaZ[(Nd+3+2*p):(Nd+2+2*p+p^2)],nrow=p)))+
        X[i,j,n]*(Srcm)
      
      ABCx<-outer(diag(M)[,i],Ax[j,])+outer(Zn0[,i],Bx[j,])+outer(Zn0[i,],Cx[j,])
      Rx<-X[,,n+1]-sweep(matrix(rep(sex,p),ncol=p),2,thetaX[1:p],"*")-
          X[,,n]%*%Ax-Zn0%*%X[,,n]%*%Bx-crossprod(Zn0,X[,,n]%*%Cx)
      Rx<-Rx+X[i,j,n]*ABCx
      
      sigma2<-1/(1/t1[j]+sum(ABCx^2)/t2+sum(Srcm^2,na.rm=T)/s2+1/var0)
      L<-sex[i]*thetaX1[j]/t1[j]+sum(ABCx*Rx)/t2+sum(Srcm*Sij,na.rm=T)/s2+mu0/var0
      ex<-sigma2*L
      sdx<-sqrt(sigma2)
      
      if (label[i,j,n]){
        X[i,j,n]<-rnorm(1,ex,sdx)
      }
      if (!label[i,j,n]){
        rankid<-which(ranksW[[j]]==W[i,j,n])
        X[i,j,n]<-rtruncnorm(1,thresX[[j]][rankid],thresX[[j]][rankid+1],ex,sdx)
      }
    }
  }
  
  for (n in 2:(N-1)){
    Zn0<-Z[,,n]; diag(Zn0)<-0
    for (i in 1:M){
      for (j in 1:p){
        Sr<-Sc<-Sm<-matrix(0,nrow=M,ncol=M)
        Sr[i,]<-thetaZ[Nd+2+j]
        Sc[,i]<-thetaZ[Nd+2+p+j]
        Sm[i,]<-rowSums(sweep(matrix(X[,,n],ncol=p),2,c(Zm[j,]),"*"))
        Sm[,i]<-rowSums(sweep(matrix(X[,,n],ncol=p),2,c(Zm[,j]),"*"))
        Srcm<-Sr+Sc+Sm
        diag(Srcm)<-NA
        Sij<-Sz[,,n]-tcrossprod(X[,,n],outer(rep(1,M),thetaZ[(Nd+3):(Nd+2+p)]))-
          tcrossprod(outer(rep(1,M),thetaZ[(Nd+3+p):(Nd+2+2*p)]),X[,,n])-
          tcrossprod(X[,,n],tcrossprod(X[,,n],matrix(thetaZ[(Nd+3+2*p):(Nd+2+2*p+p^2)],nrow=p)))+
          X[i,j,n]*(Srcm)
        
        ABCx<-outer(diag(M)[,i],Ax[j,])+outer(Zn0[,i],Bx[j,])+outer(Zn0[i,],Cx[j,])
        Rx<-X[,,n+1]-sweep(matrix(rep(sex,p),ncol=p),2,thetaX[1:p],"*")-
          X[,,n]%*%Ax-Zn0%*%X[,,n]%*%Bx-crossprod(Zn0,X[,,n]%*%Cx)
        Rx<-Rx+X[i,j,n]*ABCx
        
        Rij<-sex[i]*thetaX[j]+X[i,,n-1]%*%Ax[,j]+Zn0[,i]%*%X[,,n-1]%*%Bx[,j]+t(Zn0[,i])%*%X[,,n-1]%*%Cx[,j]
        
        sigma2<-1/(1/t2+sum(ABCx^2)/t2+sum(Srcm^2,na.rm=T)/s2+1/var0)
        L<-Rij/t2+sum(ABCx*Rx)/t2+sum(Srcm*Sij,na.rm=T)/s2+mu0/var0
        ex<-sigma2*L
        sdx<-sqrt(sigma2)
        
        if (label[i,j,n]){
          X[i,j,n]<-rnorm(1,ex,sdx)
        }
        if (!label[i,j,n]){
          rankid<-which(ranksW[[j]]==W[i,j,n])
          X[i,j,n]<-rtruncnorm(1,thresX[[j]][rankid],thresX[[j]][rankid+1],ex,sdx)
        }
      }
    }
  }
  
  n<-N
  for (i in 1:M){
    for (j in 1:p){
      Rij<-sex[i]*thetaX[j]+X[i,,n-1]%*%Ax[,j]+Zn0[,i]%*%X[,,n-1]%*%Bx[,j]+t(Zn0[,i])%*%X[,,n-1]%*%Cx[,j]
      
      sigma2<-1/(1/t2+1/var0)
      L<-Rij/t2+mu0/var0
      ex<-sigma2*L
      sdx<-sqrt(sigma2)
      
      if (label[i,j,n]){
        X[i,j,n]<-rnorm(1,ex,sdx)
      }
      if (!label[i,j,n]){
        rankid<-which(ranksW[[j]]==W[i,j,n])
        X[i,j,n]<-rtruncnorm(1,thresX[[j]][rankid],thresX[[j]][rankid+1],ex,sdx)
      }
    }
  }
  return(X)
}

rthreshold_fc<-function(Z,Y,mu0=0,sd0=100){
  ranks<-seq(min(Y,na.rm=T),max(Y,na.rm=T))
  thres<-NULL
  thres[1]<- -Inf
  for (i in 2:length(ranks)){
    a<-max(Z[Y==ranks[i-1]],na.rm=T)
    b<-min(Z[Y==ranks[i]],na.rm=T)
    u<-runif(1,pnorm((a-mu0)/sd0),pnorm((b-mu0)/sd0))
    thres[i]<- mu0 + sd0*qnorm(u)
  }
  thres[length(ranks)+1]<-Inf
  return(thres)
}

main_mcmc<-function(Y,W,Xdyad,niter=5000,odens=25,burnin=2000,thin=25,Zt=F,Zrc=F,Xt=F,print=TRUE,plot=TRUE,file.name="ord_data.RData"){
  M<-dim(Y)[1]
  N<-dim(Y)[3]
  p<-dim(W)[2]
  Nd<-dim(Xdyad)[3]
  
  Z<-array(zscores(Y),dim=c(M,M,N))
  Zmean<-apply(Z,1:2,mean,na.rm=T)
  for (n in 1:N){
    Z.obs<-Z[,,n]
    Z.obs[is.na(Z.obs)]<-Zmean[is.na(Z.obs)]
    Z.obs.row<-rowMeans(Z.obs,na.rm=T)
    Z.obs.col<-colMeans(Z.obs,na.rm=T)
    ZA<-outer(Z.obs.row,Z.obs.col,"+")/2
    Z.obs[is.na(Z.obs)]<-ZA[is.na(Z.obs)]
    Z[,,n]<-Z.obs
    diag(Z[,,n])<-NA
  }
  
  X<-W
  for (i in 1:p){
    Xi<-array(zscores(W[,i,,drop=F]),dim=c(M,1,N))
    Ximean<-apply(Xi,1:2,mean,na.rm=T)
    for (n in 1:N){
      Xi.obs<-Xi[,,n,drop=F]
      Xi.obs[is.na(Xi.obs)]<-Ximean[is.na(Xi.obs)]
      Xi.obs.row<-rowMeans(Xi.obs,na.rm=T)
      Xi.obs.col<-colMeans(Xi.obs,na.rm=T)
      XiA<-outer(Xi.obs.row,Xi.obs.col,"+")/2
      Xi.obs[is.na(Xi.obs)]<-XiA[is.na(Xi.obs)]
      X[,i,n]<-Xi.obs
    }
  }
  
  sigma2<-1
  tau2<-1
  
  sigma1<-1
  tau1<-rep(1,p)
  
  Z.mean<-X.mean<-0
  thetaZ.iter<-NULL
  thetaZ1.iter<-NULL
  thetaX.iter<-NULL
  thetaX1.iter<-NULL
  tau2.iter<-sigma2.iter<-tau1.iter<-sigma1.iter<-NULL
  thresZ.iter<-NULL
  thresX.iter<-NULL
  
  for (it in 1:burnin){
    thresZ<-rthreshold_fc(Z,Y)
    thresX<-list()
    for (j in 1:p){
      thresX[[j]]<-rthreshold_fc(X[,j,],W[,j,])
    }
    
    thetaZ<-rZcoef_fc(Z,Xdyad,X,sigma2,Zt,Zrc)
    thetaZ1<-rZ1coef_fc(Z,Xdyad,sigma1)
    thetaX<-rXcoef_fc(X,Z,Xdyad,tau2,Xt)
    thetaX1<-rX1coef_fc(X,Xdyad,tau1)
    
    tau1<-rtau1_fc(X,Xdyad,thetaX1)
    sigma1<-rsigma1_fc(Z,Xdyad,thetaZ1)
    
    Z<-rZ_ord_fc(Z,Xdyad,X,thetaZ,thetaZ1,thetaX,thetaX1,Y,thresZ,sigma1,tau2,sigma2)
    X<-rX_ord_fc(Z,Xdyad,X,thetaZ,thetaX,thetaX1,W,thresX,tau1,tau2,sigma2) 
    
    if (it %% 25 ==0){
      cat("burn in complete: ",it,"\n")
    }  
  }
  
  for (it in 1:niter){
    thresZ<-rthreshold_fc(Z,Y)
    thresX<-list()
    for (j in 1:p){
      thresX[[j]]<-rthreshold_fc(X[,j,],W[,j,])
    }
    
    thetaZ<-rZcoef_fc(Z,Xdyad,X,sigma2,Zt,Zrc)
    thetaZ1<-rZ1coef_fc(Z,Xdyad,sigma1)
    thetaX<-rXcoef_fc(X,Z,Xdyad,tau2,Xt)
    thetaX1<-rX1coef_fc(X,Xdyad,tau1)
    
    tau1<-rtau1_fc(X,Xdyad,thetaX1)
    sigma1<-rsigma1_fc(Z,Xdyad,thetaZ1)
    
    Z<-rZ_ord_fc(Z,Xdyad,X,thetaZ,thetaZ1,thetaX,thetaX1,Y,thresZ,sigma1,tau2,sigma2)
    X<-rX_ord_fc(Z,Xdyad,X,thetaZ,thetaX,thetaX1,W,thresX,tau1,tau2,sigma2) 
    
    Z.mean <- Z.mean+Z
    X.mean <- X.mean+X
    thetaZ.iter<-cbind(thetaZ.iter,c(thetaZ))
    thetaZ1.iter<-cbind(thetaZ1.iter,c(thetaZ1))
    thetaX.iter<-cbind(thetaX.iter,c(thetaX))
    thetaX1.iter<-cbind(thetaX1.iter,c(thetaX1))
    tau2.iter[it]<-tau2
    sigma2.iter[it]<-sigma2
    tau1.iter<-cbind(tau1.iter,tau1)
    sigma1.iter[it]<-sigma1
    thresZ.iter<-cbind(thresZ.iter,unlist(thresZ))
    thresX.iter<-cbind(thresX.iter,unlist(thresX))
    
    if (print & it %% odens ==0){
      cat(it,": ",max(X,na.rm=T),tau1,sigma1,tau2,sigma2,"\n","thetaX: ",thetaX,"\n", "thetaZ: ",thetaZ,"\n", "thetaX1: ",thetaX1,"\n", "thetaZ1: ",thetaZ1,"\n")
    }
    if (plot & it %% odens ==0){
      nsample<-floor(it/thin)
      par(mfrow=c(2,5))
      plot(thetaZ.iter[1,(1:nsample*thin)],type="l",ylab="Z_autoregressive")
      plot(thetaZ.iter[2,(1:nsample*thin)],type="l",ylab="Z_reciprocity")
      plot(thetaZ.iter[3,(1:nsample*thin)],type="l",ylab="Z_sex_row")
      plot(thetaZ.iter[2+Nd+2*p+1,(1:nsample*thin)],type="l",ylab="Homophily")
      plot(thetaX.iter[1,(1:nsample*thin)],type="l",ylab="X_sex")
      plot(thetaX.iter[2,(1:nsample*thin)],type="l",ylab="X_autoregressive")
      plot(thetaX.iter[3,(1:nsample*thin)],type="l",ylab="X_contagion1")
      plot(thetaX.iter[4,(1:nsample*thin)],type="l",ylab="X_contagion2")
      plot(thetaX1.iter[1,(1:nsample*thin)],type="l",ylab="X1_sex")
      plot(tau1.iter[1,(1:nsample*thin)],type="l",ylab="tau1")
    }
    if (it %% 1000 ==0){
      nsample<-floor(it/thin)
      jpeg(paste(file.name,".jpg",sep=""),w=1400,h=800)
      par(mfrow=c(2,5))
      plot(thetaZ.iter[1,(1:nsample*thin)],type="l",ylab="Z_autoregressive")
      plot(thetaZ.iter[2,(1:nsample*thin)],type="l",ylab="Z_reciprocity")
      plot(thetaZ.iter[3,(1:nsample*thin)],type="l",ylab="Z_sex_row")
      plot(thetaZ.iter[2+Nd+2*p+1,(1:nsample*thin)],type="l",ylab="Homophily")
      plot(thetaX.iter[1,(1:nsample*thin)],type="l",ylab="X_sex")
      plot(thetaX.iter[2,(1:nsample*thin)],type="l",ylab="X_autoregressive")
      plot(thetaX.iter[3,(1:nsample*thin)],type="l",ylab="X_contagion1")
      plot(thetaX.iter[4,(1:nsample*thin)],type="l",ylab="X_contagion2")
      plot(thetaX1.iter[1,(1:nsample*thin)],type="l",ylab="X1_sex")
      plot(tau1.iter[1,(1:nsample*thin)],type="l",ylab="tau1")
      dev.off()
      Zmean<-Z.mean/it
      Xmean<-X.mean/it
      save(niter,Y,W,Xdyad,Z,X,Zmean,Xmean,thetaZ.iter,thetaZ1.iter,thetaX.iter,thetaX1.iter,thresZ.iter,thresX.iter,tau1.iter,sigma1.iter,tau2.iter,sigma2.iter,file=file.name)
    }    
  }
}



