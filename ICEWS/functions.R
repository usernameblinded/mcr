library(rfunctions) 

rbetaX_fc<-function(Y,X,S=diag(nrow=ncol(X)),t2=1e2)
{
  m<-dim(Y)[1] ; n<-dim(Y)[3] ; r<-dim(X)[2]
  l<-rep(0,m*r+r^2+r^2)
  Q<-diag(m*r+r^2+r^2)/t2
  iV<-kron(diag(m),solve(S)) 

  Ir<-diag(nrow=r) 

  for(i in 2:n)
  {
    z<-c(X[,,i])
    W<-cbind( diag(m*r), kron(Ir,X[,,i-1]), kron(Ir,Y[,,i-1]%*%X[,,i-1]) )
    l<-l+t(W)%*%iV%*%z
    Q<-Q+t(W)%*%iV%*%W
  }

  iQ<-solve(Q)

  beta<-c( rmvnorm(1,iQ%*%l,iQ) )

  T<-matrix(beta[1:(m*r)],m,r)
  A<-t( matrix( beta[m*r + 1:r^2] ,r,r) )
  C<-t( matrix( beta[m*r + r^2 + 1:r^2] ,r,r) ) 

  list(T=T,A=A,C=C)
}

rS_fc<-function(Y,X,T,A,C)
{
  SS<-matrix(0,dim(X)[2],dim(X)[2]) 
  for(t in 2:dim(X)[3]) 
  {
    SS<-SS+ crossprod( X[,,t] - T - X[,,t-1]%*%t(A) - Y[,,t-1]%*%X[,,t-1]%*%t(C) ) 
  }
  solve(rwish( solve( SS+diag(nrow=dim(X)[2]) ), dim(X)[1]*(dim(X)[3]-1) + dim(X)[2] + 2 ) )
}


rX_fc<-function(Y,M,alpha,H,s2,X,T,A,C,S=diag(nrow=ncol(X)))
{

  m<-dim(Y)[1] ; n<-dim(Y)[3] ; r<-dim(X)[2]
  iS<-solve(S) 
  iV<-kron(iS,diag(m)) 

  for(t in 1:n)
  {
    for(i in 1:m)
    {

      ## -- contribution from past  
      if(t==1){ l0<-rep(0,r) ; Q0<-diag(nrow=r) }
      if(t>1)
      {
        l0<-iS%*%c( T[i,] + A%*%X[i,,t-1] + C%*%t(X[,,t-1])%*%Y[i,,t-1] )
        Q0<-iS
      }

      ## -- contribution from future
      if(t<n)
      {
        ## from attribute model
        z<-c( t( X[,,t+1] - T -  diag(m)[,-i]%*%X[-i,,t]%*%t(A) - Y[,-i,t]%*%X[-i,,t]%*%t(C) ) )
        W<-kron( diag(m)[i,], A) + kron( Y[i,,t] , C) 
        l1<-t(W)%*%iV%*%z
        Q1<-t(W)%*%iV%*%W 

        ## from network model
        z<- Y[i,-i,t+1] - M[i,-i] - alpha*Y[i,-i,t]
        W<- X[-i,,t]%*%H
        l2<-t(W)%*%z/s2
        Q2<-crossprod(W)/s2
      }

      l<-l0+l1+l2 ; Q<-Q0+Q1+Q2

      iQ<-solve(Q)
      X[i,,t]<-rmvnorm(1,iQ%*%l,iQ)
    }
  }
  X
}


#### ----

vech <- function (A, diag=FALSE){ return(A[upper.tri(A,diag=F)]) }



rbetaY_fc <- function (Y,U,s2,H_diag=TRUE,t2=1e2){ 


  N <- dim(Y)[3]
  R <- dim(U)[2]
  M <- dim(Y)[1]
  npar<-.5*M*(M-1) + 1 + R^(2-H_diag)
  mean_prior<-rep(0,npar)
  var_prior<-t2*diag(npar)

  var_prior_inv <- solve(var_prior)

  # get vectorization of Y
  y <- apply(Y,3,vech)
  U_outer_diag <- array(apply(U,2,apply,2,tcrossprod),dim=c(M,M,N,R))

  # initialize quadratic and linear term
  Q <- L <- 0

  # initialize design matrix
  if (H_diag){
    D <- matrix(nrow=(M-1)*M/2,ncol=(M-1)*M/2+R+1)
  }
  if (!H_diag){
    D <- matrix(nrow=(M-1)*M/2,ncol=(M-1)*M/2+1+R*(R+1)/2)
  }

  # design matrix for intercept
  D[,1:((M-1)*M/2)] <- diag((M-1)*M/2)

  for (i in 1:(N-1)){
    # diagonal of H
    if(R>1) {UD <- apply(U_outer_diag[,,i,],3,vech) }
    if(R==1){UD <- vech(U_outer_diag[,,i,] ) }

    # off-diagonal of H
    if (!H_diag){
      for (r in 1:(R-1)){
        for (s in (r+1):R){
          UD <- cbind(UD, vech(outer(c(U[,r,i]),c(U[,s,i]))) + vech(outer(c(U[,s,i]),c(U[,r,i]))))
                          }
                         }
                 }

    # quadratic and linear term
    if (H_diag){
      QD <- diag((M-1)*M/2+R+1)
      QD[1:((M-1)*M/2),((M-1)*M/2)+1] <- QD[((M-1)*M/2)+1,1:((M-1)*M/2)] <- y[,i]
      QD[((M-1)*M/2)+1,((M-1)*M/2)+1] <- sum(y[,i]^2)
      QD[1:((M-1)*M/2),((M-1)*M/2)+1+1:R] <- UD
      QD[((M-1)*M/2)+1+1:R,1:((M-1)*M/2)] <- t(UD)
      QD[((M-1)*M/2)+1,((M-1)*M/2)+1+1:R] <- QD[((M-1)*M/2)+1+1:R,((M-1)*M/2)+1] <- crossprod(UD,y[,i])
      QD[((M-1)*M/2)+1+1:R,((M-1)*M/2)+1+1:R] <- crossprod(UD)
    }

    if (!H_diag){
      QD <- diag((M-1)*M/2+1+R*(R+1)/2)
      QD[1:((M-1)*M/2),((M-1)*M/2)+1] <- QD[((M-1)*M/2)+1,1:((M-1)*M/2)] <- y[,i]
      QD[((M-1)*M/2)+1,((M-1)*M/2)+1] <- sum(y[,i]^2)
      QD[1:((M-1)*M/2),((M-1)*M/2)+1+1:(R*(R+1)/2)] <- UD
      QD[((M-1)*M/2)+1+1:(R*(R+1)/2),1:((M-1)*M/2)] <- t(UD)
      QD[((M-1)*M/2)+1,((M-1)*M/2)+1+1:(R*(R+1)/2)] <- QD[((M-1)*M/2)+1+1:(R*(R+1)/2),((M-1)*M/2)+1] <- crossprod(UD,y[,i])
      QD[((M-1)*M/2)+1+1:(R*(R+1)/2),((M-1)*M/2)+1+1:(R*(R+1)/2)] <- crossprod(UD)
    }

    # update quadratic and linear term
    Q <- Q + QD
    L <- L + c(y[,i+1],sum(y[,i]*y[,i+1]),crossprod(UD,y[,i+1]))
  }

  Q <- Q/s2 + var_prior_inv
  L <- L/s2 + var_prior_inv %*% mean_prior

  # get full conditional parameter
  var_pos <- chol2inv(chol(Q))
  mean_pos <- crossprod(var_pos, L)
  beta_Y<-rmvnorm(1, mean_pos, var_pos)

  # preparing output
  Mu <- matrix(0, ncol=M, nrow=M)
  H <- matrix(0, ncol=R, nrow=R)

  Mu[upper.tri(Mu,diag=F)] <- beta_Y[1:((M-1)*M/2)]
  Mu <- Mu+t(Mu)
  diag(Mu) <- NA

  alpha <- beta_Y[((M-1)*M/2)+1]

  H[lower.tri(H,diag=F)] <- beta_Y[((M-1)*M/2)+1+R+(1:(R*(R-1)/2))]
  H <- H+t(H)
  diag(H) <- beta_Y[((M-1)*M/2)+1+1:R]
  H[is.na(H)]<-0

  return(list(Mu=Mu, alpha=alpha, H=H))
}




#### ----
#rs2_fc <- function(Y,Mu,alpha,H,U,gamma_prior=c(2,2)){
rs2_fc <- function(Y,Mu,alpha,H,U,gamma_prior=c(1,1)/2){
  M <- dim(Y)[1]
  N <- dim(Y)[3]
  # calculate the multiplicative term
  Uprod <- function(U,H){
    return(tcrossprod(tcrossprod(U,H),U))
  }

  UHUT <- array(apply(U[,,-N,drop=F],3,Uprod,H),dim=c(M,M,N-1))

  # calculate the difference
  dY <- Y[,,2:N] - alpha*Y[,,2:N-1] - UHUT
  dY <- sweep(dY, 1:2, Mu, "-")

  dY_vech <- apply(dY,3,vech)

  s2 <- rgamma(1, (N-1)*M*(M-1)/4+gamma_prior[1], sum(dY_vech^2)/2+gamma_prior[2])
  s2 <- 1/s2

  return(s2)
}





