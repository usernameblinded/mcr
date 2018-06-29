source("functions.R")
# Network array
network<-array(dim=c(26,26,4))
for (t in 1:4){
  network[,,t]<-as.matrix(read.table(paste("klas12b/klas12b-net-",t,".dat",sep="")))
  diag(network[,,t])<-NA
}
network[network==9]<-NA

# Changing nodal attributes
delinquency<-as.matrix(read.table("klas12b/klas12b-delinquency.dat"))
alcohol<-cbind(NA,as.matrix(read.table("klas12b/klas12b-alcohol.dat")))
delinquency[delinquency==0]<-NA
alcohol[alcohol==0]<-NA

# Static dyadic
primary<-as.matrix(read.table("klas12b/klas12b-primary.dat"))
diag(primary)<-NA

# Static nodal attributes
demography<-as.matrix(read.table("klas12b/klas12b-demographics.dat"))
demography[demography==0]<-NA
demography<-demography[-21,]
sex<-1*(demography[,1]==1)
age<-demography[,2]
ethnicity<-1*(demography[,3]==1)
religion<-1*(demography[,4]!=2)

advice<-as.matrix(read.table("klas12b/klas12b-advice.dat"))
advice[advice==0]<-NA

# Student 21 nolonger in the class
network<-network[-21,-21,]
delinquency<-delinquency[-21,]
alcohol<-alcohol[-21,]
primary<-primary[-21,-21]
advice<-advice[-21,]

M<-dim(network)[1]
N<-dim(network)[3]
p<-1
Y<-network
W<-array(dim=c(M,p,N))
for (n in 1:N){
    W[,,n]<-cbind(delinquency[,n])
}
Nd<-3
Xdyad<-array(dim=c(M,M,Nd))
Xdyad[,,1]<-outer(sex,rep(1,M))
Xdyad[,,2]<-outer(rep(1,M),sex)
Xdyad[,,3]<-outer(sex,sex,"==")

main_mcmc(Y,W,Xdyad,niter=20000,burnin=20000,Zt=T,Xt=T,odens=25,thin=25,file.name="result.RData")

