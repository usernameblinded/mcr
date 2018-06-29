#### ---- data, functions and libraries
load("~/Dropbox/Data/ICEWS_DV/Y_247x247x112") 

for(t in 1:dim(Y)[3]) { diag(Y[,,t])<-NA } 


#### ---- top m actors
m<-50
deg<-apply(abs(Y),1,sum,na.rm=TRUE) + apply(abs(Y),2,sum,na.rm=TRUE)
bigdeg<-which(deg>=sort(deg,decreasing=TRUE)[m] )
Y<-Y[bigdeg,bigdeg,]

#### ---- normal scores and symmetrize
Z<-rfunctions::zscores(Y) 
Y<-(Z+aperm(Z,c(2,1,3)))/2
Y<-Y-median(Y,na.rm=TRUE) 


