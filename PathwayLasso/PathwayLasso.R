library("freebird")
library('tidyverse')
library('ggpubr')
   
#Data := A matrix composed by exposure matrix, mediator matrix and outcome vector.
#exposure := A scalar represent the column of data should be the exposure.
#outcome := A scalar represent the column of data should be the outcome.
#nmed := number of mediator.
#med_start := The number of column of the mediators start.
#lambda := the lambda one may want to use. It can also be a vector providing the lambda set to be choosen by cross-validation.


pathlasso.apply<-function(Data,exposure,outcome,nmed,med_start,lambda){
  A<-Data[,exposure]
  M<-Data[,((med_start):(med_start+nmed-1))]
  Y<-Data[,outcome] #replace with your outcome variable
  n=length(Y)
  q=nmed
  # approximate covariance matrix
  Sigma10<-diag(rep(1,q))
  Sigma20<-matrix(1,1,1)
  
  m.A<-mean(A)
  m.M<-apply(M,2,mean)
  m.Y<-mean(Y)
  sd.A<-sd(A)
  sd.M<-apply(M,2,sd)
  sd.Y<-sd(Y)
  
  A<-scale(A) 
  M<-scale(M) 
  Y<-scale(Y) 
  
  # Pathway Lasso method parameters
  phi<-2
  gamma<-0           # adpative lasso parameter
  rho<-1             # ADMM parameter
  max.itr<-10000     # change if needed for model convergence 
  tol<-1e-6
  thred<-1e-6
  thred2<-1e-3
  omega.p<-0.1      # omega = omega.p*lambda
  
  AB.est=B.est=A.est<-c()
  Global = c()
  lambda.criteria = rep(0,length(lambda))
  
  if(length(lambda)!=1)
  {
    for (i in 1:length(lambda)) {
      
      for(s in c(1:5))
      {
        cat('validation.lambda = ',lambda[i],'slice = ',s,'\n')
        train.A = A[-c(seq(s,n,5))]
        train.M = M[-c(seq(s,n,5)),]
        train.Y = Y[-c(seq(s,n,5))]
        
        valid.A = A[c(seq(s,n,5))]
        valid.M = M[c(seq(s,n,5)),]
        valid.Y = Y[c(seq(s,n,5))]
        
        
        Sigma10<-diag(rep(1,q))
        Sigma20<-matrix(1,1,1)
        
        n=length(train.Y)
        test <- mediation_net_ADMM_NC(train.A,train.M,train.Y,lambda=lambda[i],omega=omega.p*lambda[i],phi=phi,Phi1=diag(c(0,rep(1,q))),
                                   Phi2=diag(rep(1,q+1)),rho=rho,rho.increase=FALSE,tol=tol,max.itr=max.itr,thred=thred,
                                   Sigma1=Sigma10,Sigma2=Sigma20,trace=FALSE)
        
        lambda.criteria[i] = lambda.criteria[i] + test$BIC
        
      }
    }
    lambda.select = lambda[which(lambda.criteria==min(lambda.criteria))]
  }else{
    lambda.select = lambda
  }

  cat('lambda.select = ',lambda.select,'\n')
  
  out<-mediation_net_ADMM_NC(A,M,Y,lambda.select,omega=omega.p*lambda.select,phi=phi,Phi1=diag(c(0,rep(1,q))),
                             Phi2=diag(rep(1,q+1)),rho=rho,rho.increase=FALSE,tol=tol,max.itr=max.itr,thred=thred,
                             Sigma1=Sigma10,Sigma2=Sigma20,trace=FALSE)
  
  A.est<-out$A*(sd.M/sd.A) #exposure to mediator selection
  B.est<-out$B*(sd.Y/sd.M) #mediator to outcome selection
  AB.est<-A.est*t(B.est) #mediation effect selection
  Global<- A.est%*%B.est
  
  results<-list(A.est,B.est,AB.est,Global,out$C,lambda.select)
  return(results)
  
}
