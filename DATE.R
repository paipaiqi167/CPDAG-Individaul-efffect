######### S := Exposure vector or matrix
######### Y := Outcome vector
######### thershold := the thershold dividing exposure into two part

DATE <- function( S,Y,nfold = 2,thershold=NA,Confounding = NA)
{
  if(is.vector(S))
  {
    n = length(S)
    random = sample(n,n,replace=F)
    S = S[random]
    Y = Y[random]
    
    if(is.na(Confounding))
    {
      
      if(is.na(thershold))
      {     
        thershold = mean(S)
      }
      
      estimate.p1 = rep(0,nfold)
      estimate.p0 = rep(0,nfold)
      for(i in c(1:nfold))
      {
        
        t.S = S[seq(i,n,nfold)]
        v.S = S[-seq(i,n,nfold)]
        t.Y = Y[seq(i,n,nfold)]
        v.Y = Y[-seq(i,n,nfold)]
        t.n = length(t.S)
        v.n = length(v.S)
        
        p.1 = sum(ifelse(t.S >= thershold,1,0)) / t.n
        p.0 = 1- p.1
        
        e.y1 = mean(t.Y[which(t.S>= thershold)])
        e.y0 = mean(t.Y[which(t.S < thershold)])
        
        for(j in c(1:v.n))
        {
          estimate.p1[i] = estimate.p1[i] + (ifelse(v.S[j] >= thershold,1,0)/p.1)*(v.Y[j] - e.y1) + e.y1
          estimate.p0[i] = estimate.p0[i] + (ifelse(v.S[j] < thershold,1,0)/p.0)*(v.Y[j] - e.y0) + e.y0
          
        }
        estimate.p1[i] = estimate.p1[i] /v.n
        estimate.p0[i] = estimate.p0[i] /v.n
      }
      
      total = estimate.p1 -estimate.p0
      return(mean(total))
    
    }
  }
  
  if(is.matrix(S))
  {
    n = dim(S)[1]
    p = dim(S)[2]
    random = sample(n,n,replace=F)
    S = S[random,]
    Y = Y[random]
    
    if(is.na(Confounding))
    {
      
      if(is.na(thershold))
      {     
        thershold = c()
        thershold = apply(S,2,function(x)mean(x))
      }
      
      out = list()
      out$total = c()
      for(k in c(1:p))
      {
          estimate.p1 = rep(0,nfold)
          estimate.p0 = rep(0,nfold)
          for(i in c(1:nfold))
          {
            
            t.S = S[seq(i,n,nfold),k]
            v.S = S[-seq(i,n,nfold),k]
            t.Y = Y[seq(i,n,nfold)]
            v.Y = Y[-seq(i,n,nfold)]
            t.n = length(t.S)
            v.n = length(v.S)
            
            p.1 = sum(ifelse(t.S >= thershold[k],1,0)) / t.n
            p.0 = 1- p.1
            
            e.y1 = mean(t.Y[which(t.S>= thershold[k])])
            e.y0 = mean(t.Y[which(t.S < thershold[k])])
            
            for(j in c(1:v.n))
            {
              estimate.p1[i] = estimate.p1[i] + (ifelse(v.S[j] >= thershold[k],1,0)/p.1)*(v.Y[j] - e.y1) + e.y1
              estimate.p0[i] = estimate.p0[i] + (ifelse(v.S[j] < thershold[k],1,0)/p.0)*(v.Y[j] - e.y0) + e.y0
              
            }
            estimate.p1[i] = estimate.p1[i]/v.n
            estimate.p0[i] = estimate.p0[i]/v.n
          }
          out$total[k] = mean(estimate.p1 -estimate.p0)
      }
      
      return(out)
      
    }
  }
  
}
