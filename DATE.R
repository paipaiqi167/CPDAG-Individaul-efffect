######### S := Exposure vector 
######### Y := Outcome vector
######### thershold := the thershold dividing exposure into two part



DATE <- function( S,Y,nfold = 2,thershold=NA,Confounding = NA)
{
  if(is.vector(S))
  {
    n = length(S)
    
    if(is.na(Confounding))
    {
      if(is.na(thershold))
      {     
        thershold = mean(S)
        
        estimate.p1 = rep(0,nfold)
        estimate.p0 = rep(0,nfold)
        for(i in c(1:nfold))
        {
          set = sample(n,floor(n/nfold),replace = F)
          
          t.S = S[set]
          v.S = S[-set]
          t.Y = Y[set]
          v.Y = Y[-set]
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
          
        }
        
        total = estimate.p1-estimate.p0
        return(mean(total))
      }
    }
  }
  
  
}

DATE(Exposure[,1],Y)
lm(Y~Exposure)
