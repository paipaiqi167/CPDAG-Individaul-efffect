######### S := Exposure vector (or matrix)
######### Y := Outcome vector
######### type := the type of the exposure & Y: "binary" or "continuous
######### nfold := n-fold cross-fitting: decideing by number n

DATE <- function(S, Y, nfold = 2, type = type, Confounding = NA){
  
  if(type == "binary"){
    #print("Binary estimation")
    if(is.vector(S)){
      n = length(S)
      random = sample(n,n,replace=F)
      S = S[random]
      Y = Y[random]
      
      if(is.na(Confounding)){
        estimate.p1 = rep(0,nfold)
        estimate.p0 = rep(0,nfold)
        for(i in c(1:nfold)){
          t.S = S[seq(i,n,nfold)]
          v.S = S[-seq(i,n,nfold)]
          t.Y = Y[seq(i,n,nfold)]
          v.Y = Y[-seq(i,n,nfold)]
          t.n = length(t.S)
          v.n = length(v.S)
          
          p.1 = sum(t.S) / t.n
          p.0 = 1 - p.1
          
          e.y1 = mean(t.Y[which(t.S == 1)])
          e.y0 = mean(t.Y[which(t.S == 0)])
          
          for(j in c(1:v.n)){
            estimate.p1[i] = estimate.p1[i] + (ifelse(v.S[j] == 1,1,0)/p.1)*(v.Y[j] - e.y1) + e.y1
            estimate.p0[i] = estimate.p0[i] + (ifelse(v.S[j] == 0,1,0)/p.0)*(v.Y[j] - e.y0) + e.y0
          }
          
          estimate.p1[i] = estimate.p1[i] /v.n
          estimate.p0[i] = estimate.p0[i] /v.n
        }
        
        total = estimate.p1 - estimate.p0
      }
    }
    # matrix (not done yet)
    
  }
  
  ############## continuous
  if(type == "continuous"){
    #print("Continuous estimation")
    if(is.vector(S)){
      n = length(S)
      random = sample(n,n,replace=F)
      S = S[random]
      Y = Y[random]
      
      if(is.na(Confounding)){
        estimate.p1 = rep(0,nfold)
        estimate.p0 = rep(0,nfold)
        for(i in c(1:nfold)){
          t.S = S[seq(i,n,nfold)]
          v.S = S[-seq(i,n,nfold)]
          t.Y = Y[seq(i,n,nfold)]
          v.Y = Y[-seq(i,n,nfold)]
          t.n = length(t.S)
          v.n = length(v.S)
          
          p.1 = sum(t.S) / t.n
          p.0 = 1 - p.1
          
          e.y1 = mean(t.Y)
          e.y0 = mean(t.Y - 1)
          
          for(j in c(1:v.n)){
            estimate.p1[i] = estimate.p1[i] + (v.S[j]/p.1)*(v.Y[j] - e.y1) + e.y1
            estimate.p0[i] = estimate.p0[i] + ((v.S[j] - 1)/p.0)*(v.Y[j] - e.y0) + e.y0
          }
          
          estimate.p1[i] = estimate.p1[i] /v.n
          estimate.p0[i] = estimate.p0[i] /v.n
        }
        
        total = estimate.p1 - estimate.p0
      }
    }
    # matrix (not done yet)
    
  }
  return(mean(total))
}

