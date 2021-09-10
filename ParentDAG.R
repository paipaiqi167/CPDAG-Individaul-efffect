library(stringr)
library(pcalg)

################## x.pos: the mediator node you may concern its parent node

ParentDAG <- function(x.pos=NA,graphEst){
  
  p <- numNodes(graphEst)
  ad.g <- wgtMatrix(graphEst)
  ad <- pdag2allDags(ad.g)$dags
  n.dags <- nrow(ad)

  
  parent = list()
  parent$ndag = n.dags
  parent$node = list()
  
  for(j in c(1:length(x.pos)))
  {
    parent$node[[j]] = list()
  }
  
  for (i in 1:n.dags)
  {for(j in c(1:length(x.pos)))
  {
    wgt.unique <- t(matrix(ad[i, ], p, p))
    if(length(which(wgt.unique[x.pos[j], ] != 0))==0)
    {
      parent$node[[j]][[i]] = 0
    }else{
      parent$node[[j]][[i]] <- as.vector(which(wgt.unique[x.pos[j], ]!= 0))
    }
  }
  }
 
  return(parent)

}
################## S: the matrix of exposure
################## G: the matrix of all the mediator
################## Y: the vector of outcome
################## s: the exposure node you may concern

INDAG <- function(S,G,Y,s,Confound=NA){
  if(is.na(Confound))
  {
    
    p = dim(S)[2]
    q = dim(G)[2]
    n = dim(S)[1]
    
    s.G = apply(G, 2, function(x)x=x-mean(x))
    s.S = S - mean(S)
    s.Y = Y - mean(Y)
    
    beta_part1 = c()
    Z_part1 = matrix(0,n,q)
    
    for(j  in c(1:q))
    {
      Part1 = lm(G[,g[j]]~S[,s])
      beta_part1[j] = coef(Part1)[2]
   
      R = s.G[,g[j]] - s.S[,g[j]]*coef(Part1)[2]
      Z_part1[,j] = s.S[,g[j]]*R
    }
    
    
    suffStat <- list(C = cor(G), n = nrow(G))
    CP =  pc(suffStat, indepTest = gaussCItest,
                p = ncol(G), alpha = 0.01)
    Pr = ParentDAG(g,CP@graph)
    
    
    beta_part2 = rep(0,q)
    Z_part2 = matrix(0,n,q)
    
    for (j in c(1:q)) 
    {
      for(i in c(1:Pr$ndag))
      {
        if(Pr$node[[j]][[i]][1]!=0)
        {
          M = lm(Y~G[,g[j]]+G[,c(Pr$node[[j]][[i]])])
          beta_part2[j] = beta_part2[j]  + coef(M)[2]
          
          R = Y - s.G[,c(g[j],Pr$node[[j]][[i]])]%*%M$coefficients[-1]
          Z_part2[,j]  = Z_part2[,j]  + ( t(solve(cov(s.G[,c(g[j],Pr$node[[j]][[i]])]))%*%t(s.G[,c(g[j],Pr$node[[j]][[i]])]))[,1]*R )
            
        }else{
          
          M = lm(Y~G[,g[j]])
          beta_part2[j]  = beta_part2[j]  + coef(M)[2]
          
          R = Y - s.G[,g[j]]*M$coefficients[-1]
          Z_part2[,j]  = Z_part2[,j]  + (s.G[,g[j]]*R)
          
        }
      }
      
      beta_part2[j] = beta_part2[j] /Pr$ndag
      Z_part2[,j] = as.vector(Z_part2[,j]/Pr$ndag)
      
    }
    
    point_estimate = as.vector(beta_part1*beta_part2)
    
    standard = c()
    for(i in c(1:q))
    {
      standard[i] = sqrt(mean((beta_part2[i]*Z_part2[i] + beta_part1[i]*Z_part1[i])^2))
    }
    
    Quantity = sqrt(n)*point_estimate/standard
    p.value = 2*(1-pnorm(abs(Quantity)))
    
    upperbound = point_estimate + 1.96*standard/sqrt(n)
    lowerbound = point_estimate - 1.96*standard/sqrt(n)
    
    IND = list()
    IND$indeffect = point_estimate
    IND$Statistics = Quantity
    IND$p.value = p.value
    IND$upper = upperbound
    IND$lower = lowerbound
    
    return(IND)
    
  }else{
    
    beta_part1 = coef(lm(G~S+Confound))[2]
    suffStat <- list(C = cor(G), n = nrow(G))
    CP =  pc(suffStat, indepTest = gaussCItest,
             p = ncol(G), alpha = 0.01)
    
    Pr = ParentDAG(g,CP@graph)
    
    beta_part2 = 0
    for(j in c(1:q))
    {}
    for(i in c(1:Pr$ndag))
    {
      
      if(Pr$pnode[[i]][1]!=0)
        M = lm(Y~G[,g]+G[,c(Pr$pnode[[i]])]+Confound)
      
      beta_part2 = 0
      beta_part2 = beta_part2 + coef(M)[2]
      
    }
    beta_part2 = beta_part2/Pr$ndag
    
    IND = list()
    IND$indeffect = as.numeric(beta_part1*beta_part2)
   
    return(IND)
   
}
}
