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
  parent$pnode = list()
  

  for (i in 1:n.dags)
  {
    wgt.unique <- t(matrix(ad[i, ], p, p))
    if(length(which(wgt.unique[x.pos, ] != 0))==0)
    {
      parent$pnode[[i]] = 0
    }else{
      parent$pnode[[i]] <- c(which(wgt.unique[x.pos, ] != 0))
    }

  }
 
  return(parent)

}

################## S: the matrix of exposure
################## G: the matrix of all the mediator
################## Y: the vector of outcome
################## s: the exposure node you may concern
################## g: the mediator node you may concern

INDAG <- function(S,G,Y,s,g,Confound=NA){
  if(is.na(Confound))
  {
    p = dim(S)[2]
    q = dim(G)[2]
    n = dim(S)[1]
    
    s.G = apply(G, 2, function(x)x=x-mean(x))
    s.S = S - mean(S)
    s.Y = Y - mean(Y)
    
    Part1 = lm(G[,g]~S[,s])
    beta_part1 = coef(Part1)[2]
 
    R = s.G[,g] - s.S[,g]*coef(Part1)[2]
    Z_part1 = s.S[,g]*R
   
    suffStat <- list(C = cor(G), n = nrow(G))
    CP =  pc(suffStat, indepTest = gaussCItest,
                p = ncol(G), alpha = 0.01)
    Pr = ParentDAG(g,CP@graph)
    
    beta_part2 = 0
    Z_part2 = 0
    for(i in c(1:Pr$ndag))
    {
      if(Pr$pnode[[i]][1]!=0)
      {
        M = lm(Y~G[,g]+G[,c(Pr$pnode[[i]])])
        beta_part2 = beta_part2 + coef(M)[2]
        
        R = Y - s.G[,c(g,Pr$pnode[[i]])]%*%M$coefficients[-1]
        Z_part2 = Z_part2 + ( t(solve(cov(s.G[,c(g,Pr$pnode[[i]])]))%*%t(s.G[,c(g,Pr$pnode[[i]])]))[,1]*R )
          
      }else{
        
        M = lm(Y~G[,g])
        beta_part2 = beta_part2 + coef(M)[2]
        
        R = Y - s.G[,g]*M$coefficients[-1]
        Z_part2 = Z_part2 + (s.G[,g]*R)
        
      }
    }
    
    beta_part2 = beta_part2/Pr$ndag
    Z_part2 = as.vector(Z_part2/Pr$ndag)
    point_estimate = as.numeric(beta_part1*beta_part2)
    
    standard = sqrt(mean( (beta_part2*Z_part2 + beta_part1*Z_part1)^2))
    
    Quantity = sqrt(n)*point_estimate/standard
    p.value = 2*(1-pnorm(abs(Quantity)))
    
    IND = list()
    IND$indeffect = point_estimate
    IND$Statistics = Quantity
    IND$p.value = p.value
    IND$upper = point_estimate + 1.96*standard/sqrt(n)
    IND$lower = point_estimate - 1.96*standard/sqrt(n)
    
    return(IND)
    
  }else{
    beta_part1 = coef(lm(G~S+Confound))[2]
    suffStat <- list(C = cor(G), n = nrow(G))
    CP =  pc(suffStat, indepTest = gaussCItest,
             p = ncol(G), alpha = 0.01)
    
    Pr = ParentDAG(g,CP@graph)
    
    beta_part2 = 0
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
