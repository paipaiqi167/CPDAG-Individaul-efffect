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

################## S: the vector of exposure
################## G: the matrix of all the mediator
################## Y: the vector of outcome
################## g: the mediator node you may concern


INDAG <- function(S,G,Y,g,Confound=NA){
  if(is.na(Confound))
  {
    beta_part1 = coef(lm(G~S))[2]
    suffStat <- list(C = cor(G), n = nrow(G))
    CP =  pc(suffStat, indepTest = gaussCItest,
                p = ncol(G), alpha = 0.01)
    
    Pr = ParentDAG(g,CP@graph)
    
    beta_part2 = 0
    for(i in c(1:Pr$ndag))
    {
      if(Pr$pnode[[i]][1]!=0)
        M = lm(Y~G[,g]+G[,c(Pr$pnode[[i]])])
      
      beta_part2 = 0
      beta_part2 = beta_part2 + coef(M)[2]
      
    }
    beta_part2 = beta_part2/Pr$ndag
    
    IND = list()
    IND$indeffect = as.numeric(beta_part1*beta_part2)
    
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
