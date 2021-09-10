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

## Simulation 
##################################
times = 100
Result_Ind_PathwayLasso1 = matrix(0,times,11)
Result_Ind_PathwayLasso1 = data.frame(Result_Ind_PathwayLasso1)
colnames(Result_Ind_PathwayLasso1) = c('trail','M1','M2','M3','M4','M5','M6','M7','M8','M9','M10')

Result_Ind_CPDAG1 = matrix(0,times,31)
Result_Ind_CPDAG1 = data.frame(Result_Ind_CPDAG1)
colnames(Result_Ind_CPDAG1) =c('trail','M1','M2','M3','M4','M5','M6','M7','M8','M9','M10','M1upper','M2upper'
                               ,'M3upper','M4upper','M5upper','M6upper','M7upper','M8upper','M9upper','M10upper'
                               ,'M1lower','M2lower','M3lower','M4lower','M5lower','M6lower','M7lower','M8lower','M9lower','M10lower')

Result_Ind_Situation1_true = matrix(0,times,11)
Result_Ind_Situation1_true= data.frame(Result_Ind_Situation1_true)
colnames(Result_Ind_Situation1_true) = c('trail','M1','M2','M3','M4','M5','M6','M7','M8','M9','M10')
  
  
  
scaling = c(0.8)
for( scales in  c(1))
{
  for(t in c(1:times))
  {
    set.seed(t+1000)
    cat('ttttttttttttttttttttttttt=',t+(scales-1)*100,'\n')
    n=161
    p=38
    q=10
    
    Cor1 = pcor(raw_data[,c(1:9)])[[1]]
    Cor2 = pcor(raw_data[,c(10:18)])[[1]]
    Cor3 = pcor(raw_data[,c(19:26)])[[1]]
    Cor4 = pcor(raw_data[,c(27:38)])[[1]]
    
    Exposure = matrix(0,n,38)
    for(i in c(1:38))
    {
      Exposure[,i] = rnorm(n,0,1)
    }
    Exposure[,c(1:9)] = Exposure[,c(1:9)]%*%Cor1
    Exposure[,c(10:18)] = Exposure[,c(10:18)]%*%Cor2
    Exposure[,c(19:26)] = Exposure[,c(19:26)]%*%Cor3
    Exposure[,c(27:38)] = Exposure[,c(27:38)]%*%Cor4
    
    for(i in c(1:38))
    {
      Exposure[,i] = (Exposure[,i] - mean(Exposure[,i]))/sd(Exposure[,i])
    }
    
    exp_coef = matrix(rnorm(p*q,0,2),p,q)
    b = rnorm(q,0,2)
    for(i in c(1:q))
    {
      if(runif(1,0,1) >= scaling[scales])
      {
        exp_coef[,i] = 0
        b[i] = 0
      }
    }
    b = b*sqrt(0.5/scaling[scales])
    for(i in c(1:p))
    {
      exp_coef[i,] = exp_coef[i,]*sqrt(0.5/scaling[scales])
    }
    
    Mediator = Exposure%*%exp_coef
    
    mp = 0.5
    M_matrix = matrix(0,q,q)
    for(i in c(1:q))
    {for(j in c(1:q))
    {
      if(i >= j)
      {
        M_matrix[i,j] = 0
      }
      else if(M_matrix[i,j] <= 0)
      {
        if(i <= 4 && j <= 5)
        {
          if(runif(1,0,1) <= mp)
            M_matrix[i,j] = rnorm(1,0,1)
        }else if(i >= 6 && i <= 9 && j <= 10)
        {
          if(runif(1,0,1) <= mp)
            M_matrix[i,j] = rnorm(1,0,1)
  
        }
      }
      
    }
    }
    
    Mediator[,1] <- Mediator[,1] + rnorm(n,0,1)
    for(i in c(2:q))
    {
      if(i > 2)
        Mediator[,i] <- Mediator[,i] + Mediator[,c(1:(i-1))]%*%M_matrix[c(1:(i-1)),i] + rnorm(n,0,1)
      else  ### odd number
        Mediator[,i] <- Mediator[,i] + Mediator[,c(1:(i-1))]*M_matrix[c(1:(i-1)),i] + rnorm(n,0,1)
    }
    
    for(i in c(1:q))
    {
      Mediator[,i] = (Mediator[,i] - mean(Mediator[,i]))/sd(Mediator[,i])
    }
    
    
    
    Y <- cbind(rep(1,n),Exposure)
    Y <- Y%*%rep(1,(p+1))+Mediator%*%b
    for(i in c(1:n))
    {
      Y[i] = Y[i] + rnorm(1,0,1)
    }
    
    X = cbind(Exposure,Mediator,Y)
    
    Path = pathlasso.apply(X,1,49,10,39,lambda=c(0.01))
    Result_Ind_PathwayLasso1[t,1] = t
    Result_Ind_PathwayLasso1[t,-1] = Path[[2]]
    

    DAG = INDAG(Exposure,Mediator,Y,1)
    
    Result_Ind_CPDAG1[t,1] = t
    Result_Ind_CPDAG1[t,c(2:(q+1))] = DAG$indeffect
    Result_Ind_CPDAG1[t,c((q+2):(2*q+1))] = DAG$upper
    Result_Ind_CPDAG1[t,c((2*q+2):(3*q+1))] = DAG$lower
    
    Result_Ind_Situation1_true[t,1] = t
    Result_Ind_Situation1_true[t,-1] = exp_coef[1,]*b
  }
}
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
