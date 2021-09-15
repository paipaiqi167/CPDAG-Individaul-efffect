library(stringr)
library(partitions)
library('graph')

ParentDAG <- function(x.pos=NA,graphEst){
  
  p <- numNodes(graphEst)
  ad.g <- wgtMatrix(graphEst)
  q = dim(ad.g)[1]
  parent = list()
  parent$node = list()
  parent$ndag = rep(1,q)
  
  if(max(ad.g+t(ad.g)) <= 1)
  {
    for(j in c(1:length(x.pos)))
    {
      wgt.unique <- wgtMatrix(graphEst)
      if(length(which(wgt.unique[x.pos[j], ] != 0))==0)
      {
        parent$node[[j]] = 0
      }else{
        parent$node[[j]] <- as.vector(which(wgt.unique[x.pos[j], ]!= 0))
      }
    }
    
    cat("Caculate mannually, DAG = 1 \n")
    return(parent)
    
  }else{
    
    Skeletons = ad.g+t(ad.g)
    S.set = list()
    for(j in c(1:q))
    {
      S.set[[j]] = 0
    }
    
    M.set = list()
    v. = 0
    s = c()
    for(i in c(1:(q-1)))
    {
      for(j in c((i+1):q))
    {
      if(Skeletons[i,j]==2)
      {
        if(length(intersect(s,c(i,j)))==0)
        {
          S.set[[i]] = as.vector(as.matrix(c(S.set[[i]],j)))
          S.set[[j]] = as.vector(as.matrix(c(S.set[[j]],i)))
          
          s = c(s,i,j)
          
          S.set[[i]] = S.set[[i]][-which(S.set[[i]]==0)]
          S.set[[j]] = S.set[[j]][-which(S.set[[j]]==0)]
        }else{
          
          v. =1
          S.set[[i]] = as.vector(as.matrix(c(S.set[[i]],j)))
          S.set[[j]] = as.vector(as.matrix(c(S.set[[j]],i)))
          
          s = c(s,i,j)
          
          if(length(which(S.set[[i]]==0))>0)
            S.set[[i]] = S.set[[i]][-which(S.set[[i]]==0)]
          if(length(which(S.set[[j]]==0))>0)
            S.set[[j]] = S.set[[j]][-which(S.set[[j]]==0)]
          
        }
      }
    }
    }

    parent$node= list()
    for(j in c(1:q))
    {    
      parent$node[[j]] = list()
    }
    cluster.set = list()
    
    NewS.set = S.set
    count = 1
    for(i in c(1:q))
    {
      if( sum(NewS.set[[i]] )!= 0)
      {
        cluster.set[[count]] = list()
        un.set = unique(Searching(S.set,i))
        
        if(!is.na(un.set[1]))
        {
        cluster.set[[count]] = un.set
        for(j in c(1:length(cluster.set[[count]])))
        {
          NewS.set[[ cluster.set[[count]][j] ]] = 0
        }
        count = count+1
        }
        
      }
    }
    
    parent$ndag = rep(1,q)
    error = 0
    i=1
    wgt.unique =  ad.g
    while(i <=(count-1) && error == 0)
    {
        node.set = unlist(cluster.set[[i]])
        l = length(node.set)
        for(j in c(1:l))
        {
          if(length(as.vector(which(wgt.unique[node.set[j], ]!= 0))) >= l)
          {
            error = 1
            cat("unknown structure \n")
            return(NA)
          }
        }
        
        parent$ndag[node.set] = l
        for( j in c(1:l)) #### j is ndag
        {for(m in c(1:l)) #### m is node
        {
          
          if(j==l)
          {
            parent$node[[node.set[m]]][[j]] = 0
          }else{
            
            if(l == 2)
            {
              parent$node[[node.set[m]]][[j]] = node.set[-m]
              
            }else{
            if(j <= (m-1))
            {
              parent$node[[node.set[m]]][[j]] = node.set[(m-1)]
            }else{
              parent$node[[node.set[m]]][[j]] = node.set[(m+1)]
            }
            }
          }  
          
        } 
        }
        i=i+1
    }
    
    if(error ==0)
    {
      cat("almost done")
      wgt.unique <- wgtMatrix(graphEst)
      for(j in c(1:length(x.pos)))
      {
        if(parent$ndag[j] == 1)
        {
          if(length(which(wgt.unique[x.pos[j], ] != 0))==0)
          {
            parent$node[[j]][[1]] = 0
          }else{
            parent$node[[j]][[1]] <- as.vector(which(wgt.unique[x.pos[j], ]!= 0))
          }
        }
        
      }
      cat("Caculate mannually \n")
      return(parent)
    
    
    }
}

  ######################################################
  #cat("Caculate by package \n")
  
  #ad <- pdag2allDags(ad.g)$dags
  #n.dags <- nrow(ad)
  
  #parent$ndag = n.dags
  
  #for (i in 1:n.dags)
  #{for(j in c(1:length(x.pos)))
  #{
  #  wgt.unique <- t(matrix(ad[i, ], p, p))
  #  if(length(which(wgt.unique[x.pos[j], ] != 0))==0)
  #  {
  #    parent$node[[j]][[i]] = 0
  #  }else{
  #    parent$node[[j]][[i]] <- as.vector(which(wgt.unique[x.pos[j], ]!= 0))
  #  }
  #}
  #}
  #return(parent)
}

INDAG <- function(S,G,Y,s,Confound=NA){
  if(is.na(Confound))
  {
    
    p = dim(S)[2]
    q = dim(G)[2]
    n = dim(S)[1]
    
    g=c(1:q)
    
    s.G = apply(G, 2, function(x)x=x-mean(x))
    s.S = S - mean(S)
    s.Y = Y - mean(Y)
    
    beta_part1 = c()
    Z_part1 = matrix(0,n,q)
    
    for(j  in c(1:q))
    {
      Part1 = lm(G[,g[j]]~S[,s])
      beta_part1[j] = coef(Part1)[2]
   
      R = s.G[,g[j]] - s.S[,s]*coef(Part1)[2]
      Z_part1[,j] = s.S[,s]*R
    }
    
    
    suffStat <- list(C = cor(G), n = nrow(G))
    CP =  pc(suffStat, indepTest = gaussCItest,
                p = ncol(G), alpha = 0.01)
    cat("CP finished \n")
    Pr = ParentDAG(g,CP@graph)
    
    if(!is.list(Pr))
    {
      cat("unkown structure\n")
      return(0)
    }
    
    beta_part2 = rep(0,q)
    Z_part2 = matrix(0,n,q)
    
    for (j in c(1:q)) 
    {
      print(j)
      for(i in c(1:Pr$ndag[j]))
      {
        if(Pr$node[[j]][[i]][1]!=0)  ###　this node has parent nodes
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
      
      beta_part2[j] = beta_part2[j] /Pr$ndag[j]
      Z_part2[,j] = as.vector(Z_part2[,j]/Pr$ndag[j])
      
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

Searching <- function(Set,k,avoid=NA)
{
  n.k = as.vector(t(as.matrix(k)))
  n.set = unlist(Set[k][1])
  
  if(is.na(avoid[1])&&length(n.set) >= 2)   ### reject the node which is not on the head
    return(NA)
    
  if(!is.na(avoid)) 
  {
    if(sum(intersect(n.set,avoid)-n.set)==0)
    {
      return(NULL)
    }
  }
  out = c(n.k,n.set)
  
  if(is.na(avoid)){
    out = c(out,Searching(Set,n.set,n.k))
  }else{
    out = c(out,Searching(Set,n.set[-which(n.set==avoid)],n.k))
    
  }
  return(out)
}

