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
        if(sum(intersect(un.set,-1)) < 0)
        {
          un.set = un.set[-which(un.set==-1)]
          un.set = c(-1,un.set)
          type = 1
        }
        
        if(!is.na(un.set[1]))
        {
          cluster.set[[count]] = un.set
          for(j in c((1+type):(length(cluster.set[[count]]))))
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

    while(i <=(count-1) && error == 0)
    {
        node.set = unlist(cluster.set[[i]])
        l = length(node.set)
        wgt.unique <- wgtMatrix(graphEst)
        
        if(node.set[1]==-1){
        node.set = node.set[-1]
        l = l-1
        }
        
        target.set = node.set
        for(m  in c(1:l))
        {
          target.set = c(target.set,as.vector(which(wgt.unique[node.set[m], ]!= 0)))
          target.set = c(target.set,as.vector(which(wgt.unique[,node.set[m]]!= 0)))
        }
        target.set = unique(target.set)
        target.set =  sort(target.set) 
        target.matrix = wgt.unique[target.set,target.set]
        tp = length(target.set)
        
        Target <- pdag2allDags(target.matrix)$dags
        number <- nrow(Target)
        if(is.null(number))
        {
          parent$ndag[node.set[j]] = 0
          for(j in c(1:l))
          {
              parent$node[[node.set[j]]][[1]] = 0
          }
          i=i+1
          
        }else{
          
        parent$ndag[node.set] <- nrow(Target)
        
        for (m in 1:number)
        {for(j in c(1:l))
        {
          wgt.unique <- t(matrix(Target[m,],tp,tp))
          if(length(which(wgt.unique[order(target.set)[which(target.set==node.set[j])], ] != 0))==0)
          {
            parent$node[[node.set[j]]][[m]] = 0
          }else{
            parent$node[[node.set[j]]][[m]] <- target.set[as.vector(which(wgt.unique[order(target.set)[which(target.set==node.set[j])], ]!=0 ))]
          }
        }
        }
        i=i+1
        }
    }
    

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
      
      if(parent$ndag[j] == 0)
        parent$ndag[j] = 1
    }
    return(parent)
  
  

}

}

INDAG <- function(S,G,Y,s,Confound=NA,True.Parent = NA){
  if(is.na(Confound))
  {
    S = as.matrix(S)
    
    p = dim(S)[2]
    q = dim(G)[2]
    n = length(Y)
    
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
    
    if(!anyNA(True.Parent))
    {
      Pr = ParentDAG2(True.Parent)
    }else{
      
      suffStat <- list(C = cor(G), n = nrow(G))
      CP =  pc(suffStat, indepTest = gaussCItest,
                  p = ncol(G), alpha = 0.01)
      Pr = ParentDAG(g,CP@graph)
    }
    
    beta_part2 = rep(0,q)
    Z_part2 = matrix(0,n,q)
    
    for (j in c(1:q)) 
    {

      for(i in c(1:Pr$ndag[j]))
      {
        #print(c(j,i))
        if(Pr$node[[j]][[i]][1]!=0)  ###　this node has parent nodes
        {
          
          M = lm(Y~G[,g[j]]+G[,c(Pr$node[[j]][[i]])]+S[,s])
          beta_part2[j] = beta_part2[j]  + coef(M)[2]
          
          R = Y - s.G[,c(g[j],Pr$node[[j]][[i]])]%*%M$coefficients[-c(1,length(M$coefficients))]
          Z_part2[,j]  = Z_part2[,j]  + ( t(solve(cov(s.G[,c(g[j],Pr$node[[j]][[i]])]))%*%t(s.G[,c(g[j],Pr$node[[j]][[i]])]))[,1]*R )
            
        }else{
          
          M = lm(Y~G[,g[j]]+S[,s])
          beta_part2[j]  = beta_part2[j]  + coef(M)[2]
          
          R = Y - s.G[,g[j]]*M$coefficients[-c(1,length(M$coefficients))]
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
      standard[i] = sqrt(mean((beta_part2[i]*Z_part2[,i] + beta_part1[i]*Z_part1[,i])^2))
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
  n.set = setdiff(n.set,intersect(n.set,avoid))
  
  if(anyNA(avoid)&&length(n.set) >= 2)   ### reject the node which is not on the head
    return(NA)
   
  if(!anyNA(avoid)) 
  {
    if(sum(n.set)==0)
    {

      return(NULL)
    }
  }
  out = c(n.k,n.set)
  
  if(length(setdiff(n.set,intersect(n.set,avoid))) <=1)  ## whether is a core
  {
    if(anyNA(avoid)){
      
      avoid.set = c(n.k,n.set)
      out = c(out,Searching(Set,n.set,avoid.set))
      
    }else{
      avoid.set = c(n.k,n.set,avoid)
      out = c(out,Searching(Set,n.set,avoid.set))
      
    }
    return(out)
    
  }else{
    
    target.set = setdiff(n.set,intersect(n.set,avoid))
    if(anyNA(avoid)){avoid.set = c(n.k,n.set)}else{avoid.set = c(n.k,n.set,avoid)}
    m = length(target.set)
    for(i in c(1:m))
    {
      out = c(out,Searching(Set,target.set[i],avoid.set))
      avoid.set = unique(c(out,avoid.set))
      
    }
    out = unique(out)
    return(c(-1,out))
  }
}

ParentDAG2 <- function(M){
  
  q <- dim(M)[1]

  parent = list()
  parent$node = list()
  parent$ndag = rep(1,q)

  for(j in c(1:q))
  {
    if(length(which(M[j,-j] != 0))==0)
    {
      parent$node[[j]] = 0
    }else{
      parent$node[[j]] <- setdiff(as.vector(which(M[j,]!= 0)),j)
    }
  }
  
  return(parent)
    
}




