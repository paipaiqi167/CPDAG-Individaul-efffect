library(ppcor)

n=161
p=38
q=61

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

for(i in c(1:p))
{
  Exposure[,i] = (Exposure[,i] - mean(Exposure[,i]))/sd(Exposure[,i])
}

exp_coef = matrix(rnorm(p*q,0,1),p,q)
for(i in c(1:p))
  exp_coef[i,] = exp_coef[i,]*ifelse(runif(q,0,1)<= 0.2,1,0)

b <- rnorm(q,0,1)
b <- b*ifelse(runif(length(b),0,1) <= 0.2,1,0)
Mediator = Exposure%*%exp_coef

dp = 0.067
M_matrix = M_generate(q,dp)
Mediator = Mediator%*%M_matrix + matrix(rnorm(1,0,1),n,q)


Y <- cbind(rep(1,n),Exposure)
Y <- Y%*%rep(1,(p+1))+Mediator%*%b
for(i in c(1:n))
{
  Y[i] = Y[i] + rnorm(1,0,1)
}

X = cbind(Exposure,Mediator,Y)


M_generate <- function(q,p){
  M_matrix = diag(rep(1,q),)
  for(i in c(1:(q-1)))
  {
    for(j in c((i+1):q))
    {
      if(runif(1,0,1) <= p)
      {
        if(runif(1,0,1) <= 0.5)
        {
          M_matrix[i,j] = ifelse(runif(1,0,1) <= 0.5, 1,-1)*runif(1,0.5,2)
        }else{
          M_matrix[j,i] = ifelse(runif(1,0,1) <= 0.5, 1,-1)*runif(1,0.5,2)
        }
      }
    }
  }
  return(M_matrix)
}
