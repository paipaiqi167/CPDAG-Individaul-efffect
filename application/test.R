set.seed(t+1000)
n=161
p=1
q=60

Exposure = rnorm(n,0,1)

exp_coef <- runif(q,0.5,1)
exp_coef <- exp_coef*ifelse(runif(q,0,1) <= 0.5,1,-1)
exp_coef <- exp_coef*ifelse(runif(q,0,1)<= 0.1,1,0)

b <- runif(q,0.5,1)
b <- b*ifelse(runif(q,0,1) <= 0.5,1,-1)
b <- b*ifelse(runif(length(b),0,1) <= 0.1,1,0)
Mediator = Exposure%*%t(exp_coef)

dp = 5/choose(61,2) #0.00273224
M_matrix = M_generate(q,dp)
Mediator = Mediator%*%M_matrix + matrix(rnorm(n*q,0,1),n,q)

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
M1= lm(Y~Exposure)
mean(abs(M1$residuals))
