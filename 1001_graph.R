library(foreach)
library(ncvreg)
library(pcalg)

raw_data <- read.csv("data_example.csv")
ver1_data = apply(raw_data,2,function(x)((x-mean(x))/sd(x)))
Exposure = ver1_data [,c(1:42)]
Mediator = ver1_data [,c(43:103)]
Y = ver1_data[,104]
Confounding = ver1_data[,c(105:107)]

# ### Cyclooxygenase Pathway
# c(2,17,18,19,20,21,25,26,27,24,1,3,35,22)
# 
# ### Cytochrome p450 Pathway
# c(30,31,33,37,38,39,40,41,42,16,7,8,9,10,11,32,36,48)
# 
# ### Lipoxygenase Pathway
# c(28, 29, 6, 43, 46, 47, 34, 49, 50, 51, 4, 5, 12, 45)
# 
# ### Parent Compound
# c(15,14,13,52,53)
# 
# ### Oxidative Stress
# c(94,95)
# 
# ### Protein Damage
# c(96,97,98)
# 
# ### Inflammatory
# c(99:104)


Pathway_list = list()
Pathway_list[[1]] = c(2,17,18,19,20,21,25,26,27,24,1,3,35,22)
Pathway_list[[2]] = c(30,31,33,37,38,39,40,41,42,16,7,8,9,10,11,32,36,48)
Pathway_list[[3]] = c(28, 29, 6, 43, 46, 47, 34, 49, 50, 51, 4, 5, 12, 45)
Pathway_list[[4]] = c(15,14,13,52,53)
Pathway_list[[5]] = c(94,95)
Pathway_list[[6]] = c(96,97,98)
Pathway_list[[7]] = c(99:103)
for(i in c(1:4))
  Pathway_list[[i]][which(Pathway_list[[i]]>23)] = Pathway_list[[i]][which(Pathway_list[[i]]>23)] -1

for(i in c(1:4))
  Pathway_list[[i]][which(Pathway_list[[i]]>43)] = Pathway_list[[i]][which(Pathway_list[[i]]>43)] -1

for(i in c(5:7))
  Pathway_list[[i]] = Pathway_list[[i]] -42


Mediation_list1 = list()
Mediation_list1[[1]] = list()
Mediation_list1[[2]] = list()
Mediation_list1[[3]] = list()
Mediation_list1[[4]] = list()
for(i in c(1:4))
{for(j in c(1:7))
{
  
  if(length(which(DAG_list[[i]][[j]]$p.value <= 0.05)) >= 1)
    Mediation_list1[[i]][[j]] = colnames(raw_data[,(Pathway_list[[j]]+42)])[which(DAG_list[[i]][[j]]$p.value <= 0.05)]
  else
    Mediation_list1[[i]][[j]] = 0
  
}
}

########################################################################
P_list = list()
for(j in c(1:7))
{
  cat("i=",i,"j=",j,'\n')
  
  P_list[[j]] = hilma(Y,Mediator[,Pathway_list[[j]]],as.matrix(cbind(Exposure[,c(39:42)], Confounding)))
}

record <- matrix(nrow=4, ncol=7)
for(i in 1: 7){
  record[, i] <- P_list[[i]]$pvalue_beta_hat[1: 4]
}
colnames(record) = c("Cyclooxygenase", "Cytochrome", "Lipoxygenase", "Parent", "Oxidative Stress",
                     "Protein", "Inflammatory")

record
df <- matrix(nrow = 28, ncol = 3)
t = 1
for(i in 1: 7){
  for(j in 1: 4){
    print(t)
    df[t, 1] <- colnames(record)[i]
    df[t, 2] <- j
    df[t, 3] <- record[j, i]
    t = t + 1
  }
}
df <- data.frame(df)
colnames(df) <- c("Med", "Exp_gp", "p.value")
df$Exp_gp <- ifelse(df$Exp_gp == 1, "Phthalates", ifelse(df$Exp_gp == 2, "Phenols", ifelse(df$Exp_gp == 3, "Polycyclic", "Metals")))
df$Med <- as.factor(df$Med)
df$p.value <- as.numeric(df$p.value)
library(ggplot2)
ggplot(df, aes(y = Exp_gp, color = Med)) + 
  geom_point(aes(x = p.value)) + 
  scale_x_continuous(limits = c(0, 1)) + geom_vline(xintercept=0.05, linetype="dashed", col = "red") + 
  ylab("Exposure group") + xlab("p-value")

###########################################################################################
P_list2 = list()
for(j in c(1:7))
{
  cat("i=",i,"j=",j,'\n')
  
  P_list2[[j]] = hilma(Y,Mediator[,Pathway_list[[j]]],as.matrix(cbind(Exposure[,c(1:38)], Confounding)))
}
P_list2
record2 <- matrix(nrow=38, ncol=7)
for(i in 1: 7){
  record2[, i] <- P_list2[[i]]$pvalue_beta_hat[1: 38]
}
colnames(record2) = c("Cyclooxygenase", "Cytochrome", "Lipoxygenase", "Parent", "Oxidative Stress",
                      "Protein", "Inflammatory")

row.names(record2) <- colnames(Exposure)[1: 38]
record2
df2 <- matrix(nrow = 7*38, ncol = 3)
t = 1
for(i in 1: 7){
  for(j in 1: 38){
    print(t)
    df2[t, 1] <- colnames(record2)[i]
    df2[t, 2] <- row.names(record2)[j]
    df2[t, 3] <- record2[j, i]
    t = t + 1
  }
}
df2 <- data.frame(df2)
colnames(df2) <- c("Med", "Exp_gp", "p.value")
df2
df2$Med <- as.factor(df2$Med)
df2$p.value <- as.numeric(df2$p.value)
library(ggplot2)
ggplot(df2, aes(y = Exp_gp, color = Med)) + 
  geom_point(aes(x = p.value)) + 
  scale_x_continuous(limits = c(0, 1)) + geom_vline(xintercept=0.05, linetype="dashed", col = "red") + 
  ylab("Exposure") + xlab("p-value")
