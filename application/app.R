raw_data <- read.csv("data_example.csv")
########### Dummy ###############
Cov = raw_data[,c(104:109)]

race_1 = ifelse(Cov$race == 1, 1, 0)
race_2 = ifelse(Cov$race == 2, 1, 0)
race_3 = ifelse(Cov$race == 3, 1, 0)
age = Cov$age
t1bmi = Cov$t1bmi
sg = Cov$sg
insurance_0 = ifelse(Cov$insurance == 0, 1, 0)
insurance_1 = ifelse(Cov$insurance == 1, 1, 0)
education_0 = ifelse(Cov$education == 0, 1, 0)
education_1 = ifelse(Cov$education == 1, 1, 0)
education_2 = ifelse(Cov$education == 2, 1, 0)
education_3= ifelse(Cov$education == 3, 1, 0)
Confounding = data.frame(race_1, race_2, race_3,
                      age, t1bmi, sg,
                      insurance_0, insurance_1,
                      education_0, education_1,
                      education_2, education_3)
Confounding <- as.matrix(Confounding)
Confounding <- as.matrix(cbind(t1bmi,sg))
#################################
ver1_data = apply(raw_data[,-104:-109],2,function(x)((x-mean(x))/sd(x)))
Exposure = ver1_data [,c(1:42)]
Mediator = ver1_data [,c(43:103)]
Y = ver1_data[,104]
#Cofounding = ver1_data [,c(104:109)]



for(i in c(1:38))
{
  cat("i =",i,'\n')
  Direct_MCP_dummy[[i]] = hima(Exposure[,i],Y,Mediator,COV.XM = Confounding,COV.MY = Confounding)
  #Direct_Pathway_dummy[[i]] = pathlasso.apply(raw_data,i,110,61,43,0.01)
  Direct_DAG_dummy[[i]] = INDAG(Exposure, Mediator,Y,i,Confound = Confounding)
}
M1 = hilma(Y = Y, G = Mediator,S=cbind(Exposure))


for(i in c(39:42))
{
  cat("i =",i,'\n')
  Direct_MCP_dummy[[i]] = hima(Exposure[,i],Y,Mediator,COV.XM = Confounding,COV.MY = Confounding)
  Direct_Pathway_dummy[[i]] = pathlasso.apply(raw_data,i,110,61,43,0.01)
  Direct_DAG_dummy[[i]] = INDAG(Exposure[,i], Mediator,Y,1,Confound = Confounding)
}
M1 = hilma(Y = Y, G = Mediator,S=cbind(Exposure))
M2 = hilma(Y = Y, G = Mediator,S=cbind(Exposure[,39]))
M3 = hilma(Y = Y, G = Mediator,S=cbind(Exposure[,40]))
M4 = hilma(Y = Y, G = Mediator,S=cbind(Exposure[,41]))
M5 = hilma(Y = Y, G = Mediator,S=cbind(Exposure[,42]))

Global_Table = data.frame(matrix(0,42,6))
Global_Table[,1] = colnames(raw_data[,1:42])
Global_Table[,2] = c(M1$beta_hat[1:38],M2$beta_hat[1],M3$beta_hat[1],M4$beta_hat[1],M5$beta_hat[1])

for(i in c(1:38))
{
  if(M1$pvalue_beta_hat[i] <= 0.05)
    Global_Table[i,3] = 1
  Global_Table[i,4] = M1$lower[i]
  Global_Table[i,5] = M1$upper[i]
}

Global_Table[39,4] = M2$lower
Global_Table[39,5] = M2$upper
Global_Table[40,4] = M3$lower
Global_Table[40,5] = M3$upper
Global_Table[41,4] = M4$lower
Global_Table[41,5] = M4$upper
Global_Table[42,4] = M5$lower
Global_Table[42,5] = M5$upper


colnames(Global_Table) <- c("Name","Value", "Siginificant", "CI.lower", "CI.upper")
Global_Table

##################

write.csv(Global_Table, file = "Global_table.csv")

###############
saveRDS(object = Direct_DAG_dummy, file = "Direct_DAG_dummy.rds")
saveRDS(object = Direct_MCP_dummy, file = "Direct_MCP_dummy.rds")
saveRDS(object = Direct_Pathway_dummy, file = "Direct_Pathway_dummy.rds")
saveRDS(object = Direct_Debias, file = "Direct_Debias.rds")


###################
library(ggplot2)
df <- Global_Table[1:38,]
df$Siginificant <- factor(df$Siginificant)
df <- df[which((df$CI.upper-df$CI.lower)<150),]
ggplot(df, aes(y = Name)) + 
  geom_errorbar(aes(xmin=CI.lower, xmax=CI.upper,col=Siginificant)) + 
  geom_point(aes(x=Value,col=factor(Siginificant))) +
  xlab("Global Indirect Effect") +
  ylab("Exposures") +
  guides(guide_legend("Significant test"))
###################

oh_phe_23 <- which(colnames(Exposure) == "OH_PHE_23")

Direct_MCP_dummy[[oh_phe_23]]$results$Bonferroni.p
Direct_MCP_dummy[[oh_phe_23]]$ID
mediator_oh_phe_23 <- list()
mediator_oh_phe_23

MCP_mediator = list()
DAG_Mediator = list()
for(j in c(1:38))
{
  if(length(which(Direct_MCP_dummy[[j]]$results$Bonferroni.p <= 0.05)) >= 1)
    MCP_mediator[[j]] = row.names(Direct_MCP_dummy[[j]]$results)[which(Direct_MCP_dummy[[j]]$results$Bonferroni.p <= 0.05)]
  else
    MCP_mediator[[j]] = 0
  
  if(length(which(Direct_DAG_dummy[[j]]$p.value <= 0.05)) >= 1)
    DAG_Mediator[[j]] = colnames(raw_data[,1:38])[which(Direct_DAG_dummy[[j]]$p.value <= 0.05)]
  else
    DAG_Mediator[[j]] = 0
}

colnames(raw_data)



colnames(Exposure)[37]
colnames(Mediator)[33]
