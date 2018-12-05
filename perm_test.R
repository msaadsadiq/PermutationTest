
library(ggplot2)
library(reshape)

numReps = 1000 
sd.pite <- data.frame(sd=numeric(numReps))
rmse.pite <- data.frame(mse=numeric(500))

##### Calculate SD

for (j in 1:numReps){
  sd.pite$sd[j] <- sd(resultList[[j]]$simData$rf.Pite)
}

pVal.sd <- length(which(sd.pite$sd >= sd.pite$sd[1]))/(numReps)

ggplot(sd.pite, aes(sd))+ geom_density() +
  geom_vline(aes(xintercept=sd.pite$sd[1]), color="red", linetype="dashed") +
  labs(x="SD of PITE_null", title="RF PITE - ")



###  calculate RMSE

for (j in 1:numReps){
  rmse.pite$mse[j] <- sqrt(mean( (resultList[[j]]$simData$avgTE - resultList[[j]]$simData$rf.Pite)^2   ))
}

C <- sort(rmse.pite[rmse.pite$mse > quantile(rmse.pite$mse,prob=1-5/100),1])[1]
  
pVal.mse <- length(which(rmse.pite$mse < C))/(numReps)

ggplot(rmse.pite, aes(mse))+ geom_density() +
  geom_vline(aes(xintercept=rmse.pite$mse[1]), color="red", linetype="dashed")+
  labs(x="RMSE of PITE_null with No Retraining", title="RF PITE Null density with True Reference")





#### Combined Plot
sd.melt <- melt(sd.pite)
rmse.melt <- melt(rmse.pite)

ggplot(sd.melt,aes(x=value, fill=variable)) + geom_density(alpha=0.25)+
  labs(x="SD of PITE_null", title="Sampled cases from Null Dist with True Reference Null")

ggplot(rmse.melt,aes(x=value, fill=variable)) + geom_density(alpha=0.25)+
  labs(x="RMSE of PITE_null", title="Sampled cases from Null Dist with True Reference Null")




#### Power Test

rmse.100 <-rmse.250<-rmse.500<-rmse.1k<-rmse.5k<-rmse.10k<-data.frame(val=numeric(numReps))

for (j in 1:numReps){
  rmse.10k$val[j] <- sqrt(mean( (resultList[[j]]$simData$avgTE - resultList[[j]]$simData$rf.Pite)^2   ))
}

p100 <- length( which(rmse.100 >= rmse.10k ) ) / 1000
p250 <- length( which(rmse.250 >= rmse.10k ) ) / 1000
p500 <- length( which(rmse.500 >= rmse.10k ) ) / 1000
p1k  <- length( which(rmse.1k  >= rmse.10k ) ) / 1000
p5k  <- length( which(rmse.5k  >= rmse.10k ) ) / 1000



data = data.frame(matrix(ncol=59, nrow=1000))

for (i in 1:12){
  
  temp <- read.csv(file=paste("result",i,".csv",sep=""), header = TRUE)
  write.csv(x=temp, file=paste("result",i+48,".csv",sep=""))
  
}


for (i in 1:59){
  
  data[,i] <- read.csv(file=paste("result",i,".csv",sep=""), header = TRUE)[,2]
  
  
}


pval.mse <- data.frame(p=numeric(59))

for (i in 1:59){
  pval.mse[i,1] <- length(which(data[,i] < data[1,i]))/(numReps)
}

Power= length(which(pval.mse$p < .05)) / numReps



ggplot(pval.mse, aes(p))+ geom_density() +
  geom_vline(aes(xintercept=0.05), color="red", linetype="dashed")+
  labs(x="RMSE of PITE_null with No Retraining", title="RF PITE Null density with True Reference")



