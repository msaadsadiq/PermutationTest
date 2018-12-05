
library(randomForestSRC)
library(ggplot2)

"
# Programming Convections
# Variable / function def. use camel casing: numCov, calcMean 
# Class / object def. use first Caps casing: UserClass, UserObject
#
# Models to try 
# 1. GLM
# 2. BART
# 3. Random Forest
#
# Number of Nuisance variables
# 0, 75 and 150
"

#Variable definitions

options(rf.cores = -1L, mc.cores= -1L) #setting libraries to use all available Cores

numReps = 1000         
n = 100
NumCov=5       # Total number of Covariates
NumBi=2        # Number of Binary Covariates
NumNuiCov=c(75,150)[1]   # Number of Nuisance Covariates
effTrt = effCtrl = c(0.406, -0.239, 0.703, -0.090, -0.299)  # Beta parameters

datasets <- 600

## forest parameters
ntree            <- c(10, 50, 100, 250, 500, 1000)[6]
nsplit           <- c(0, 1, 2, 10)[4]
nodesize.grw     <- c(1, 3, 5, 10)[1]

mtry <- 40
## verbose flag
verbose <- FALSE

ind.pite <- data.frame(one=numeric(datasets),two=numeric(datasets),five=numeric(datasets),thousand=numeric(datasets))
rmse.pite <- data.frame(one=numeric(datasets),two=numeric(datasets),five=numeric(datasets),thousand=numeric(datasets))
rmse <- data.frame(one=numeric(datasets),two=numeric(datasets),five=numeric(datasets),thousand=numeric(datasets))



#### Data generation function

generateData <- function(n = 5000, NumCov=5,NumBi=2,NumNuiCov=75,sigma = 1, seed=seed, loopIndex,outLoop) {
  
  set.seed(outLoop*2398);
  x <- matrix(rep(0), nrow = n, ncol = NumCov + NumNuiCov)
  x[, 1:NumBi] <- matrix(rbinom(n*NumBi, 1, 0.5), ncol = NumBi)
  x[, (NumBi+1):NumCov] <- matrix(rnorm(n*(NumCov-NumBi), 0,1), ncol = NumCov - NumBi)
  
  if(NumNuiCov == 75) {
    x[, (NumCov + 1):(NumCov + 35)] <- matrix(rbinom(n*35, 1, 0.5), ncol = 35)
    x[, (NumCov + 36):(NumCov + NumNuiCov)] <- matrix(rnorm(n*40, 0, 1), ncol = 40)
  }
  
  if(NumNuiCov == 150) {
    x[, (NumCov + 1):(NumCov + 70)] <- matrix(rbinom(n*70, 1, 0.5), ncol = 70)
    x[, (NumCov + 71):(NumCov + NumNuiCov)] <- matrix(rnorm(n*80, 0, 1), ncol = 80)
  }
  
  trtAssgn <-rep(c(1,0), each = n/2)   #  outiside seed, Randomly permute the trt assignment
  if(loopIndex > 1){
    trtAssgn <- sample(trtAssgn)
  }
  
  
  seed=set.seed(seed)
  eps = rnorm(n/2)
  
  
  # the small difference in yc, yt apart from 0.5 is random noise representing differences due to chance
  yt <-  0.5 + 0.406*x[(1:(n/2)),1]   - 0.239*x[(1:(n/2)),2]   + 0.703*x[(1:(n/2)),3]   - 0.090*x[(1:(n/2)),4]   - 0.299*x[(1:(n/2)),5]   + rnorm(n/2)
  yc <-        0.406*x[((n/2+1):n),1] - 0.239*x[((n/2+1):n),2] + 0.703*x[((n/2+1):n),3] - 0.090*x[((n/2+1):n),4] - 0.299*x[((n/2+1):n),5] + rnorm(n/2)
  
  retData <- data.frame(x)
  retData$trtAssgn <- trtAssgn
  retData$yc <- retData$yt <- rep(NA, n)
  retData$yc[retData$trtAssgn == 0] <- yc
  retData$yt[retData$trtAssgn == 1] <- yt
  retData$Y[retData$trtAssgn == 0] <- yc
  retData$Y[retData$trtAssgn == 1] <- yt
  retData$avgTE <-  mean(yt)-mean(yc)
  
  
 
  
  return(retData)
}


#### Algorithm to examine Type1 and Type2 error

for(outLoop in 1:500){  # repeat over 500 datasets
  
  cat("outLoop:", outLoop,"\n")
  
  resultList<- lapply(1:numReps, function(loopIndex) {
    
    cat("Replication:", loopIndex,"\n")
    
    set.seed(740*outLoop)
    seedList=round(runif(numReps,0,n))
    
    simData <- generateData(n = n,NumCov=5,NumBi=2,NumNuiCov=NumNuiCov, sigma = .1, seed = seedList[loopIndex],loopIndex,outLoop )
    
    dataTrt  <- simData[ which(simData$trtAssgn ==1) ,c( (1:(NumCov+NumNuiCov)),(NumCov+NumNuiCov+4))]     # significant covariates for the conrol group
    dataCtrl <- simData[ which(simData$trtAssgn ==0) ,c( (1:(NumCov+NumNuiCov)),(NumCov+NumNuiCov+4))]       # significant covariates for the treatment group 
    
    rfCtrl <- rfsrc(Y ~ ., dataCtrl,  importance = "none", 
                    ntree = ntree, nsplit = nsplit, nodesize = nodesize.grw,
                    mtry = mtry, verbose=verbose)
    
    # Regress trt column values
    rfTrt <- rfsrc(Y ~ ., dataTrt,  importance = "none",
                   ntree = ntree, nsplit = nsplit, nodesize = nodesize.grw, 
                   mtry = mtry, verbose= verbose)
    
    # I calculate all four Yhats because synthetic forest have different method for prediction
    rf.Yhat.0.0 <- rfCtrl$predicted.oob   #predict(rfCtrl,dataCtrl[,1:(NumCov+NumNuiCov)])$predicted.oob
    rf.Yhat.1.1 <- rfTrt$predicted.oob    #predict(rfTrt,dataTrt[,1:(NumCov+NumNuiCov)])$predicted.oob
    
    rf.Yhat.1.0 <- predict(rfTrt, dataCtrl[,1:(NumCov+NumNuiCov)])$predicted
    rf.Yhat.0.1 <- predict(rfCtrl, dataTrt[,1:(NumCov+NumNuiCov)])$predicted
    
    rf.Pite.1 =   rf.Yhat.1.1 - rf.Yhat.0.1
    rf.Pite.0 =   rf.Yhat.1.0 - rf.Yhat.0.0
    
    simData$rf.Pite = c(NA)
    simData$rf.Pite[simData$trtAssgn == 1] <-   rf.Pite.1
    simData$rf.Pite[simData$trtAssgn == 0] <-   rf.Pite.0
    
    simData <- simData[,85:86]
    
    
    # seed to all the covariates, so the patients remain the same
    
    # Repeat loop and calculate pite again 
    
    #### return the goodies
    list(
      simData = simData
    )
    
  }
  ) # lapply closed
  
  # Gather the 1000 permutation results
  rmse.pite <- data.frame(mse=numeric(numReps))
  
  for (unListing in 1:numReps){
    rmse.pite$mse[unListing] <- sqrt(mean( (resultList[[unListing]]$simData$avgTE - resultList[[unListing]]$simData$rf.Pite)^2   ))
  }
  write.csv(rmse.pite$mse, file=paste("result",outLoop,".csv",sep=""))
  
  # save.image("rfsrc_nu75_n5k_rep1k_Y_outLoop500.RData")
} #for loop closed





for (i in 1:600){
  dummy <-  read.csv(paste("C:/Users/mss255/Box Sync/SAAD/Permutation_Test/PowerLoop_n100_mtry40/result",i,".csv",sep=""),header=TRUE)[2]
  rmse.pite[i,1] <- length(which(dummy[2:1000,1]<dummy[1,1]))/numReps
  ind.pite[i,1] <- (rmse.pite[i,1] < 0.05)*1
  
  dummy <- read.csv(paste("C:/Users/mss255/Box Sync/SAAD/Permutation_Test/PowerLoop_n250_mtry40/",i,".csv",sep=""),header=TRUE)[2]
  rmse.pite[i,2] <- length(which(dummy[2:1000,1]<dummy[1,1]))/numReps
  ind.pite[i,2] <- (rmse.pite[i,2] <  0.05)*1
  
  dummy <- read.csv(paste("C:/Users/mss255/Box Sync/SAAD/Permutation_Test/PowerLoop_n500_mtry40/",i,".csv",sep=""),header=TRUE)[2]
  rmse.pite[i,3] <- length(which(dummy[2:1000,1]<dummy[1,1]))/numReps
  ind.pite[i,3] <- (rmse.pite[i,3] < 0.05)*1

#  dummy <-  read.csv(paste("C:/Users/mss255/Box Sync/SAAD/Permutation_Test/n1000_mtry40/",i,".csv",sep=""),header=TRUE)[3]
#  rmse.pite[i,4] <- length(which(dummy[1,1]>dummy[2:1000,1]))/numReps
#  ind.pite[i,4] <- (rmse.pite[i,4] > 0.05)*1
}

x1 <- mean(ind.pite[,1])
x2 <- mean(ind.pite[,2])
x3 <- mean(ind.pite[,3])
x4 <- mean(ind.pite[,4])

sampleSize <- c(100,250,500)#,1000)
Power <- c(x1,x2,x3)#,x4)


plot(sampleSize,Power,type="o")


for (i in 1:600){
  rmse[i,1] <-  read.csv(paste("C:/Users/mss255/Box Sync/SAAD/Permutation_Test/PowerLoop_n100_mtry40/result",i,".csv",sep=""),header=TRUE)[2]
  rmse[i,2]<- read.csv(paste("C:/Users/mss255/Box Sync/SAAD/Permutation_Test/PowerLoop_n250_mtry40/",i,".csv",sep=""),header=TRUE)[2]
  rmse[i,3]<- read.csv(paste("C:/Users/mss255/Box Sync/SAAD/Permutation_Test/PowerLoop_n500_mtry40/",i,".csv",sep=""),header=TRUE)[2]
}