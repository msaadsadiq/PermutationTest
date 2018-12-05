library(randomForestSRC)

options(rf.cores = 48, mc.cores= 48)

#Variables
numReps = 1000          #append more repitions to resultList if necessary 
n = 5000
NumCov=5
NumBi=2
NumNuiCov=75
effTrt = effCtrl = c(0.406, -0.239, 0.703, -0.090, -0.299)


## forest parameters
ntree            <- c(10, 50, 100, 250, 500, 1000)[6]
nsplit           <- c(0, 1, 2, 10)[4]
nodesize.grw     <- c(1, 3, 5, 10)[1]

## synthetic forests
mtry <- NumNuiCov
## verbose flag
verbose <- FALSE




generateData <- function(n = 5000, NumCov=5,NumBi=2,NumNuiCov=75,sigma = 1, seed=seed, loopIndex) {
  
  set.seed(2398);
  set.seed(2398);
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
  
  seed=set.seed(seed)
  eps = rnorm(n/2)
  if(loopIndex > 1){
    trtAssgn <- sample(trtAssgn)
  }
  
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


##### Baseline RMSE 
set.seed(740)
seedList=round(runif(numReps,0,n))

baseData <- generateData(n = n,NumCov=5,NumBi=2,NumNuiCov=NumNuiCov, sigma = .1, seed = seedList[10000],1 )

dataTrt  <- baseData[ which(baseData$trtAssgn ==1) ,c( (1:(NumCov+NumNuiCov)),(NumCov+NumNuiCov+4))]     # significant covariates for the conrol group
dataCtrl <- baseData[ which(baseData$trtAssgn ==0) ,c( (1:(NumCov+NumNuiCov)),(NumCov+NumNuiCov+4))]       # significant covariates for the treatment group 

rfCtrl <- rfsrc(Y ~ ., dataCtrl[,1:81],  importance = "none", 
                ntree = ntree, nsplit = nsplit, nodesize = nodesize.grw,
                mtry = mtry, verbose=verbose)

rfTrt <- rfsrc(Y ~ ., dataTrt[,1:81],  importance = "none",
               ntree = ntree, nsplit = nsplit, nodesize = nodesize.grw, 
               mtry = mtry, verbose= verbose)

rf.Yhat.0.0 <- rfCtrl$predicted.oob   #predict(rfCtrl,dataCtrl[,1:(NumCov+NumNuiCov)])$predicted.oob
rf.Yhat.1.1 <- rfTrt$predicted.oob    #predict(rfTrt,dataTrt[,1:(NumCov+NumNuiCov)])$predicted.oob

rf.Yhat.1.0 <- predict(rfTrt, dataCtrl[,1:(NumCov+NumNuiCov)])$predicted
rf.Yhat.0.1 <- predict(rfCtrl, dataTrt[,1:(NumCov+NumNuiCov)])$predicted

rf.Pite.1 =   rf.Yhat.1.1 - rf.Yhat.0.1
rf.Pite.0 =   rf.Yhat.1.0 - rf.Yhat.0.0

baseData$rf.Pite = c(NA)
baseData$rf.Pite[baseData$trtAssgn == 1] <-   rf.Pite.1
baseData$rf.Pite[baseData$trtAssgn == 0] <-   rf.Pite.0



####### Permuted estimate; Training only once

permutedData <- function(n = 5000, NumCov=5,NumBi=2,NumNuiCov=75,sigma = 1, seed=seed, loopIndex) {
  
  set.seed(2398);
  set.seed(2398);
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
  
  seed=set.seed(seed)
  eps = rnorm(n/2)
  trtAssgn <- sample(trtAssgn)
  
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
  retData <- retData[order(retData$trtAssgn, decreasing = TRUE),]
  retData$index <- c(1:n)
  
  return(retData)
}


permData <- permutedData(n = n,NumCov=5,NumBi=2,NumNuiCov=NumNuiCov, sigma = .1, seed = seedList[9999],1 )

dataTrtperm  <- permData[ which(permData$trtAssgn ==1) ,c( (1:(NumCov+NumNuiCov)),(NumCov+NumNuiCov+4),86)]    
dataCtrlperm <- permData[ which(permData$trtAssgn ==0) ,c( (1:(NumCov+NumNuiCov)),(NumCov+NumNuiCov+4),86)]      

permRfCtrl <- rfsrc(Y ~ ., dataCtrlperm[,1:81],  importance = "none", 
                ntree = ntree, nsplit = nsplit, nodesize = nodesize.grw,
                mtry = mtry, verbose=verbose)

permRfTrt <- rfsrc(Y ~ ., dataTrtperm[,1:81],  importance = "none", 
               ntree = ntree, nsplit = nsplit, nodesize = nodesize.grw, 
               mtry = mtry, verbose= verbose)

rf.Yhat.0.0 <- permRfCtrl$predicted.oob  
rf.Yhat.1.1 <- permRfTrt$predicted.oob   

rf.Yhat.1.0 <- predict(permRfTrt, dataCtrlperm[,1:(NumCov+NumNuiCov)])$predicted
rf.Yhat.0.1 <- predict(permRfCtrl, dataTrtperm[,1:(NumCov+NumNuiCov)])$predicted

rf.Pite.1 =   rf.Yhat.1.1 - rf.Yhat.0.1
rf.Pite.0 =   rf.Yhat.1.0 - rf.Yhat.0.0

permData$rf.Pite = c(NA)
permData$rf.Pite[permData$trtAssgn == 1] <-   rf.Pite.1
permData$rf.Pite[permData$trtAssgn == 0] <-   rf.Pite.0


##### Looping without Training

resultList <- lapply(1:numReps, function(loopIndex) {
  
  cat("Replication:", loopIndex,"\n")
  simData <- permData
  simData$trtAssgn <- sample(simData$trtAssgn)
    
  dataTrt  <- simData[ 1:2500    ,c( (1:(NumCov+NumNuiCov+1)),(NumCov+NumNuiCov+4))]    
  dataCtrl <- simData[ 2501:5000 ,c( (1:(NumCov+NumNuiCov+1)),(NumCov+NumNuiCov+4))]       
  
  rf.Yhat.1.1 <- rf.Yhat.0.0 <- rf.Yhat.1.0 <- rf.Yhat.0.1 <- NA

  rf.Yhat.1.1[dataTrt$trtAssgn==1] <- permRfTrt$predicted.oob[dataTrt$trtAssgn==1]   
  rf.Yhat.1.1[dataTrt$trtAssgn==0] <- predict(permRfTrt, dataTrt[,1:(NumCov+NumNuiCov)])$predicted[dataTrt$trtAssgn==0]
  
  rf.Yhat.0.0[dataCtrl$trtAssgn==0] <- permRfTrt$predicted.oob[dataCtrl$trtAssgn==0]  
  rf.Yhat.0.0[dataCtrl$trtAssgn==1] <- predict(permRfCtrl, dataCtrl[,1:(NumCov+NumNuiCov)])$predicted[dataCtrl$trtAssgn==1]
  
  
  rf.Yhat.1.0[dataCtrl$trtAssgn==0] <- predict(permRfTrt, dataCtrl[,1:(NumCov+NumNuiCov)])$predicted[dataCtrl$trtAssgn==0]
  rf.Yhat.1.0[dataCtrl$trtAssgn==1] <- permRfTrt$predicted.oob[dataCtrl$trtAssgn==1]
  
  rf.Yhat.0.1[dataTrt$trtAssgn==1] <- predict(permRfCtrl, dataTrt[,1:(NumCov+NumNuiCov)])$predicted[dataTrt$trtAssgn==0]
  rf.Yhat.0.1[dataTrt$trtAssgn==0] <- permRfCtrl$predicted.oob[dataTrt$trtAssgn==1]
  
  ### subtract
  rf.Pite.1 =   rf.Yhat.1.1 - rf.Yhat.0.1
  rf.Pite.0 =   rf.Yhat.1.0 - rf.Yhat.0.0
    
  simData$rf.Pite = c(NA)
  simData$rf.Pite[simData$trtAssgn == 1] <-   rf.Pite.1
  simData$rf.Pite[simData$trtAssgn == 0] <-   rf.Pite.0

  #### return the goodies
  list(
    #bartPITE = bartPITE,
    #bartSdPITE = bartSdPITE,
    simData = simData
  )
  
})

