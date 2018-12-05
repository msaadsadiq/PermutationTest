### Power runs


library(randomForestSRC)

options(rf.cores = 48, mc.cores= 48)

#Variables
numReps = 1000          #append more repitions to resultList if necessary 
n = 5000 #100,250,500,1000,5000
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
  if(loopIndex > 1){
    retData$trtAssgn <- sample(retData$trtAssgn)
  }
  
  return(retData)
}



resultList <- lapply(1:numReps, function(loopIndex) {
  
  cat("Replication:", loopIndex,"\n")
  
  set.seed(740)
  seedList=round(runif(numReps,0,n))
  
  simData <- generateData(n = n,NumCov=5,NumBi=2,NumNuiCov=NumNuiCov, sigma = .1, seed = seedList[loopIndex],loopIndex )
  
  dataTrt  <- simData[ which(simData$trtAssgn ==1) ,c( (1:(NumCov+NumNuiCov)),(NumCov+NumNuiCov+4))]     # significant covariates for the conrol group
  dataCtrl <- simData[ which(simData$trtAssgn ==0) ,c( (1:(NumCov+NumNuiCov)),(NumCov+NumNuiCov+4))]       # significant covariates for the treatment group 
  
  rfCtrl <- rfsrc(Y ~ ., dataCtrl,  importance = "none", na.action = "na.omit",
                  ntree = ntree, nsplit = nsplit, nodesize = nodesize.grw,
                  mtry = mtry, verbose=verbose)
  
  # Regress trt column values
  rfTrt <- rfsrc(Y ~ ., dataTrt,  importance = "none", na.action = "na.omit",
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
  
  
  # seed to all the covariates, so the patients remain the same
  
  # Repeat loop and calculate pite again 
  
  #### return the goodies
  list(
    #bartPITE = bartPITE,
    #bartSdPITE = bartSdPITE,
    simData = simData
  )
  
})

save.image("power_rf_n5000_rep1k.RData")




