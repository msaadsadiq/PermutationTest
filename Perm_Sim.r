
library(parallel)
library(randomForestSRC)
library(BayesTree)
# Programming Convections
# Variable / function def. use camel casing: numCov, calcMean 
# Class / object def. use first Caps casing: UserClass, UserObject


# Models to try 
# 1. GLM
# 2. BART
# 3. Random Forest

# Number of Nuisance variables
# 0, 75 and 150



options(rf.cores = 48, mc.cores= 48)

#Variables
numReps = 1000          #append more repitions to resultList if necessary 
n = 100
NumCov=5
NumBi=2
NumNuiCov=75
effTrt = effCtrl = c(0.406, -0.239, 0.703, -0.090, -0.299)


## forest parameters
ntree            <- c(10, 50, 100, 250, 500, 1000)[6]
nsplit           <- c(0, 1, 2, 10)[4]
nodesize.grw     <- c(1, 3, 5, 10)[1]

## synthetic forests
mtry <- 40
## verbose flag
verbose <- FALSE




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


for(outLoop in 1:500){

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


outer <- 1
index <- seq(1,1000, by=1)
f <-  function(index) {
  #libloc = c("C:/Users/mss255/Documents/R/win-library/3.3")
  #library("BayesTree", lib.loc = libloc)
  library(BayesTree)
  NumCov=5
  NumBi=2
  NumNuiCov=75
  n = 500
  
  
  set.seed(.GlobalEnv$outer);
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
  
  
  seedList=round(runif(1000,0,99999))
  seed=set.seed(seedList[index])
  eps = rnorm(n/2)
  
  # the small difference in yc, yt apart from 0.5 is random noise representing differences due to chance
  yt <-  0.5 + 0.406*x[(1:(n/2)),1]   - 0.239*x[(1:(n/2)),2]   + 0.703*x[(1:(n/2)),3]   - 0.090*x[(1:(n/2)),4]   - 0.299*x[(1:(n/2)),5]   + rnorm(n/2)
  yc <-        0.406*x[((n/2+1):n),1] - 0.239*x[((n/2+1):n),2] + 0.703*x[((n/2+1):n),3] - 0.090*x[((n/2+1):n),4] - 0.299*x[((n/2+1):n),5] + rnorm(n/2)
  
  simData <- data.frame(x)
  simData$trtAssgn <- trtAssgn
  simData$yc <- simData$yt <- rep(NA, n)
  simData$yc[simData$trtAssgn == 0] <- yc
  simData$yt[simData$trtAssgn == 1] <- yt
  simData$Y[simData$trtAssgn == 0] <- yc
  simData$Y[simData$trtAssgn == 1] <- yt
  simData$avgTE <-  mean(yt)-mean(yc)
  
  if(x > 1){
    simData$trtAssgn <- sample(simData$trtAssgn)
  }
  
  W <- as.numeric(simData$trt)
  X <- simData[, 1:80]
  Y <- simData$Y
  len <- length(Y)
  
  
  bartFit <- bart(data.frame(W=W, X),
                  Y,
                  verbose = FALSE,
                  x.test = rbind(data.frame(W=1,X), data.frame(W=0, X)),
                  ntree = 1000)
  
  
  simData$bt.pite <- bartFit$yhat.test.mean[1:len] - bartFit$yhat.test.mean[-(1:len)]
  
  return(data.frame(simData$bt.pite, simData$avgTE))
}


for(outLoop in 1:500){
  
  cat("outLoop:", outLoop,"\n") 

  cl <- makeCluster(16)
  bartRes <- parLapply(cl, index, f)
  stopCluster(cl)

  outer <- outer+1
  rmse.pite <- data.frame(mse=numeric(1000))

  for (unListing in 1:1000){
    rmse.pite$mse[unListing] <- sqrt(mean( (bartRes[[unListing]]$simData.avgTE - bartRes[[unListing]]$simData.bt.pite)^2   ))
  }
  write.csv(rmse.pite$mse, file=paste("result",outLoop,".csv",sep=""))

} #for loop closed










