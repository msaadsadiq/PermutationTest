
generateData <- function(n = 100000, NumCov=5,NumBi=2,NumNuiCov=0,sigma = 1, seed=seed, loopIndex) {
  
  set.seed(2398);
  x <- matrix(rep(0), nrow = n, ncol = NumCov + NumNuiCov)
  x[, 1:NumBi] <- matrix(rbinom(n*NumBi, 1, 0.5), ncol = NumBi)
  x[, (NumBi+1):NumCov] <- matrix(rnorm(n*(NumCov-NumBi), 0,1), ncol = NumCov - NumBi)
  
  if(NumNuiCov == 75) {
    x[, (NumCov + 1):(NumCov + 35)] <- matrix(rbinom(n*35, 1, 0.5), ncol = 35)
    x[, (NumCov + 36):(NumCov + NumNuiCov)] <- matrix(rnorm(n*40, 0, 1), ncol = 40)
  }
  
  trtAssgn <-rep(c(1,0), each = n/2)   #  outiside seed, Randomly permute the trt assignment
  
  seed=set.seed(seed)
  eps = rnorm(n/2)
  if(loopIndex != 1){
    trtAssgn <- sample(trtAssgn)
  }
  
  
  yt <- 0.5 + 0.406*x[((n/2+1):n),1] - 0.239*x[((n/2+1):n),2] + 0.703*x[((n/2+1):n),3] - 0.090*x[((n/2+1):n),4] - 0.299*x[((n/2+1):n),5] + eps
  
  yc <-   0.406*x[(1:(n/2)),1] - 0.239*x[(1:(n/2)),2] + 0.703*x[(1:(n/2)),3] - 0.090*x[(1:(n/2)),4] - 0.299*x[(1:(n/2)),5] + eps
  
  
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
