########################################################################
# pite.gen function is used to generate data for PITE approach.        #
#    n: Total sample size                                              #
#    seeds: a vector of seeds                                          #
#    NumNuiCov: Number of Nuisance covariates (can be 0 or 75 or 150)  #                        
#    NumCov: the total number of non-nuisance covariates               #
#    NumBi: Number of non-nuisance binary covariates having effects    #
#    (NumCov - NumBi) : Number of non-nuisance continuous covariates   # 
#    EfTRT: Covariate effects of the treatment group                   #
#    EfCntl: Covariate effects of the control group                    #
########################################################################

# Programming Convections
# Variable / function def. use camel casing: numCov, calcMean 
# Class / object def. use first Caps casing: UserClass, UserObject


# ==== test =====

n=10000
numReps=400
set.seed(940)
seedList=round(runif(numReps,0,n))
NumNuiCov=75
NumCov=5
NumBi=2
ATE=0.5
EfTRT = EfCntl = c(0.406, -0.239, 0.703, -0.090, -0.299)



r1  <- lapply(1:numReps, function(loopIndex) {
  cat("Replication:", loopIndex,"\n")
  set.seed(740)
  
  x <- matrix(rep(0), nrow = n, ncol = NumCov + NumNuiCov)
  x[, 1:NumBi] <- matrix(rbinom(n*NumBi, 1, 0.5), ncol = NumBi)
  x[, (NumBi+1):NumCov] <- matrix(rnorm(n*(NumCov-NumBi), 0,1), ncol = NumCov - NumBi)
  
  if(NumNuiCov == 75) {
    x[, (NumCov + 1):(NumCov + 35)] <- matrix(rbinom(n*35, 1, 0.5), ncol = 35)
    x[, (NumCov + 36):(NumCov + NumNuiCov)] <- matrix(rnorm(n*40, 0, 1), ncol = 40)
  }
  
  xt <- x[(1:(n/2)),(1:NumCov)]       # significant covariates for the treatment group 
  xc <- x[((n/2+1):n),(1:NumCov)]     # significant covariates for the conrol group
  
  trt <- rep(c(1,0), each = n/2)     # treatment condition variable: half treatment, half control
  
  seed=set.seed(seedList[loopIndex])
 
  y <- numeric(n)
  y[1:(n/2)] <- yt <- ATE*trt[1:(n/2)] + xt%*%EfTRT + rnorm(n/2)       # observed y in the treatment group
  y[(n/2+1):n] <- yc <- ATE*trt[(n/2+1):n] + xc%*%EfCntl + rnorm(n/2)  # observed y in the control group 
  

  trt <- sample(trt)
  
  
  data <- data.frame(y, trt, x)
  
  mt <- lm(yt ~ x[1:(n/2),])                        # fit the model for the treatment group
  mc <- lm(yc ~ x[((n/2+1):n),])                    # fit the model for the control group 
  
  et <- coefficients(mt)           # extract coefficients from the treated model
  ec <- coefficients(mc)           # extract coefficients from the control model 
  
  pred_t <-  c(et[1] + (x[((n/2+1):n),]%*%et[2:(NumCov+NumNuiCov+1)]), et[1] + (x[1:(n/2),]%*%et[2:(NumCov+NumNuiCov+1)]))  # the predicted value to the treatment condition    pite <- pred_e-pred_c 
  pred_c <-  c(ec[1] + (x[((n/2+1):n),]%*%ec[2:(NumCov+NumNuiCov+1)]), ec[1] + (x[1:(n/2),]%*%ec[2:(NumCov+NumNuiCov+1)]))  # the predicted value to the control condition
  pite.ind <- pred_t - pred_c
  pite.sd <- sd(pite.ind)
  
  
  list("pite.ind" = pite.ind, 
       "pite.sd" = pite.sd,
       "gen.data" = data
       ) 
})


save.image("lm_nu75_n10k_rep500.RData")


##### Plotting the Results

sd.pite <- data.frame(sd=numeric(399))
rmse.pite <- data.frame(mse=numeric(300))

pVal =0


for (j in 1:399){
  sd.pite$sd[j] <- rO[[j]]$pite.sd
}

pVal.sd <- length(which(sd.pite$sd > sd.pite$sd[1]))/(399)

ggplot(sd.pite, aes(sd))+
  geom_density() +   geom_vline(aes(xintercept=sd.pite$sd[1]), color="red", linetype="dashed") + 
  labs(x="SD of PITE_null", title="LM PITE Null density with True Reference")

  +geom_vline(aes(xintercept=as.numeric(summary(sd.pite$sd)[1])), color="black", linetype="dashed") +
    geom_vline(aes(xintercept=as.numeric(summary(rmse.pite$mse)[4])), color="black", linetype="dashed") +
    geom_vline(aes(xintercept=as.numeric(summary(rmse.pite$mse)[5])), color="black", linetype="dashed") +


  
  
  
  
  
#################################### 

  for (j in 1:300){
    avgTE <- mean(rO[[j]]$gen.data$y[rO[[j]]$gen.data$trt==1]) - mean(rO[[j]]$gen.data$y[rO[[j]]$gen.data$trt==0])
    rmse.pite$mse[j] <- sqrt(mean( (avgTE - rO[[j]]$pite.ind)^2 ))
  }

pVal.mse <- length(which(rmse.pite$mse > rmse.pite$mse[1]))/(300)



ggplot(rmse.pite, aes(mse))+
  geom_density() +geom_vline(aes(xintercept=as.numeric(summary(rmse.pite$mse)[2])), color="black", linetype="dashed") +
  geom_vline(aes(xintercept=as.numeric(summary(rmse.pite$mse)[4])), color="black", linetype="dashed") +
  geom_vline(aes(xintercept=as.numeric(summary(rmse.pite$mse)[5])), color="black", linetype="dashed") +
  geom_vline(aes(xintercept=pVal.mse), color="green", linetype="dashed")

