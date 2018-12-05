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



pite_gen<-function(n, seeds, NumNuiCov, NumCov, NumBi, ATE, EfTRT, EfCntl){   # n: total sample size
  set.seed(seeds)
  if(n <= (NumCov + NumNuiCov)*2+2) 
    stop("degrees of freedom of residual is negative, sample size (n) needs to be larger than twice of the number of total covariates (NumCov + NumNuiCov)")
  if(length(EfTRT)!=NumCov) 
    stop("the length of treatment effect (EfTRT) vector is not the same as the number of non-nuisance covariates (NumCov)")
  if(length(EfCntl)!= NumCov) 
    stop("the length of control effect(EfCntl) vector is not the same as the number of covariates (NumCov)")
  if(NumBi > NumCov)
    stop("Number of non-nuisance binary variables is larger than the total number of non-nuisance variables(NumCov)")

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
 
  xt <- x[(1:(n/2)),(1:NumCov)]       # significant covariates for the treatment group 
  xc <- x[((n/2+1):n),(1:NumCov)]     # significant covariates for the conrol group
  
  trt <- rep(c(1,0), each = n/2)     # treatment condition variable: half treatment, half control
  
  y <- numeric(n)
  y[1:(n/2)] <- yt <- ATE*trt[1:(n/2)] + xt%*%EfTRT + rnorm(n/2)       # observed y in the treatment group
  y[(n/2+1):n] <- yc <- ATE*trt[(n/2+1):n] + xc%*%EfCntl + rnorm(n/2)  # observed y in the control group 
  
  data <- data.frame(y, trt, x)
  
  mt <- lm(yt ~ x[1:(n/2),])                        # fit the model for the treatment group
  mc <- lm(yc ~ x[((n/2+1):n),])                    # fit the model for the control group 
  
  et <- coefficients(mt)           # extract coefficients from the treated model
  ec <- coefficients(mc)           # extract coefficients from the control model 
  
  pred_t <-  c(et[1] + (x[((n/2+1):n),]%*%et[2:(NumCov+NumNuiCov+1)]), et[1] + (x[1:(n/2),]%*%et[2:(NumCov+NumNuiCov+1)]))  # the predicted value to the treatment condition    pite <- pred_e-pred_c 
  pred_c <-  c(ec[1] + (x[((n/2+1):n),]%*%ec[2:(NumCov+NumNuiCov+1)]), ec[1] + (x[1:(n/2),]%*%ec[2:(NumCov+NumNuiCov+1)]))  # the predicted value to the control condition
  pite.ind <- pred_t - pred_c
  pite.sd <- sd(pite.ind)
  return(list("pite.ind" = pite.ind, 
              "pite.sd" = pite.sd, 
              "pite.trt" = pite.ind[1:(n/2)], 
              "pite.cntl" = pite.ind[(n/2+1):n],
              "gen.data" = data, 
              "NumCov" = NumCov,
              "NumNuiCov" = NumNuiCov,
              "gen.x" = x))
}

# ==== test =====

EfTRT = EfCntl = c(0.406, -0.239, 0.703, -0.090, -0.299)


