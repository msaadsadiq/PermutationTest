

### Simulation that BART fails and CF and MI performs well

library(BayesTree)
library(randomForestSRC)
library(mi)
library(betareg)
library(MASS)
library(ggplot2)
library(gridExtra)

### Alot of nuiscance variables affect them equally
### Compared all simulations that we have done

N=10000
n = 10000
sigma = .1
case = 2
options(rf.cores=48, mc.cores=48)

nrep <- c(100)
n.iter = 30
n.chains = 100
## forest parameters
ntree            <- c(10, 50, 100, 250, 500, 1000)[6]
nsplit           <- c(0, 1, 2, 10)[3]
nodesize.grw     <- c(1, 3, 5, 10)[2]

## synthetic forests
nodesizeSeq      <- c(1:10, 20, 30, 50, 100)
mtrySeq          <- c(1, 10, 12,15,18)
mtry <-16
## verbose flag
verbose <- FALSE

xt <- data.frame(mvrnorm(10000, mu = rep(0, 16), diag(1, 16)) )
colnames(xt) <- paste("W.", 1:16, sep = "")
XL2 <- ones <- rep(1,N)
XL1 <- zeros <- rep(0,N)

XL1[xt$W.3 >= -0.75 & xt$W.3 <=0.75] <- ones[xt$W.3 >= -0.75 & xt$W.3 <=0.75]
XL2[xt$W.4 >= -0.75 & xt$W.4 <=0.75] <- zeros[xt$W.4 >= -0.75 & xt$W.4 <=0.75]



x.c <- matrix(rnorm(n * 11), n, 11)
x.d <- matrix(rbinom(n * 9, size = 1, prob = .5), n, 9)
x <- data.frame(x.c, x.d)
colnames(x) <- paste("W.",1:ncol(x),sep="")

## exposure variables
F <- (-2.0 + .028 * x[, 1] - 0.374 * x[, 2] -.03 * x[, 3] + .118 * x[, 4] 
      - 0.394 * x[, 11] + .875 * x[, 12] + .9 * x[, 13])
prob.true <- plogis(F)
x$trt  <- rbinom(n, size = 1, prob = prob.true)


## noise
eps <- rnorm(n, sd = sigma)

x$y0.true <- (2.455
              - sin(.40 * x[, 1] + .154 * x[, 2] * XL2 - .152 * x[, 11] - .126 * x[, 12] * XL1) )
x$y1.true <- (2.455
              - 1 * ((.254 * (x[, 2] ^ 2) * XL2 - .152 * x[, 11] - .4 * x[, 11]^2 - .126 * x[, 12] * XL1) > 0) )
x$y <- (2.455
        - (1 - x$trt) * sin(.4 * x[, 1] + .154 * x[, 2] * XL2 - .152 * x[, 11] - .126 * x[, 12] )
        - x$trt * (1 * (.254 * (x[, 2] ^ 2) * XL2 - .152 * x[, 11] - .4 * x[, 11]^2 - .126 * x[, 12] * XL1) > 0)
        + eps)      


x$y0 <- x$y1 <- rep(NA, n)
x$y0[x$trt == 0] <- x$y[x$trt == 0]
x$y1[x$trt == 1] <- x$y[x$trt == 1]

x$TE <- x$y1.true - x$y0.true




dataTrt  <- x[x$trt ==1, c(25,1:20)]
dataCtrl <- x[x$trt ==0, c(26,1:20)]

dataTrtMI <- x[, c(25,1:20)]
dataCtrlMI <- x[, c(26,1:20)]

##### MI

mdf.cntrl <- missing_data.frame(dataCtrlMI)
mdf.trt <- missing_data.frame(dataTrtMI)

# Impute missing values using mi in control column
mi.trt.impute  <- mi(mdf.cntrl, n.chains =n.chains, n.iter= n.iter, verbose=FALSE) #n.iter=100
# Impute missing values by mi in treatment column
mi.cntrl.impute  <- mi(mdf.trt, n.chains =n.chains, n.iter= n.iter, verbose=FALSE) #max.minutes=20000

pool.cntrl.coeff <-  pool(y0 ~  W.1 + W.2 + W.3 + W.4 + W.5 + W.6 + W.7 + W.8 + W.9 + W.10 +
                            W.11 + W.12 + W.13 + W.14 + W.15 + W.16 + W.17 + W.18 + W.19 + W.20 
                          , mi.trt.impute)

pool.trt.coeff   <-  pool(y1 ~ W.1 + W.2 + W.3 + W.4 + W.5 + W.6 + W.7 + W.8 + W.9 + W.10 +
                            W.11 + W.12 + W.13 + W.14 + W.15 + W.16 + W.17 + W.18 + W.19 + W.20 
                          , mi.cntrl.impute)

cntrl.coefficients <- summary(pool.cntrl.coeff)$coefficients[2:21,1]
trt.coefficients   <- summary(pool.trt.coeff)$coefficients[2:21,1]




mi.Yhat.0.1 = as.matrix(x[5001:10000,1:20]) %*% as.vector(cntrl.coefficients)
mi.Yhat.0.0 = as.matrix(x[1:5000,1:20]) %*% as.vector(cntrl.coefficients)

mi.Yhat.1.1 = as.matrix(x[5001:10000,1:20]) %*% as.vector(trt.coefficients)
mi.Yhat.1.0 = as.matrix(x[1:5000,1:20]) %*% as.vector(trt.coefficients)

mi.pite.0 <- mi.Yhat.1.0 - mi.Yhat.0.0
mi.pite.1 <- mi.Yhat.1.1 - mi.Yhat.0.1

# combining the two pites together

x$mi.pite <- c(mi.pite.0,mi.pite.1)

x$mi.bias <-  x$mi.pite - x$TE[5001:10000] # bias = TE - PITE
x$RMSE.mi <-  (x$mi.pite - x$TE[5001:10000]) ^2


##### Syn Forest

rf.cntrl.est  <-   rfsrcSyn(y0 ~ ., dataCtrl, importance = "none", na.action = "na.omit",
                            ntree = ntree, nsplit = nsplit, nodesizeSeq = nodesizeSeq,
                            mtrySeq = mtrySeq, nodesize = nodesize.grw, oob=TRUE,
                            verbose = verbose)

# Regress trt column values
rf.trt.est    <-   rfsrcSyn(y1 ~ ., dataTrt, importance = "none", na.action= "na.omit",
                            ntree = ntree, nsplit = nsplit, nodesizeSeq = nodesizeSeq,
                            mtrySeq = mtrySeq, nodesize = nodesize.grw, oob=TRUE,
                            verbose = verbose)


Syn.rf.Yhat.0.0 <-  rf.cntrl.est$rfSyn$predicted.oob #rfsrcSyn(object = rf.cntrl.est, newdata = dataCtrl)$rfSyn$predicted.oob
Syn.rf.Yhat.1.1 <-  rf.trt.est$rfSyn$predicted.oob

Syn.rf.Yhat.1.0   <- rfsrcSyn(object = rf.trt.est, newdata = dataCtrl)$rfSynPred$predicted
Syn.rf.Yhat.0.1   <- rfsrcSyn(object = rf.cntrl.est, newdata = dataTrt)$rfSynPred$predicted

Syn.rf.pite.0 = Syn.rf.Yhat.1.0 - Syn.rf.Yhat.0.0
Syn.rf.pite.1 = Syn.rf.Yhat.1.1 - Syn.rf.Yhat.0.1

x$Syn.rf.pite <- c(NA)
x$Syn.rf.pite[x$trt ==0 ] <- Syn.rf.pite.0
x$Syn.rf.pite[x$trt ==1 ] <- Syn.rf.pite.1

x$Syn.rf.bias     <-  x$Syn.rf.pite - x$TE
x$Syn.RMSE.rf     <-  (x$Syn.rf.pite - x$TE) ^2

#### Forest
rf.cntrl.est <- rfsrc(y0 ~ ., dataCtrl,  importance = "none", na.action = "na.omit",
                      ntree = ntree, nsplit = nsplit, nodesize = nodesize.grw,
                      mtry = mtry, verbose=verbose, oob=TRUE)

# Regress trt column values
rf.trt.est <- rfsrc(y1 ~ ., dataTrt,  importance = "none",
                    ntree = ntree, nsplit = nsplit, nodesize = nodesize.grw, na.action = "na.omit",
                    mtry = mtry, verbose= verbose, oob=TRUE)

# I calculate all four Yhats because synthetic forest have different method for prediction
rf.Yhat.0.0 <- rf.cntrl.est$predicted.oob  #predict(rf.cntrl.est,data.0.test[,1:20])$predicted.oob
rf.Yhat.1.1 <- rf.trt.est$predicted.oob #predict(rf.trt.est,data.1.test[,1:20])$predicted.oob

rf.Yhat.1.0 <- predict(rf.trt.est, dataCtrl)$predicted
rf.Yhat.0.1 <- predict(rf.cntrl.est, dataTrt)$predicted


rf.pite.0 = rf.Yhat.1.0 - rf.Yhat.0.0
rf.pite.1 = rf.Yhat.1.1 - rf.Yhat.0.1

x$rf.pite <- c(NA)
x$rf.pite[x$trt==0] <- rf.pite.0
x$rf.pite[x$trt==1] <- rf.pite.1

x$rf.bias     <-  x$rf.pite - x$TE
x$RMSE.rf     <-  (x$rf.pite - x$TE ) ^2

### BART

dta <- x
W <- as.numeric(dta$trt)
X <- dta[, 1:20]
Y <- dta$y
len <- length(Y)


bartFit <- bart(data.frame(W=W, X),
                Y,
                verbose = FALSE,
                x.test = rbind(data.frame(W=1,X), data.frame(W=0, X)),
                ntree = 1000)


x$bt.pite <- bartFit$yhat.test.mean[1:len] - bartFit$yhat.test.mean[-(1:len)]
x$bt.bias <- x$bt.pite - x$TE
x$bt.mse <- (x$bt.pite - x$TE)^2




#### Group plots, all three methods in same window

### Plot Bias Density

firstplot <-   ggplot(x, aes(x$rf.bias)) +
  geom_density(alpha = 0.1) + 
  labs(x = "RF.Bias", y="Density")

secondplot <-   ggplot(x, aes(x$Syn.rf.bias)) +
  geom_density(alpha = 0.1) + 
  labs(x = "Syn.RF.Bias", y="Density")


thirdplot <-    ggplot(x, aes(x$mi.bias)) +
  geom_density(alpha = 0.1) + 
  labs(x = "MI.Bias", y="Density")

fourthplot <-    ggplot(x, aes(x$bt.bias)) +
  geom_density(alpha = 0.1) + 
  labs(x = "BT.Bias", y="Density")

grid.arrange(firstplot, secondplot, thirdplot,fourthplot)


### Plot MSE Density

firstplot <-   ggplot(x, aes(x$RMSE.rf)) +
  geom_density(alpha = 0.1) + 
  labs(x = "RF.RMSE", y="Density")

secondplot <-   ggplot(x, aes(x$Syn.RMSE.rf)) +
  geom_density(alpha = 0.1) + 
  labs(x = "Syn.RF.RMSE", y="Density")


thirdplot <-    ggplot(x, aes(x$RMSE.mi)) +
  geom_density(alpha = 0.1) + 
  labs(x = "MI.RMSE", y="Density")

fourthplot <-    ggplot(x, aes(x$bt.mse)) +
  geom_density(alpha = 0.1) + 
  labs(x = "BT.RMSE", y="Density")

grid.arrange(firstplot, secondplot, thirdplot,fourthplot)

### Pite vs TE

firstplot =  ggplot(x, aes(x$TE, x$rf.pite)) +
  geom_point(alpha = 0.1) + geom_smooth(method = "loess",level=0.9, span=0.3) +
  theme(legend.position="none") + labs(y = "RF.PITE", x="TE")

secondplot =  ggplot(x, aes(x$TE, x$Syn.rf.pite)) +
  geom_point(alpha = 0.1) + geom_smooth(method = "loess",level=0.9, span=0.3) +
  theme(legend.position="none") + labs(y = "Syn.RF.PITE", x="TE")


thirdplot =  ggplot(x, aes(x$TE, x$mi.pite)) +
  geom_point(alpha = 0.1) + geom_smooth(method = "loess",level=0.9, span=0.3) +
  theme(legend.position="none") + labs(y = "MI.PITE", x="TE")

fourthplot =  ggplot(x, aes(x$TE, x$bt.pite)) +
  geom_point(alpha = 0.1) + geom_smooth(method = "loess",level=0.9, span=0.3) +
  theme(legend.position="none") + labs(y = "BT.PITE", x="TE")

grid.arrange(firstplot, secondplot, thirdplot,fourthplot)


### Bias vs TE

firstplot =  ggplot(x, aes(x$TE, x$rf.bias)) +
  geom_point(alpha = 0.1) + geom_smooth(method = "loess",level=0.9, span=0.3) +
  theme(legend.position="none") + labs(y = "RF.Bias", x="TE")

secondplot =  ggplot(x, aes(x$TE, x$Syn.rf.bias)) +
  geom_point(alpha = 0.1) + geom_smooth(method = "loess",level=0.9, span=0.3) +
  theme(legend.position="none") + labs(y = "Syn.RF.Bias", x="TE")


thirdplot =  ggplot(x, aes(x$TE, x$mi.bias)) +
  geom_point(alpha = 0.1) + geom_smooth(method = "loess",level=0.9, span=0.3) +
  theme(legend.position="none") + labs(y = "MI.Bias", x="TE")

fourthplot =  ggplot(x, aes(x$TE, x$bt.bias)) +
  geom_point(alpha = 0.1) + geom_smooth(method = "loess",level=0.9, span=0.3) +
  theme(legend.position="none") + labs(y = "BT.Bias", x="TE")

grid.arrange(firstplot, secondplot, thirdplot,fourthplot)



