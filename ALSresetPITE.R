als = imp_tot2[imp_tot2$.imp == 2, !names(imp_tot2) %in% c(".imp", ".id")]
library(car)

als$Study_Arm = recode(als$Study_Arm, "'Active' = 1; 'Placebo'= 0", as.factor.result = F)
als$DeltaFlag = recode(als$DeltaFlag, "'TRUE' = 1; 'FALSE' = 0", as.factor.result = F)
als$Gender = recode(als$Gender, "'Female' = 1; 'Male' = 0", as.factor.result = F)

write.csv(als, "C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\MERGED&CLEAN\\als.csv", row.names = F)
#
y = als$ALSFRS_slope
cov = als[,!names(als) %in% c("subject_id", "ALSFRS_slope", "Study_Arm")]
trt = als$Study_Arm

data <- data.frame(y, trt)
# estimate the ATE
if (var.test(y~trt)$p.value <= 0.05) {
    tt <- t.test(y~trt, data = data, var.equal = F)$estimate
    ATE.est <- tt[2] - tt[1]
}

if(var.test(y~trt)$p.value > 0.05) {
    tt <- t.test(y~trt, data = data, var.equal = T)$estimate
    ATE.est <- tt[2] - tt[1]
}


# NO INTERACTION 

inter = F
dataxy <- data.frame(y, cov)
# create the interaction between pair predictors x, add in to the dataset
# dataint now have y, treatment condition, x, and the pairwise interaction between x
dataint2 <- model.matrix(y ~ .^2, data = dataxy)
dataint <- data.frame(y, trt, x = dataint2[,-1])

xt <- dataxy[dataint$trt ==1, -1]
xc <- dataxy[dataint$trt ==0, -1]

# remove the estiamted ATE from the treatment group
yt <- y[trt == 1] - ATE.est      # remove the estimated average treatment effect from the observed y in the treatment group
yc <- y[trt == 0]                # subset control group y

# seperate the dataset (dataint) to treatment group and control group
TRTmxi <- cbind(yt, xt)
CNTLmxi <- cbind(yc, xc)

# if the interaction of predictors are not taken into considertion
if (inter == FALSE) {
    mt <- lm(yt ~., data = TRTmxi)
    mc <- lm(yc ~., data = CNTLmxi)
    
    et <- coefficients(mt)                                               # extract coefficients from the treated model
    ec <- coefficients(mc)                                               # extract coefficients from the control model
    pred_t <-  c(et[1] + (as.matrix(CNTLmxi[,2:(ncol(cov)+1)]) %*% as.matrix(et[2:(length(et))])), et[1] + (as.matrix(TRTmxi[,2:(ncol(cov)+1)])%*%as.matrix(et[2:(length(et))])))  # the predicted value to the treatment condition
    pred_c <-  c(ec[1] + (as.matrix(CNTLmxi[,2:(ncol(cov)+1)]) %*% as.matrix(ec[2:(length(ec))])), ec[1] + (as.matrix(TRTmxi[,2:(ncol(cov)+1)])%*%as.matrix(ec[2:(length(ec))])))  # the predicted value to the control condition
}

pite.ind <- pred_t - pred_c
pite.sd <- sd(pite.ind); pite.sd
#write.csv(cbind(et, ec), "C:/Users/Dawn/Dropbox/PITEPermutation/ApplicationALS/CoefNoInt.csv")


#======================Permutation Test

# extract the pite result and the information from the dataset
n <- length(pite.ind)
pite.ind <- pite.ind
pite.sd <- pite.sd
TOTNumCov <- ncol(cov)

# extract y
y <- y

# extract x  
x <- cov

# permute nperm times
nperm = 1000
prmt <- matrix(NA, nrow = n, ncol = nperm)
warn.count = 0

set.seed(9999)
for(i in 1:nperm){
    trt <- sample(trt, replace = F)    # shuffle (permute) the treatment conditions
    xtp <- as.matrix(x[trt == 1, ])                       # get treatment groups' covariates after shuffling
    xcp <- as.matrix(x[trt == 0, ])                       # get control groups' covariates after shuffling
    ytp <- y[trt == 1]                                    # get the y of the treatment
    ycp <- y[trt == 0]                                    # get the y of the control 
    mt <- lm(ytp ~ xtp)                                   # regress the observed y on the shuffled x of the treatment group 
    mc <- lm(ycp ~ xcp)                                   # regress the observed y on the suffled x of the control group 
    et <- coefficients(mt)                                # extract the coefficients from the regression model of the treatment group 
    ec <- coefficients(mc)                                # extract the coefficients from the regression model of the control group 
    
    # if singular condition happened in the dataset, the effects of covariates can be zero
    # record how many times sigular happened, if sigular happened, the loop will skip this permutation 
    # ==========  this is NOT an ideal way to extract sigular issue. 
    # ==========  It might be better to check sigular issue righ after the treatment were shuffled, rather than after the estimation
    if(sum(is.na(et))>0 | sum(is.na(ec))>0){               
        warn.count = warn.count + 1
        next
    }
    
    # calculate PITE for nperm times
    pred_t <-  c(et[1] + (xcp%*%et[2:(length(et))]), et[1] + (xtp%*%et[2:(length(et))]))  # calculate individual's predicted value to treatment
    pred_c <-  c(ec[1] + (xcp%*%ec[2:(length(ec))]), ec[1] + (xtp%*%ec[2:(length(ec))]))  # calcualte individual's predicted valued to control
    prmt[,i] <- pred_t-pred_c                                                             # calculate PITE for the ith permutation replicate
}

# calculate the standard deviatoin of PITEs for each permutation
prmt.sd <- apply(prmt, 2, sd, na.rm = T) 
#par(mfrow = c(1, 1))
par(mfrow = c(1, 2))
hist(prmt.sd, xlim = c(0, 0.5), ylim = c(0, 80), breaks = 50, main = "Histogram of prmt.sd \n No interaction")
abline(v = pite.sd, col = 'red')
text(pite.sd, 80, paste("pite.sd = ", round(pite.sd, 3)), col = 'blue')
text(0.3, 30, paste("# of permutation = ", nperm))
text(0.3, 35, paste("# of singularity occurred = ", warn.count))


############### two way interaction ########################################################################################
inter = T
covNR <- cov[,!names(cov)%in%c("DeltaFlag", "BulbarOnly")]
dataxy <- data.frame(y, covNR)
# create the interaction between pair predictors x, add in to the dataset
# dataint now have y, treatment condition, x, and the pairwise interaction between x
dataint2 <- model.matrix(y ~ .^2, data = dataxy)
dataint2 = as.data.frame(dataint2)
dataint2 = cbind(dataint2, cov$DeltaFlag, cov$BulbarOnly)


xt <- dataint2[trt == 1, -1]
xc <- dataint2[trt == 0, -1]

# remove the estiamted ATE from the treatment group
yt <- y[trt == 1] - ATE.est      # remove the estimated average treatment effect from the observed y in the treatment group
yc <- y[trt == 0]                # subset control group y

# seperate the dataset (dataint) to treatment group and control group
TRTmxi <- cbind(yt, xt)
CNTLmxi <- cbind(yc, xc)


# if the pairwise interaction of predictors are taken into consideration
if (inter == TRUE) {
    mt <- lm(yt ~., data = TRTmxi)
    mc <- lm(yc ~., data = CNTLmxi)
    et <- coefficients(mt)                                                        # extract coefficients from the treated model
    ec <- coefficients(mc)                                                        # extract coefficients from the control model
    pred_t <-  c(et[1] + (as.matrix(CNTLmxi[,2:(ncol(CNTLmxi))])%*%as.matrix(et[2:(length(et))])), et[1] + (as.matrix(TRTmxi[,2:(ncol(TRTmxi))])%*%as.matrix(et[2:(length(et))])))  # the predicted value to the treatment condition
    pred_c <-  c(ec[1] + (as.matrix(CNTLmxi[,2:(ncol(CNTLmxi))])%*%as.matrix(ec[2:(length(ec))])), ec[1] + (as.matrix(TRTmxi[,2:(ncol(TRTmxi))])%*%as.matrix(ec[2:(length(ec))])))  # the predicted value to the control condition
    
}

# calcualte the estimate: individual pite estimate and the standard deviation of pite estimates among individuals
pite.ind.int <- pred_t - pred_c
pite.sd.int <- sd(pite.ind.int) ; pite.sd.int


####   2-way interaction: Permutation test

# extract the pite result and the information from the dataset
n <- length(pite.ind.int)
#pite.ind <- pite.ind
#pite.sd <- pite.sd.int
TOTNumCov <- ncol(dataint2) -1

# extract y
y <- y

# extract x  
x <- dataint2[, !names(dataint2)%in%c("(Intercept)")]


# permute nperm times
nperm = 1000
prmt.int <- matrix(NA, nrow = n, ncol = nperm)
warn.count = 0

set.seed(9990)
for(i in 1:nperm){
    trt <- sample(trt, replace = F)    # shuffle (permute) the treatment conditions
    xtp <- as.matrix(x[trt == 1, ])                       # get treatment groups' covariates after shuffling
    xcp <- as.matrix(x[trt == 0, ])                       # get control groups' covariates after shuffling
    ytp <- y[trt == 1]                                    # get the y of the treatment
    ycp <- y[trt == 0]                                    # get the y of the control 
    mt <- lm(ytp ~ xtp)                                   # regress the observed y on the shuffled x of the treatment group 
    mc <- lm(ycp ~ xcp)                                   # regress the observed y on the suffled x of the control group 
    et <- coefficients(mt)                                # extract the coefficients from the regression model of the treatment group 
    ec <- coefficients(mc)                                # extract the coefficients from the regression model of the control group 
    
    # if singular condition happened in the dataset, the effects of covariates can be zero
    # record how many times sigular happened, if sigular happened, the loop will skip this permutation 
    # ==========  this is NOT an ideal way to extract sigular issue. 
    # ==========  It might be better to check sigular issue righ after the treatment were shuffled, rather than after the estimation
    if(sum(is.na(et))>0 | sum(is.na(ec))>0){               
        warn.count = warn.count + 1
        next
    }
    
    # calculate PITE for nperm times
    pred_t <-  c(et[1] + (xcp%*%et[2:(length(et))]), et[1] + (xtp%*%et[2:(length(et))]))  # calculate individual's predicted value to treatment
    pred_c <-  c(ec[1] + (xcp%*%ec[2:(length(ec))]), ec[1] + (xtp%*%ec[2:(length(ec))]))  # calcualte individual's predicted valued to control
    prmt.int[,i] <- pred_t-pred_c                                                             # calculate PITE for the ith permutation replicate
}

# calculate the standard deviatoin of PITEs for each permutation
prmt.sd.int <- apply(prmt.int, 2, sd, na.rm = T) 
hist(prmt.sd.int, xlim = c(0, 0.5), ylim = c(0, 80), breaks = 50, main = "Histogram of prmt.sd \n with 2-way interaction")
abline(v = pite.sd.int, col = 'red')
text(pite.sd.int, 10, paste("pite.sd = ", round(pite.sd.int, 3)), col = 'blue')
text(0.3, 30, paste("# of permutation = ", nperm))
text(0.3, 35, paste("# of singularity occurred = ", warn.count))




