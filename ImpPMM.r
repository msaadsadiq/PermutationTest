### predictve mean matching

#install.packages('mice', dependencies = T)
library(mice)
library(mice)
library(VIM)
library(lattice)
library(ggplot2)

set.seed(9999)
# md.pattern(NFinal)
# p = md.pairs(NFinal) 
# pbox(NFinal, pos = 2)
# pbox(NFinal, pos = 3)
# pbox(NFinal, pos = 4)
# pbox(NFinal, pos = 5)
# pbox(NFinal, pos = 6)
imp1 <- mice(NFinal, m = 5)
imp_tot2 <- complete(imp1, "long", inc = TRUE)
origin = imp_tot2[imp_tot2$.imp == 0,]

#dev.off()
col <- rep(c("blue", "red")[1 + as.numeric(is.na(imp1$data$Temperature))], 6)
stripplot(Temperature ~ .imp, data = imp_tot2, jit = TRUE, ylim = c(34.5, 39), col = col, xlab = "imputation Number")

col <- rep(c("blue", "red")[1 + as.numeric(is.na(imp1$data$Height.cm))], 6)
stripplot(Height.cm ~ .imp, data = imp_tot2, jit = TRUE, col = col, xlab = "imputation Number")

col <- rep(c("blue", "red")[1 + as.numeric(is.na(imp1$data$Diagnosis_Delta))], 6)
stripplot(Diagnosis_Delta ~ .imp, data = imp_tot2, jit = TRUE, col = col, xlab = "imputation Number")

col <- rep(c("blue", "red")[1 + as.numeric(is.na(imp1$data$DeltaFlag))], 6)
stripplot(DeltaFlag ~ .imp, data = imp_tot2, jit = TRUE, col = col, xlab = "imputation Number")

col <- rep(c("blue", "red")[1 + as.numeric(is.na(imp1$data$Start_Delta))], 6)
stripplot(Start_Delta ~ .imp, data = imp_tot2, jit = TRUE, col = col, xlab = "imputation Number")
