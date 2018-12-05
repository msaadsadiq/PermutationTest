
# all ALS datasets
adv_Events = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\AdverseEvents.csv", header = T)
alsfrs = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\alsfrs.csv", header = T)
alsHist = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\AlsHistory.csv", header = T)
ConMeds = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\ConMeds.csv", header = T)
death = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\DeathData.csv", header = T)
dmg = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\demographics.csv", header = T)
FamHist = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\FamilyHistory.csv", header = T)
Fvc = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\Fvc.csv", header = T)
labs = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\Labs.csv", header = T)
Rilu = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\Riluzole.csv", header = T)
Svc = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\Svc.csv", header = T)
TRT = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\Treatment.csv", header = T)
VitSign = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\VitalSigns.csv", header = T)


## Dependent Varaible: datasets in PRO-ACT data folder (Var:subjectID, ALSFRS_slope)
ldbrd = read.table("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\PRO-ACT data\\PROACT_ progression\\ALSFRS_slope_PROACT_leaderboard.txt", header = T, sep = "|")
train = read.table("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\PRO-ACT data\\PROACT_ progression\\ALSFRS_slope_PROACT_training.txt", header = T, sep = "|")
train2 = read.table("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\PRO-ACT data\\PROACT_ progression\\ALSFRS_slope_PROACT_training2.txt", header = T, sep = "|")
valid = read.table("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\10172017\\Dataset\\PRO-ACT data\\PROACT_ progression\\ALSFRS_slope_PROACT_validation.txt", header = T, sep = "|")
total = rbind(ldbrd, train, train2, valid)
names(total)[1] = "subject_id"

# Dataset from Jack
alsJack = read.csv("C:\\Users\\Dawn\\Dropbox\\PITEPermutation\\ApplicationALS\\ALS.csv", header = T)



#============ Clean VitSign Dataset
# remove subjects' measures at time > 0 ( only keep the measured time at time 0 and time negative)
# sort the subjects' measures, keep the one's time that is the largest, if they have multiple time points measured 0 (trial start) and before
# 
VitSign.B = VitSign[which(VitSign$Vital_Signs_Delta <= 0),] 
VitSign.BO = VitSign.B[order(VitSign.B$Vital_Signs_Delta, decreasing = T),]
VitSign.NonDup = VitSign.BO[!duplicated(VitSign.BO$subject_id),]



# clean weight #

w = VitSign.NonDup[!is.na(VitSign.NonDup$Weight),names(VitSign.NonDup) %in% c('subject_id', 'Weight', "Weight_Units")]
wp = w[which(w$Weight_Units == 'Pounds'& !is.na(w$Weight)) ,]
wk = w[which(w$Weight_Units == 'Kilograms'& !is.na(w$Weight)),]
wn = w[which(w$Weight_Units == ''& !is.na(w$Weight)),]


library(sm)
w.u <- factor(w$Weight_Units, levels = c('', 'Kilograms', 'Pounds'), 
              labels = c("None", "KG", "Pounds"))
sm.density.compare(w$Weight, w.u, xlab = "Weight")
colfill<-c(2:(2+length(levels(w.u)))) 
title(main = "Density Plots of Different Weight Units")
legend(locator(1), levels(w.u), fill=colfill)

# recode all having weight values but with unit ==  '' to 'Kilograms' (careful!)
VitSign.NonDup[!is.na(VitSign.NonDup$Weight) & VitSign.NonDup$Weight_Units =='',]$Weight_Units <- 'Kilograms'

# Use Baseline_weight as the priority
VitSign.NonDup$Weight.kg <- VitSign.NonDup$Baseline_Weight

# convert pounds to kg
VitSign.NonDup$Weight.kg[is.na(VitSign.NonDup$Weight.kg)] = VitSign.NonDup$Weight[is.na(VitSign.NonDup$Weight.kg)] 
VitSign.NonDup$Weight.kg[VitSign.NonDup$Weight_Units == 'Pounds'] <- VitSign.NonDup$Weight[VitSign.NonDup$Weight_Units == 'Pounds']*0.453592



# clean height # 

h = VitSign.NonDup[!is.na(VitSign.NonDup$Height),names(VitSign.NonDup) %in% c('subject_id', 'Height', "Height_Units")]
hcm = h[which(h$Height_Units == 'Centimeters'& !is.na(h$Height)) ,]
hin = h[which(h$Height_Units == 'Inches'& !is.na(h$Height)),]
hn = h[which(h$Height_Units == ''& !is.na(h$Height)),]

h.u <- factor(h$Height_Units, levels = c('', 'Centimeters', 'Inches'), 
              labels = c("None", "cm", "inches"))
sm.density.compare(h$Height, h.u, xlab = "Height")
colfill<-c(2:(2+length(levels(h.u)))) 
title(main = "Density Plots of Different Height Units")
legend(locator(1), levels(h.u), fill=colfill)

# recode all subjects ahving height values but with unit == '' to "Centimeters' (careful!)
VitSign.NonDup[!is.na(VitSign.NonDup$Height) & VitSign.NonDup$Height_Units =='',]$Height_Units <- 'Centimeters'

# convert Inches to Centimeters
VitSign.NonDup$Height.cm = VitSign.NonDup$Height
VitSign.NonDup$Height.cm[VitSign.NonDup$Height_Units == 'Inches'] <- VitSign.NonDup$Height[VitSign.NonDup$Height_Units == 'Inches']*2.54



# clean temparature #
t = VitSign.NonDup[!is.na(VitSign.NonDup$Temperature),names(VitSign.NonDup) %in% c('subject_id', 'Temperature', "Temperature_Units")]
tn = t[which(t$Temperature_Units == ''& !is.na(t$Temperature)),]
t.u <- factor(t$Temperature_Units, levels = c('', 'C'), 
              labels = c("None", "Celsius"))
sm.density.compare(t$Temperature, t.u, xlab = "Temperature")
colfill<-c(2:(2+length(levels(t.u)))) 
title(main = "Density Plots of Different Temperature Units")
legend(locator(1), levels(t.u), fill=colfill)
# remove 2 subjects temparature < 25
VitSign.NonDup$Temperature[VitSign.NonDup$Temperature < 34.8] <- NA


# clean blood pressure 

# combine bp diastolic and supine diastolic 
VitSign.NonDup$BP_Diastolic = VitSign.NonDup$Blood_Pressure_Diastolic
VitSign.NonDup$BP_Diastolic[is.na(VitSign.NonDup$BP_Diastolic)] <- VitSign.NonDup$Baseline_Supine_BP_Diastolic[is.na(VitSign.NonDup$BP_Diastolic)]

VitSign.NonDup$BP_Systolic = VitSign.NonDup$Blood_Pressure_Systolic
VitSign.NonDup$BP_Systolic[is.na(VitSign.NonDup$BP_Systolic)] <- VitSign.NonDup$Baseline_Supine_BP_Systolic[is.na(VitSign.NonDup$BP_Systolic)]


# clean pulse 
VitSign.NonDup$Pulse_C = VitSign.NonDup$Baseline_Supine_Pulse # 7038 missing in Pulse_C
VitSign.NonDup$Pulse_C[is.na(VitSign.NonDup$Pulse_C)] <- VitSign.NonDup$Supine_Pulse[is.na(VitSign.NonDup$Pulse_C)] #6278 missing
VitSign.NonDup$Pulse_C[is.na(VitSign.NonDup$Pulse_C)] <- VitSign.NonDup$Pulse[is.na(VitSign.NonDup$Pulse_C)] #340 missing


# Select variables in VitSign.NonDup  
VS = VitSign.NonDup[,names(VitSign.NonDup) %in% c("subject_id", "BP_Diastolic", "BP_Systolic","Pulse_C","Respiratory_Rate","Temperature",
    "Pulse_C", "Height.cm","Weight.kg")]

#dim(VS[!is.na(VS$Baseline_Standing_BP_Diastolic) & !is.na(VS$Baseline_Supine_BP_Diastolic),])

#============ Clean treatment dataset 
# treatment group delta (the time duration between first time the patient was assessed during the trial and the time medication was first given))
# flag the treatment group delta if the duration is larger than 365

TRT$DeltaFlag = as.factor(TRT$'Treatment_Group_Delta' > 365)
TRT = TRT[,!names(TRT) == 'Treatment_Group_Delta']

#============ Clean Demographic dataset
dmg_C = dmg[,names(dmg) %in% c("subject_id","Demographics_Delta", "Age","Date_of_Birth", "Race_Caucasian", "Sex")]

# Sex
dmg_C$Gender = as.factor(ifelse(dmg_C$Sex=='Female', 'Female', ifelse(dmg_C$Sex =='Male', 'Male', NA)))

# AGE
dmg_C$Age_C = dmg_C$Age #3021 missing
dmg_C$Age_C[is.na(dmg_C$Age_C)] <- dmg_C$Date_of_Birth[is.na(dmg_C$Age_C)]/-365 # 2077 missing

# Race
dmg_C$White = ifelse(is.na(dmg_C$Race_Caucasian),0,1) # white: 6716, NonWhite: 4007 (in the nonwhite section, there might be a portion of missing)

# (duplication checked! not shown)
DMG = dmg_C[, names(dmg_C) %in% c("subject_id", "Age_C", "White", "Gender")]



#=========== Clean adv_Events dataset
adv_Events$severityBin = ifelse(adv_Events$Severity %in% c('Fatal', 'Life Threatening','Marked', 'Moderate','MODERATE', 'Severe',  'SEVERE', 'Yes'), 1, 
                                ifelse(adv_Events$Severity %in% c('Mild', 'MILD', 'No' ),0, NA))
ADV = adv_Events[!duplicated(adv_Events$subject_id), names(adv_Events) %in% c('subject_id', 'severityBin')]

#=========== Clean alsHist dataset
ALSHIST_C = alsHist[, names(alsHist) %in% c("subject_id", "Site_of_Onset" , "Diagnosis_Delta")]
ALSHIST= ALSHIST_C[!is.na(ALSHIST_C$Site_of_Onset) & ALSHIST_C$Site_of_Onset!='',]
# recode dummy 
ALSHIST$LimbOnly = ifelse(ALSHIST$Site_of_Onset=='Onset: Limb', 1, 0)
ALSHIST$BulbarOnly = ifelse(ALSHIST$Site_of_Onset == 'Onset: Bulbar', 1, 0)
ALSHIST = ALSHIST[,names(ALSHIST) %in% c('subject_id','LimbOnly', 'BulbarOnly','Diagnosis_Delta')]



#=========== Clean ConMeds dataset

# link ConMeds and ALSHIST
ConMeds$UseRiluzole  = ifelse(ConMeds$Medication_Coded=='RILUZOLE', 1, 0)

x = ConMeds[order(ConMeds$subject_id, ConMeds$Start_Delta),]
RUser = x[which(x$UseRiluzole ==1),]  # 2525 10
RUB = RUser[which(RUser$Start_Delta<=0), ]  # 2280 10
RUBU = RUB[!duplicated(RUB$subject_id),] # 2252 10 

RLHist = RUBU[, names(RUBU) %in% c("subject_id", "Start_Delta","UseRiluzole")]


CMlist = data.frame(table(ConMeds$subject_id))$Var1
RLlist = data.frame(table(RLHist$subject_id))$Var1
OTlist = CMlist[!CMlist %in% RLlist]
OTHist = data.frame(subject_id = OTlist, Start_Delta = 0, UseRiluzole = 0 )

CM = rbind(RLHist, OTHist)



# MERGE THEM!

test = Reduce(function(...) merge(..., by='subject_id', all=TRUE), list(total, TRT, VS, DMG, ADV, ALSHIST, CM))
NFinal = test[!is.na(test$ALSFRS_slope) & !is.na(test$Study_Arm),] #2910
summary(NFinal)
str(NFinal)


