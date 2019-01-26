## import the data
library(readr)
###qualdat<- read_csv("E:/whole traits/alltraits1.csv")
qualdat <- read.csv("~/Documents/whole traits/Copy of Copy of alltraits1.csv")
str(qualdat)
### remove the sur with "1" 
qualdat <- subset(qualdat, !(qualdat$Sur==1))

#######changing class of the date

qualdat$Headingt= as.Date(qualdat$headingt,format = "%m/%d/%Y")
qualdat$Flowert= as.Date(qualdat$flowert,format = "%m/%d/%Y")
qualdat$Halfht= as.Date(qualdat$halfht,format = "%m/%d/%Y")
qualdat$Hflowert= as.Date(qualdat$hflowert,format = "%m/%d/%Y")
qualdat$datest= as.Date(qualdat$datest,format = "%m/%d/%Y")

#####changing some date
###heading time
qualdat$Headingt[qualdat$Headingt=="2017-8-17"] <- as.Date("2017-8-16")
qualdat$Headingt[qualdat$Headingt=="2017-8-19"] <- as.Date("2017-8-18")
qualdat$Headingt[qualdat$Headingt=="2017-9-10"] <- as.Date("2017-9-11")

###Flower time
qualdat$Flowert[qualdat$Flowert=="2017-9-10"] <- as.Date("2017-9-11")
qualdat$Flowert[qualdat$Flowert=="2017-9-16"] <- as.Date("2017-9-15")


###half heading time
qualdat$Halfht[qualdat$Halfht=="2017-7-29"] <- as.Date("2017-7-28")
qualdat$Halfht[qualdat$Halfht=="2017-8-17"] <- as.Date("2017-8-16")
qualdat$Halfht[qualdat$Halfht=="2017-8-19"] <- as.Date("2017-8-18")
qualdat$Halfht[qualdat$Halfht=="2017-9-16"] <- as.Date("2017-9-15")


#####calculate days
qualdat$hday <- as.numeric(as.Date(qualdat$headingt,format = "%m/%d/%Y")-as.Date(qualdat$datest,format = "%m/%d/%Y"))
qualdat$fday <-as.numeric(as.Date(qualdat$flowert,format = "%m/%d/%Y")-as.Date(qualdat$datest,format = "%m/%d/%Y"))
qualdat$hhday <- as.numeric(as.Date(qualdat$halfht,format = "%m/%d/%Y")-as.Date(qualdat$datest,format = "%m/%d/%Y"))
qualdat$hfday <- as.numeric(as.Date(qualdat$hflowert,format = "%m/%d/%Y")-as.Date(qualdat$datest,format = "%m/%d/%Y"))

###write the data out 
write.csv(qualdat, file = "~/Documents/whole traits/Copy of alltraitsflowerday.csv",row.names = T, na = ".")


#### Scatterplot Matrices ( can hep to check some outlies)
pairs(~hday+fday+hhday+hfday,data=qualdat, 
      main="Simple Scatterplot Matrix for flowering traits")

## Examine distribution of brix data
par(mfrow=c(2,2))
boxplot(hday~ Rep, data=qualdat, xlab="Rep", ylab="heading days", main="boxplot of heading time in each rep", col="pink")

boxplot(fday~ Rep, data=qualdat, xlab="Rep", ylab="flowering days", main="boxplot of flowering time in each Rep", col="pink")

boxplot(hhday~ Rep, data=qualdat, xlab="Rep", ylab="half heading days", main="boxplot of half heading time in each Rep", col="pink")

boxplot(hfday~ Rep, data=qualdat, xlab="Rep", ylab="half flowering days", main="boxplot of half flowering time in each Rep", col="pink")

### checking normality 
###qq plot
par(mfrow=c(2,2))
qqnorm(qualdat$hday)
qqnorm(qualdat$fday)
qqnorm(qualdat$hhday)
qqnorm(qualdat$hfday)
###histogram
hist(qualdat$hday)
hist(qualdat$fday)
hist(qualdat$hhday)
hist(qualdat$hfday)

####changing class of variables
Entry= as.factor(qualdat$Entry)
REP = as.factor(qualdat$Rep)

## Calculate variance components
# requires lme4 package
install.packages("lme4")
install.packages("Matrix")
library(lme4)
library(Matrix)

##########################   Heading time days####################
##########################   Heading time days####################
##########################   Heading time days####################
####
####Heritibility

Misvarcomp <- lmer(hday~ (1|Entry)+ (1|Rep) + (1|Year), data=qualdat)
#####Extract variance components
summary(Misvarcomp)
###calulate the heritibility
####get all of the variance (variance in additive effect of each haploid genome)
Misvarcomp_Entry= unlist(VarCorr(Misvarcomp))[[1]];
Misvarcomp_Rep= unlist(VarCorr(Misvarcomp))[[2]];
Misvarcomp_Year= unlist(VarCorr(Misvarcomp))[[3]];
Var_Residual <- attr(VarCorr(Misvarcomp), "sc")^2
####heritability. I used 8 as the number of Var_residual devided, but i am not sure, please let me know if it is not right. Thanks.
h2 =(Misvarcomp_Entry)/((Misvarcomp_Entry)+(Var_Residual)/8) #####  
h2
###BLUP
## BLUPS
# fit the model
Mismodel <- lmer(hday~ (1|Entry)+ (1|Rep) + (1|Year), data=qualdat)
# estimate BLUPS
Misblup = ranef(Mismodel)
# look at output structure
str(Misblup)
# extract blup for line
Mislineblup = Misblup$Entry
# see the structure of the blup for each line
str(Mislineblup)
# save the brixlineblup output to a separate .csv file
write.csv(Mislineblup, file="~/Documents/whole traits/hdaysBLUPS.csv")
## Creating plots with the BLUPs
# Create a numeric vector with the BLUP for each line
LINEBLUP = Mislineblup[,1]

# Create a histogram with the BLUP for each line
par(mfcol=c(1,2))
hist(LINEBLUP, col="brown",main="histogram of Line BLUP of heading time")

qqnorm(LINEBLUP,main="Q-Q of Line blUP of heading time")
qqline(LINEBLUP)
## Compare BLUP to line averages on a scatterplot
lmean = tapply(qualdat$hday, qualdat$Entry, na.rm=T, mean,data=qualdat)
plot(LINEBLUP, lmean[row.names(Mislineblup)], col="blue",main="Line mean VS LineBlup",ylab="line mean") ### subsetted lmean to match LINEBLUP

#2 : Normality
mod1<-lmer(hday~ (1|Entry)+ (1|Rep) + (1|Year), data=qualdat)
par(mfcol=c(1,2))
hist(residuals(mod1), col="brown",main = "Residues distribution", xlab = "Residues", 
     ylab = "Frequency")

qqnorm(residuals(mod1),main ="Q-Q plot of residuals of heading time")
qqline(residuals(mod1))
plot(fitted(mod1), 
     residuals(mod1))

###########Flower time days############################
###########Flower time days############################
###########Flower time days############################
Misvarcomp <- lmer(fday~ (1|Entry)+ (1|Rep) + (1|Year), data=qualdat)
#####Extract variance components
summary(Misvarcomp)

###calulate the heritibility
####get all of the variance (variance in additive effect of each haploid genome)

Misvarcomp_Entry= unlist(VarCorr(Misvarcomp))[[1]];
Misvarcomp_Rep= unlist(VarCorr(Misvarcomp))[[2]];
Misvarcomp_Year= unlist(VarCorr(Misvarcomp))[[3]];
Var_Residual <- attr(VarCorr(Misvarcomp), "sc")^2
####heritability. I used 8 as the number of Var_residual devided, but i am not sure, please let me know if it is not right. Thanks.
h2 =(Misvarcomp_Entry)/((Misvarcomp_Entry)+(Var_Residual)/8) #####  
h2

###BLUP
## BLUPS
# fit the model
Mismodel <- lmer(fday~ (1|Entry)+ (1|Rep) + (1|Year), data=qualdat)
# estimate BLUPS
Misblup = ranef(Mismodel)
# look at output structure
str(Misblup)
# extract blup for line
Mislineblup = Misblup$Entry
# see the structure of the blup for each line
str(Mislineblup)
# save the brixlineblup output to a separate .csv file
write.csv(Mislineblup, file="~/Documents/whole traits/fdaysBLUPS.csv")
## Creating plots with the BLUPs
# Create a numeric vector with the BLUP for each line
LINEBLUP = Mislineblup[,1]
# Create a histogram with the BLUP for each line
par(mfcol=c(1,2))
hist(LINEBLUP, col="brown",main="histogram of Line BLUP of heading time")

qqnorm(LINEBLUP,main="Q-Q of Line blUP of heading time")
qqline(LINEBLUP)
## Compare BLUP to line averages on a scatterplot
lmean = tapply(qualdat$fday, qualdat$Entry, na.rm=T, mean,data=qualdat)
plot(LINEBLUP, lmean[row.names(Mislineblup)], col="blue",main="Line mean VS LineBlup",ylab="line mean") ### subsetted lmean to match LINEBLUP

#2 : Normality
mod1<-lmer(fday~ (1|Entry)+ (1|Rep) + (1|Year), data=qualdat)
par(mfcol=c(1,2))
hist(residuals(mod1), col="brown",main = "Residues distribution", xlab = "Residues", 
     ylab = "Frequency")

qqnorm(residuals(mod1),main ="Q-Q plot of residuals of heading time")
qqline(residuals(mod1))
plot(fitted(mod1), 
     residuals(mod1))


##########################   Half heading time days#######################################
##########################   Half heading time days#######################################
##########################   Half heading time days#######################################
####Heritibility

Misvarcomp <- lmer(hhday~ (1|Entry)+ (1|Rep) + (1|Year), data=qualdat)

#####Extract variance components
summary(Misvarcomp)

###calulate the heritibility
####get all of the variance (variance in additive effect of each haploid genome)

Misvarcomp_Entry= unlist(VarCorr(Misvarcomp))[[1]];
Misvarcomp_Rep= unlist(VarCorr(Misvarcomp))[[2]];
Misvarcomp_Year= unlist(VarCorr(Misvarcomp))[[3]];
Var_Residual <- attr(VarCorr(Misvarcomp), "sc")^2
####heritability. I used 8 as the number of Var_residual devided, but i am not sure, please let me know if it is not right. Thanks.
h2 =(Misvarcomp_Entry)/((Misvarcomp_Entry)+(Var_Residual)/8) #####  
h2


###BLUP
## BLUPS
# fit the model
Mismodel <- lmer(hhday~ (1|Entry)+ (1|Rep) + (1|Year), data=qualdat)
# estimate BLUPS
Misblup = ranef(Mismodel)
# look at output structure
str(Misblup)
# extract blup for line
Mislineblup = Misblup$Entry
# see the structure of the blup for each line
str(Mislineblup)
# save the brixlineblup output to a separate .csv file
write.csv(Mislineblup, file="~/Documents/whole traits/hhdaysBLUPS.csv")
## Creating plots with the BLUPs
# Create a numeric vector with the BLUP for each line
LINEBLUP = Mislineblup[,1]

# Create a histogram with the BLUP for each line
par(mfcol=c(1,2))
hist(LINEBLUP, col="brown",main="histogram of Line BLUP of heading time")

qqnorm(LINEBLUP,main="Q-Q of Line blUP of heading time")
qqline(LINEBLUP)
## Compare BLUP to line averages on a scatterplot
lmean = tapply(qualdat$hhday, qualdat$Entry, na.rm=T, mean,data=qualdat)
plot(LINEBLUP, lmean[row.names(Mislineblup)], col="blue",main="Line mean VS LineBlup",ylab="line mean") ### subsetted lmean to match LINEBLUP

#2 : Normality
mod1<-lmer(hhday~ (1|Entry)+ (1|Rep) + (1|Year), data=qualdat)
par(mfcol=c(1,2))
hist(residuals(mod1), col="brown",main = "Residues distribution", xlab = "Residues", 
     ylab = "Frequency")

qqnorm(residuals(mod1),main ="Q-Q plot of residuals of heading time")
qqline(residuals(mod1))
plot(fitted(mod1), 
     residuals(mod1))


##########################   Half flowering time days##########################################

##########################   Half flowering time days##########################################

##########################   Half flowering time days##########################################

Misvarcomp <- lmer(hfday~ (1|Entry)+ (1|Rep) + (1|Year), data=qualdat)

#####Extract variance components
summary(Misvarcomp)

###calulate the heritibility
####get all of the variance (variance in additive effect of each haploid genome)

Misvarcomp_Entry= unlist(VarCorr(Misvarcomp))[[1]];
Misvarcomp_Rep= unlist(VarCorr(Misvarcomp))[[2]];
Misvarcomp_Year= unlist(VarCorr(Misvarcomp))[[3]];
Var_Residual <- attr(VarCorr(Misvarcomp), "sc")^2
####heritability. I used 8 as the number of Var_residual devided, but i am not sure, please let me know if it is not right. Thanks.
h2 =(Misvarcomp_Entry)/((Misvarcomp_Entry)+(Var_Residual)/8) #####  
h2


###BLUP
## BLUPS
# fit the model
Mismodel <- lmer(hfday~ (1|Entry)+ (1|Rep) + (1|Year), data=qualdat)
# estimate BLUPS
Misblup = ranef(Mismodel)
# look at output structure
str(Misblup)
# extract blup for line
Mislineblup = Misblup$Entry
# see the structure of the blup for each line
str(Mislineblup)
# save the brixlineblup output to a separate .csv file
write.csv(Mislineblup, file="~/Documents/whole traits/hfdaysBLUPS.csv")
## Creating plots with the BLUPs
# Create a numeric vector with the BLUP for each line
LINEBLUP = Mislineblup[,1]

# Create a histogram with the BLUP for each line
par(mfcol=c(1,2))
hist(LINEBLUP, col="brown",main="histogram of Line BLUP of heading time")

qqnorm(LINEBLUP,main="Q-Q of Line blUP of heading time")
qqline(LINEBLUP)
## Compare BLUP to line averages on a scatterplot
lmean = tapply(qualdat$hfday, qualdat$Entry, na.rm=T, mean,data=qualdat)
plot(LINEBLUP, lmean[row.names(Mislineblup)], col="blue",main="Line mean VS LineBlup",ylab="line mean") ### subsetted lmean to match LINEBLUP

#2 : Normality
mod1<-lmer(hfday~ (1|Entry)+ (1|Rep) + (1|Year), data=qualdat)
par(mfcol=c(1,2))
hist(residuals(mod1), col="brown",main = "Residues distribution", xlab = "Residues", 
     ylab = "Frequency")

qqnorm(residuals(mod1),main ="Q-Q plot of residuals of heading time")
qqline(residuals(mod1))
plot(fitted(mod1), 
     residuals(mod1))

##########PCA analysis with ppca######################
##########PCA analysis with ppca######################
##########PCA analysis with ppca######################
#######
#######
#### 1: combining four blup files together######
####import the data
#######hday
ranhday <- read.csv("~/Documents/whole traits/hdaysBLUPS.csv",header=T,na.strings=c("","NA"))
ranhday$ID <- ranhday$X
ranhday$hd <- ranhday$X.Intercept.
str(ranhday)
ranhday <- ranhday[,c(3:4)]
str(ranhday)
#########fday
ranfday <- read.csv("~/Documents/whole traits/fdaysBLUPS.csv",header=T,na.strings=c("","NA"))
ranfday$ID <- ranfday$X
ranfday$fd <- ranfday$ X.Intercept.
str(ranfday)
ranfday <- ranfday[,c(3:4)]
str(ranfday)
########hhday
ranhhday <- read.csv("~/Documents/whole traits/hhdaysBLUPS.csv",header=T,na.strings=c("","NA"))
ranhhday$ID <- ranhhday$X
ranhhday$hhd <- ranhhday$X.Intercept.
str(ranhhday)
ranhhday <- ranhhday[,c(3:4)]
str(ranhhday)
########hfday
ranhfday <- read.csv("~/Documents/whole traits/hfdaysBLUPS.csv",header=T,na.strings=c("","NA"))
ranhfday$ID <- ranhfday$X
ranhfday$hf <- ranhfday$X.Intercept.
str(ranhfday)
ranhfday <- ranhfday[,c(3:4)]
str(ranhfday)

######merge four data sets ###############################
######merge four data sets ###############################
######merge four data sets ###############################
install.packages("dplyr")
library(dplyr)

########combine with missing values 
####
ftdaysblup <- join_all(list(ranfday,ranhday,ranhfday,ranhhday), by='ID')
str(ftdaysblup)
write.csv(ftdaysblup, file = "~/Documents/whole traits/hfhhfdaysblupm.csv",row.names = T)

##########PCA analysis with ppca
###########
source("https://bioconductor.org/biocLite.R") ###intall.packages("pcaMethods")
biocLite("pcaMethods")
library(pcaMethods)
#citation("pcaMethods")
pc <- pca(ftdaysblup, nPcs=1, method="ppca",center = TRUE)
fblupimputed <- data.frame(completeObs(pc))
write.csv(fblupimputed, file = "~/Documents/whole traits/fmissblupimputed.csv",row.names = T)
fblupPC1 <- scores(pc)
write.csv(fblupPC1, file = "~/Documents/whole traits/fblupPC1.csv",row.names = T)




