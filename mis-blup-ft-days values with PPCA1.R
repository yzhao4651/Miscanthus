####import the data####
library(readr)
qualdat <- read.csv("Copy of Copy of alltraits1.csv")

####remove the sur with "1" 

qualdat <- subset(qualdat, !(qualdat$Sur==1))

####calculate days
qualdat$hday <- as.numeric(as.Date(qualdat$headingt,format = "%m/%d/%Y")-as.Date(qualdat$datest,format = "%m/%d/%Y"))
qualdat$fday <-as.numeric(as.Date(qualdat$flowert,format = "%m/%d/%Y")-as.Date(qualdat$datest,format = "%m/%d/%Y"))
qualdat$hhday <- as.numeric(as.Date(qualdat$halfht,format = "%m/%d/%Y")-as.Date(qualdat$datest,format = "%m/%d/%Y"))
qualdat$hfday <- as.numeric(as.Date(qualdat$hflowert,format = "%m/%d/%Y")-as.Date(qualdat$datest,format = "%m/%d/%Y"))

####write the data out 
write.csv(qualdat, file = "~/Documents/whole traits/Copy of alltraitsflowerday.csv",row.names = T, na = ".")


####Scatterplot Matrices ( can hep to check some outlies)
pairs(~hday+fday+hhday+hfday,data=qualdat, 
      main="Simple Scatterplot Matrix for flowering traits")

####Examine distribution of flowering data
par(mfrow=c(2,2))
boxplot(hday~ Rep, data=qualdat, xlab="Rep", ylab="heading days", main="boxplot of heading time in each rep", col="pink")

boxplot(fday~ Rep, data=qualdat, xlab="Rep", ylab="flowering days", main="boxplot of flowering time in each Rep", col="pink")

boxplot(hhday~ Rep, data=qualdat, xlab="Rep", ylab="half heading days", main="boxplot of half heading time in each Rep", col="pink")

boxplot(hfday~ Rep, data=qualdat, xlab="Rep", ylab="half flowering days", main="boxplot of half flowering time in each Rep", col="pink")

####checking normality ####
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

# set up a working directory, to make it easy to go back and forth between our two computers
workingdir <- "~/Documents/whole traits"
# workingdir <- "C:/Users/lvclark/Documents/Yongli" # for Lindsay

####calculate BLUP for flowering traits
###requires lme4 package
install.packages("lme4")
install.packages("Matrix")
library(lme4)
library(Matrix)
source("hetability fucntion .R") # load heritability function
####changing class of variables
qualdat$Entry= as.factor(qualdat$Entry)
qualdat$Rep = as.factor(qualdat$Rep)
qualdat$Year = as.factor(qualdat$Year)

#### Heading time days#####
####Heritibility
h2 <- heritability(qualdat$hday, qualdat)
h2
####BLUP of flowering time

### fit the model
Mismodel <- lmer(hday~ (1|Entry)+ (1|Rep) + (1|Year) + (1|Entry:Year), data=qualdat)
### estimate BLUPS
Misblup = ranef(Mismodel)
### look at output structure
str(Misblup)
### extract blup for line
Mislineblup = Misblup$Entry
### see the structure of the blup for each line
str(Mislineblup)
### save the brixlineblup output to a separate .csv file
write.csv(Mislineblup, file= file.path(workingdir, "hdaysBLUPS.csv"))
### Creating plots with the BLUPs
### Create a numeric vector with the BLUP for each line
LINEBLUP = Mislineblup[,1]

### Create a histogram with the BLUP for each line
par(mfcol=c(1,2))
hist(LINEBLUP, col="brown",main="histogram of Line BLUP of heading time")

qqnorm(LINEBLUP,main="Q-Q of Line blUP of heading time")
qqline(LINEBLUP)
###Compare BLUP to line averages on a scatterplot
lmean = tapply(qualdat$hday, qualdat$Entry, na.rm=T, mean,data=qualdat)
plot(LINEBLUP, lmean[row.names(Mislineblup)], col="blue",main="Line mean VS LineBlup",ylab="line mean") ### subsetted lmean to match LINEBLUP

####Flower time days#####

####calulate the heritibility

h2 <- heritability(qualdat$fday, qualdat)
h2

###BLUP
## BLUPS
# fit the model
Mismodel <- lmer(fday~ (1|Entry)+ (1|Rep) + (1|Year) + (1|Entry:Year), data=qualdat)
# estimate BLUPS
Misblup = ranef(Mismodel)
# look at output structure
str(Misblup)
# extract blup for line
Mislineblup = Misblup$Entry
# see the structure of the blup for each line
str(Mislineblup)
# save the brixlineblup output to a separate .csv file
write.csv(Mislineblup, file= file.path(workingdir, "fdaysBLUPS.csv"))
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
mod1<-lmer(fday~ (1|Entry)+ (1|Rep) + (1|Year) + (1|Entry:Year), data=qualdat)
par(mfcol=c(1,2))
hist(residuals(mod1), col="brown",main = "Residues distribution", xlab = "Residues", 
     ylab = "Frequency")

qqnorm(residuals(mod1),main ="Q-Q plot of residuals of heading time")
qqline(residuals(mod1))
plot(fitted(mod1), 
     residuals(mod1))


####  Half heading time days####

#### Heritibility

h2 <- heritability(qualdat$hhday, qualdat)
h2


###BLUP
## BLUPS
# fit the model
Mismodel <- lmer(hhday~ (1|Entry)+ (1|Rep) + (1|Year) + (1|Entry:Year), data=qualdat)
# estimate BLUPS
Misblup = ranef(Mismodel)
# look at output structure
str(Misblup)
# extract blup for line
Mislineblup = Misblup$Entry
# see the structure of the blup for each line
str(Mislineblup)
# save the brixlineblup output to a separate .csv file
write.csv(Mislineblup, file= file.path(workingdir, "hhdaysBLUPS.csv"))
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
mod1<-lmer(hhday~ (1|Entry)+ (1|Rep) + (1|Year) + (1|Entry:Year), data=qualdat)
par(mfcol=c(1,2))
hist(residuals(mod1), col="brown",main = "Residues distribution", xlab = "Residues", 
     ylab = "Frequency")

qqnorm(residuals(mod1),main ="Q-Q plot of residuals of heading time")
qqline(residuals(mod1))
plot(fitted(mod1), 
     residuals(mod1))

#### Half flowering time days####

###calulate the heritibility

h2 <- heritability(qualdat$hfday, qualdat) 
h2

###BLUP
## BLUPS
# fit the model
Mismodel <- lmer(hfday~ (1|Entry)+ (1|Rep) + (1|Year) + (1|Entry:Year), data=qualdat)
# estimate BLUPS
Misblup = ranef(Mismodel)
# look at output structure
str(Misblup)
# extract blup for line
Mislineblup = Misblup$Entry
# see the structure of the blup for each line
str(Mislineblup)
# save the brixlineblup output to a separate .csv file
write.csv(Mislineblup, file= file.path(workingdir, "hfdaysBLUPS.csv"))
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
mod1<-lmer(hfday~ (1|Entry)+ (1|Rep) + (1|Year) + (1|Entry:Year), data=qualdat)
par(mfcol=c(1,2))
hist(residuals(mod1), col="brown",main = "Residues distribution", xlab = "Residues", 
     ylab = "Frequency")

qqnorm(residuals(mod1),main ="Q-Q plot of residuals of heading time")
qqline(residuals(mod1))
plot(fitted(mod1), 
     residuals(mod1))


#### PCA analysis with ppca####
####
####
#### 1: combining four blup files together######
####import the data
####hday
ranhday <- read.csv(file.path(workingdir, "hdaysBLUPS.csv"),header=T,na.strings=c("","NA"))
ranhday$ID <- ranhday$X
ranhday$hd <- ranhday$X.Intercept.
str(ranhday)
ranhday <- ranhday[,c(3:4)]
str(ranhday)
####fday
ranfday <- read.csv(file.path(workingdir, "fdaysBLUPS.csv"),header=T,na.strings=c("","NA"))
ranfday$ID <- ranfday$X
ranfday$fd <- ranfday$ X.Intercept.
str(ranfday)
ranfday <- ranfday[,c(3:4)]
str(ranfday)
####hhday
ranhhday <- read.csv(file.path(workingdir, "hhdaysBLUPS.csv"),header=T,na.strings=c("","NA"))
ranhhday$ID <- ranhhday$X
ranhhday$hhd <- ranhhday$X.Intercept.
str(ranhhday)
ranhhday <- ranhhday[,c(3:4)]
str(ranhhday)
####hfday
ranhfday <- read.csv(file.path(workingdir, "hfdaysBLUPS.csv"),header=T,na.strings=c("","NA"))
ranhfday$ID <- ranhfday$X
ranhfday$hf <- ranhfday$X.Intercept.
str(ranhfday)
ranhfday <- ranhfday[,c(3:4)]
str(ranhfday)

####merge four data sets ####
install.packages("dplyr")
library(dplyr)

####combine with missing values 
####
ftdaysblup <- join_all(list(ranfday,ranhday,ranhfday,ranhhday), by='ID')
str(ftdaysblup)
write.csv(ftdaysblup, file = file.path(workingdir, "hfhhfdaysblupm.csv"),row.names = T)

####PCA analysis with ppca
####
# install.packages("BiocManager")
BiocManager::install("pcaMethods") ###intall.packages("pcaMethods")
library(pcaMethods)
#citation("pcaMethods")
pc <- pca(ftdaysblup, nPcs=1, method="ppca",center = TRUE)
fblupimputed <- data.frame(completeObs(pc))
write.csv(fblupimputed, file = file.path(workingdir,"fmissblupimputed.csv"),row.names = T)
fblupPC1 <- scores(pc)
write.csv(fblupPC1, file = file.path(workingdir, "fblupPC1.csv"),row.names = T)




