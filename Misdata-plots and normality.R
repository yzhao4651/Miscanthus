####import the data####

library(readr)
# set up working directory
mywd <- "~/Documents/whole traits" # for Yongli
# mywd <- "." # for Lindsay

qualdat <- read.csv(file.path(mywd,"Copy of Copy of alltraits1.csv"),na.strings = c("",".","NA"))
#check the data formate
str(qualdat)
####remove the sur with "1" 
qualdat <- subset(qualdat, !(qualdat$Sur==1))
####calculate days
qualdat$hday <- as.numeric(as.Date(qualdat$headingt,format = "%m/%d/%Y")-as.Date(qualdat$datest,format = "%m/%d/%Y"))
qualdat$fday <-as.numeric(as.Date(qualdat$flowert,format = "%m/%d/%Y")-as.Date(qualdat$datest,format = "%m/%d/%Y"))
qualdat$hhday <- as.numeric(as.Date(qualdat$halfht,format = "%m/%d/%Y")-as.Date(qualdat$datest,format = "%m/%d/%Y"))
qualdat$hfday <- as.numeric(as.Date(qualdat$hflowert,format = "%m/%d/%Y")-as.Date(qualdat$datest,format = "%m/%d/%Y"))
#check the data formate
str(qualdat)
####Scatterplot Matrices ( can hep to check some outlies)
pairs(~hday+fday+hhday+hfday,data=qualdat, 
      main="Simple Scatterplot Matrix for flowering traits")
####Examine distribution of brix data
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
##
####Scatterplot Matrices ( can hep to check some outlies)
pairs(~FW_kg+length_cm+Outdi_cm+Indi_cm+nodes+BC_cm+FN,data=qualdat, 
      main="Simple Scatterplot Matrix for all traits")
####Examine distribution of brix data
par(mfrow=c(2,4))
boxplot(FW_kg~ Rep, data=qualdat, xlab="Rep", ylab="Fresh weight", main="boxplot of Fresh weight in each rep", col="pink")
boxplot(length_cm~ Rep, data=qualdat, xlab="Rep", ylab="Length", main="boxplot of Length in each Rep", col="pink")
boxplot(Outdi_cm~ Rep, data=qualdat, xlab="Rep", ylab="Out diameter", main="boxplot of diameter at Basal Internode in each Rep", col="pink")
boxplot(Indi_cm~ Rep, data=qualdat, xlab="Rep", ylab="", main="boxplot of diameter at Last Internode in each Rep", col="pink")
boxplot(nodes~ Rep, data=qualdat, xlab="Rep", ylab="nodes", main="boxplot of nodes in each Rep", col="pink")
boxplot(BC_cm~ Rep, data=qualdat, xlab="Rep", ylab="Basal Circumference", main="boxplot of Basal Circumference in each Rep", col="pink")
boxplot(FN~ Rep, data=qualdat, xlab="Rep", ylab="flower number", main="boxplot of flower number in each Rep", col="pink")
####checking normality ####
###qq plot
par(mfrow=c(2,4))
qqnorm(qualdat$FW_kg)
qqnorm(qualdat$length_cm)
qqnorm(qualdat$Outdi_cm)
qqnorm(qualdat$Indi_cm)
qqnorm(qualdat$nodes)
qqnorm(qualdat$BC_cm)
qqnorm(qualdat$FN)
###histogram
par(mfrow=c(2,5))
hist(qualdat$FW_kg)
hist(qualdat$length_cm)
hist(qualdat$Outdi_cm)
hist(qualdat$Indi_cm)
hist(qualdat$nodes)
hist(qualdat$BC_cm)
hist(qualdat$FN)

#### fresh weight and Flower number need data transform
## Third way. 
qualdatFW<-subset(qualdat, !is.na(qualdat$FW_kg))
qualdatFN<-subset(qualdat, !is.na(qualdat$FN))
qualdatfday<-subset(qualdat, !is.na(qualdat$fday))
qualdathfday<-subset(qualdat, !is.na(qualdat$hfday))
source("~/Documents/R-corde for miscanthus project/bcplot function.txt")
par(mfrow=c(2,3))
bcplot(qualdatFW$FW_kg) ####-0.025
bcplot(qualdatFN$FN)####0.175
bcplot(qualdatfday$fday) ####-0.1
bcplot(qualdathfday$hfday)####0.15

# Lindsay's code to get the transformed data
bc <- function(x, lda){ # function to transform data after lambda is determined
  while(min(x, na.rm = TRUE) <= 0){
    x <- x + 1 # positive numbers required
  }
  if(lda == 0){
    x1 <- log(x)
  } else {
    x1 <- (x^lda - 1)/lda
  }
  return(x1)
}
qualdat$FW_kg<- bc(qualdat$FW_kg,lda=-0.025)
qualdat$FN <- bc(qualdat$FN,lda=0.175)
qualdat$fday<- bc(qualdat$fday,lda=-0.1)
qualdat$hfday <- bc(qualdat$hfday,lda=0.15)
str(qualdat)
####write the data out 
# mywd <- "~/Yongli" # for Lindsay
write.csv(qualdat, file = file.path(mywd, "alltraitsnormalited.csv"),row.names = T, na = ".")
