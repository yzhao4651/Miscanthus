####import the data####
library(readr)
mywd <- "~/Documents/whole traits" # on Yongli's computer
# mywd <- "." # on Lindsay's computer
qualdat <- read.csv(file.path(mywd, "trait1718.csv"), na.strings = c("",".","NA"))
#check the data formate
str(qualdat)

###change character to number N=1,B=2,F=3,P=4
levels(qualdat$GS)[levels(qualdat$GS)=="N"] <- "1"
levels(qualdat$GS)[levels(qualdat$GS)=="n"] <- "1"
levels(qualdat$GS)[levels(qualdat$GS)=="B"] <- "2"
levels(qualdat$GS)[levels(qualdat$GS)=="b"] <- "2"
levels(qualdat$GS)[levels(qualdat$GS)=="F"] <- "3"
levels(qualdat$GS)[levels(qualdat$GS)=="p"] <- "4"
levels(qualdat$GS)[levels(qualdat$GS)=="P"] <- "4"

####calculate days
qualdat$SRD <- as.numeric(as.Date(qualdat$SRD,format = "%m/%d/%Y")-as.Date(qualdat$datest2,format = "%m/%d/%Y"))
qualdat$ADD <-as.numeric(as.Date(qualdat$ADD,format = "%m/%d/%Y")-as.Date(qualdat$datest2,format = "%m/%d/%Y"))
qualdat$HD_1 <- as.numeric(as.Date(qualdat$HD_1,format = "%m/%d/%Y")-as.Date(qualdat$datest1,format = "%m/%d/%Y"))
qualdat$FD_1 <-as.numeric(as.Date(qualdat$FD_1,format = "%m/%d/%Y")-as.Date(qualdat$datest1,format = "%m/%d/%Y"))
qualdat$HD_50. <- as.numeric(as.Date(qualdat$HD_50.,format = "%m/%d/%Y")-as.Date(qualdat$datest1,format = "%m/%d/%Y"))
qualdat$FD_50. <- as.numeric(as.Date(qualdat$FD_50.,format = "%m/%d/%Y")-as.Date(qualdat$datest1,format = "%m/%d/%Y"))
####check the data formate
str(qualdat)

#### change several variables numeric format
indx <- sapply(qualdat[,c(16,19:21,23)], is.integer)
qualdat[,c(16,19:21,23)][indx] <- lapply(qualdat[,c(16,19:21,23)][indx], function(x) as.numeric(as.character(x)))
qualdat$GS <- as.numeric(qualdat$GS)
####check the data formate
str(qualdat)

write.csv(qualdat,file="~/Documents/whole traits/data1718updated.csv",row.names = T, na = ".")
names(qualdat)

####Scatterplot Matrices for flower traits ( can hep to check some outlies)
pairs(~HD_1+FD_1+HD_50.+FD_50.,data=qualdat, 
      main="Simple Scatterplot Matrix for flowering traits")
####scater plots for other traits 
pairs(~Yld_kg+SDW_kg+Cml_cm+CmD_BI_mm+CmD_LI_mm+CmN.+CmDW_g+Bcirc_cm+FNMain+FNsmall+TFN.+FD+CCirc_cm+GS+Lg,data=qualdat, 
      main="Simple Scatterplot Matrix for flowering traits")

####Examine distribution of brix data
par(mfrow=c(2,2))
boxplot(HD_1~ Rep, data=qualdat, xlab="Rep", ylab="heading days", main="boxplot of heading time in each rep", col="pink")
boxplot(FD_1~ Rep, data=qualdat, xlab="Rep", ylab="flowering days", main="boxplot of flowering time in each Rep", col="pink")
boxplot(HD_50.~ Rep, data=qualdat, xlab="Rep", ylab="half heading days", main="boxplot of half heading time in each Rep", col="pink")
boxplot(FD_50.~ Rep, data=qualdat, xlab="Rep", ylab="half flowering days", main="boxplot of half flowering time in each Rep", col="pink")
####checking normality ####
###qq plot
par(mfrow=c(2,2))
qqnorm(qualdat$HD_1)
qqnorm(qualdat$FD_1)
qqnorm(qualdat$HD_50.)
qqnorm(qualdat$FD_50.)
###histogram
par(mfrow=c(2,2))
hist(qualdat$HD_1)
hist(qualdat$FD_1)
hist(qualdat$HD_50.)
hist(qualdat$FD_50.)
names(qualdat)
####scater plots for other traits 
pairs(~Yld_kg+SDW_kg+Cml_cm+CmD_BI_mm+CmD_LI_mm+CmN.+CmDW_g+Bcirc_cm+FNMain+FNsmall+TFN.+FD+CCirc_cm+GS+Lg+SRD+ADD,data=qualdat, 
      main="Simple Scatterplot Matrix for flowering traits")

####Examine distribution of brix data
par(mfrow=c(4,5))
boxplot(Yld_kg~ Rep, data=qualdat, xlab="Rep", ylab="Length", main="boxplot of Biomass Yeild in each Rep", col="pink")
boxplot(SDW_kg~ Rep, data=qualdat, xlab="Rep", ylab="Fresh weight", main="boxplot of sample dry weight  in each rep", col="pink")
boxplot(Cml_cm~ Rep, data=qualdat, xlab="Rep", ylab="Length", main="boxplot of Culm Length in each Rep", col="pink")
boxplot(CmD_BI_mm~ Rep, data=qualdat, xlab="Rep", ylab="Culm Outer Diameter at Basal Internode", main="boxplot of Culm Outer Diameter at Basal Internode in each Rep", col="pink")
boxplot(CmD_LI_mm~ Rep, data=qualdat, xlab="Rep", ylab="Culm Outer Diameter at Last Internode", main="boxplot of Culm Outer Diameter at Last Internode in each Rep", col="pink")
boxplot(CmN.~ Rep, data=qualdat, xlab="Rep", ylab="Culm Node Number", main="boxplot of Culm Node Number in each Rep", col="pink")
boxplot(CmDW_g~ Rep, data=qualdat, xlab="Rep", ylab=" Culm Dry Weight ", main="boxplot of Culm Dry Weight in each Rep", col="pink")
boxplot(Bcirc_cm~ Rep, data=qualdat, xlab="Rep", ylab="Basal Circumference", main="boxplot of Basal Circumference in each Rep", col="pink")
boxplot(FNMain~ Rep, data=qualdat, xlab="Rep", ylab="Compressed PLant Circumference", main="boxplot of Main flower number in each Rep", col="pink")
boxplot(FNsmall~ Rep, data=qualdat, xlab="Rep", ylab="secondry flower number", main="boxplot of secondry flower number in each Rep", col="pink")
boxplot(TFN.~ Rep, data=qualdat, xlab="Rep", ylab="flower number", main="boxplot of flower number in each Rep", col="pink")
boxplot(Lg~ Rep, data=qualdat, xlab="Rep", ylab="Lodging", main="boxplot of Lodging in each Rep", col="pink")
boxplot(FD~ Rep, data=qualdat, xlab="Rep", ylab="Frost Damage", main="boxplot of Frost Damage in each Rep", col="pink")
boxplot(CCirc_cm~ Rep, data=qualdat, xlab="Rep", ylab="Compressed PLant Circumference", main="boxplot of Compressed PLant Circumference in each Rep", col="pink")
boxplot(SRD~ Rep, data=qualdat, xlab="Rep", ylab=" Re-growth Date", main="boxplot of Re-growth Date in each Rep", col="pink")
boxplot(ADD~ Rep, data=qualdat, xlab="Rep", ylab=" Autumn Dormancy Date", main="boxplot of Autumn Dormancy Date in each Rep", col="pink")
boxplot(GS~ Rep, data=qualdat, xlab="Rep", ylab=" Growth Stage", main="boxplot of Growth Stage in each Rep", col="pink")

####checking normality ####
###qq plot

par(mfrow=c(2,3))
qqnorm(qualdat$Yld_kg)
qqnorm(qualdat$SDW_kg)
qqnorm(qualdat$Cml_cm)
qqnorm(qualdat$CmD_BI_mm)
qqnorm(qualdat$CmN.)
qqnorm(qualdat$CmDW_g)
par(mfrow=c(2,3))
qqnorm(qualdat$Bcirc_cm)
qqnorm(qualdat$FNMain)
qqnorm(qualdat$FNsmall)
qqnorm(qualdat$TFN.)
qqnorm(qualdat$Lg)
qqnorm(qualdat$FD)
par(mfrow=c(2,3))
qqnorm(qualdat$CCirc_cm)
qqnorm(qualdat$SRD)
qqnorm(qualdat$ADD)
qqnorm(qualdat$GS)

###histogram
par(mfrow=c(2,3))
hist(qualdat$Yld_kg)
hist(qualdat$SDW_kg)
hist(qualdat$Cml_cm)
hist(qualdat$CmD_BI_mm)
hist(qualdat$CmN.)
hist(qualdat$CmDW_g)
par(mfrow=c(2,3))
hist(qualdat$Bcirc_cm)
hist(qualdat$FNMain)
hist(qualdat$FNsmall)
hist(qualdat$TFN.)
hist(qualdat$Lg)
hist(qualdat$FD)
par(mfrow=c(2,4))
hist(qualdat$CCirc_cm)
hist(qualdat$SRD)
hist(qualdat$ADD)
hist(qualdat$GS)

####Normality
#### Boxcox function R codes from online
source("~/Documents/R-corde for miscanthus project/bcplot function.txt")
par(mfrow=c(2,2))
bcplot(na.omit(qualdat$SDW_kg))####0.175
bcplot(na.omit(qualdat$CmDW_g))####0.15
bcplot(na.omit(qualdat$FD_1))#####-0.075
bcplot(na.omit(qualdat$FD_50.))####0.175
par(mfrow=c(2,2))
bcplot(na.omit(qualdat$HD_50.))####0.65
bcplot(na.omit(qualdat$HD_1))#####0.3
bcplot(na.omit(qualdat$HD_50.))####0.65
bcplot(na.omit(qualdat$TFN.)) #####0.15
par(mfrow=c(2,2))
bcplot(na.omit(qualdat$FNMain))####0.225
bcplot(na.omit(qualdat$FD)) #######0.55
bcplot(na.omit(qualdat$SRD))#######0.4
bcplot(na.omit(qualdat$ADD)) #######-1


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
qualdat$SDW_kg<- bc(qualdat$SDW_kg,lda=0.175)
qualdat$CmDW_g <- bc(qualdat$CmDW_g,lda=0.15)####0.15
qualdat$FD_1 <- bc(qualdat$FD_1,lda=-0.075)#####-0.075
qualdat$FD_50. <- bc(qualdat$FD_50.,lda=0.175)####0.175
qualdat$HD_50. <- bc(qualdat$HD_50.,lda=0.65)####0.65
qualdat$HD_1 <- bc(qualdat$HD_1,lda=0.3)#####0.3
qualdat$HD_50. <- bc(qualdat$HD_50.,lda=0.65)####0.65
qualdat$TFN. <- bc(qualdat$TFN.,lda=0.15) #####0.15
qualdat$FNMain <- bc(qualdat$FNMain,lda=0.225)####0.225
qualdat$FD <- bc(qualdat$FD,lda=0.55) #######0.55
qualdat$SRD <- bc(qualdat$SRD,lda=0.4)#######0.4
qualdat$ADD <- bc(qualdat$ADD,lda=-1) #######-1
str(qualdat)
####write the data out 
write.csv(qualdat, file = "~/Documents/whole traits/traits1718normalited.csv",row.names = T, na = ".")