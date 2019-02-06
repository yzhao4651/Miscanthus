####import the data####
library(readr)
mywd <- "~/Documents/whole traits" # on Yongli's computer
# mywd <- "." # on Lindsay's computer
qualdat <- read.csv(file.path(mywd, "trait1718.csv"), na.strings = c("",".","NA"))
#check the data formate
str(qualdat)
###change character of Growth Stage to number N=1,B=2,F=3,P=4
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
# Note from Lindsay -- converting from integer to float is unnecessary, but should not cause problems
indx <- sapply(qualdat[,c(16,19:23)], is.integer)
qualdat[,c(16,19:23)][indx] <- lapply(qualdat[,c(16,19:23)][indx], function(x) as.numeric(as.character(x)))

qualdat$GS <- as.numeric(as.character(qualdat$GS))
qualdat$Entry=as.factor(qualdat$Entry)
qualdat$Rep=as.factor(qualdat$Rep)
qualdat$Year=as.factor(qualdat$Year)
####check the data formate
str(qualdat)
#write.csv(qualdat1,file="~/Documents/whole traits/data1718updated.csv",row.names = T, na = ".")
####Scatterplot Matrices for flower traits ( can hep to check some outlies)
pairs(~HD_1+FD_1+HD_50.+FD_50.,data=qualdat, 
      main="Simple Scatterplot Matrix for flowering traits")
####scater plots for other traits 
pairs(~Yld_kg+SDW_kg+Cml_cm+CmD_BI_mm+CmD_LI_mm+CmN.+CmDW_g+Bcirc_cm+FNMain+FNsmall+TFN.+FD+CCirc_cm+GS+Lg,data=qualdat, 
      main="Simple Scatterplot Matrix for flowering traits")
####using function PlotboxFunc loop to get the box for each variable 
####out_start: the first variable for boxplot, 
####out_end: the last variable for boxplot
plotboxFunc <- function(out_start,out_end, mydata, na.rm = TRUE) {
  pdf(paste("Boxplot", 2 ,".pdf",sep="")) 
  par(mar = c(5,4,1,1))
  par(mfrow=c(3,3))
  for (i in out_start:out_end){
    boxplot(qualdat[,i] ~ qualdat$Rep, col=(c("pink")), 
            xlab="", ylab="", main=paste("Boxplot of ", names(qualdat[i]), sep = ""), cex.main=1)
    mtext(text = "Rep",side = 1, line = 1.5, cex=0.7)
    mtext(text = paste("var_", names(qualdat[i]), sep = ""), side = 2, line = 2,cex=0.7) #side 2 = left
  }
  dev.off()
}
plotboxFunc(7,27,qualdat)
####qq plot
qualdat1 <- qualdat[7:27]
pdf(paste("Q-Q normal plot and Histogram", 2 ,".pdf",sep="")) 
par(mar=rep(2,4))
par(mfrow=c(4,4))
for(colName in names(qualdat1)){ 
hist(qualdat1[,colName],100,col="lightblue",xlab=colName,main=paste0("Histogram of ",colName),cex.main=1) 
qqnorm(qualdat1[,colName],main=paste0("Q-Q normal plot of ",colName),cex.main=1) 
qqline(qualdat1[,colName])
} 
dev.off()

#####the function above works fine

#####
####Normality
#### Boxcox function R codes from online
####
source("bcplot function.txt")
pdf(paste("lambda", 2 ,".pdf",sep=""))####question: only produce one image --> Did we solve this earlier?  seems to product all images now
par(mar=rep(2,4))
par(mfrow=c(4,4))
out_start=7
out_end=27
for (i in out_start:out_end){ 
  bcplot(na.omit(qualdat[i]))
} 
dev.off()
###question: how i can get a csv file with column name with lambda?
### Unfortunately the bcplot function does not return any value, so it can't 
### be used programatically to send values to a vector that can be written to a
### file.  You could possible edit the function with a statement like
### `return(lambda.hat)`


### get the lambdal for the variable need do data transformation 
##quesiton: trying to write a loop but it does not work, do you have any idea?
## The loop was not constructed correctly, since there were not curly brackets
## after it.  I have fixed it, and also indexed by trait name to make things
## easier.
lda <- read.csv(file.path(mywd, "lambda1.csv"), row.name=1)
out_start=7
out_end=27
bc1 <- function(x, lda){ # function to transform data after lambda is determined
  for(trait in rownames(lda)){
    while(min(x[[trait]], na.rm = TRUE) <= 0){
      x[[trait]] <- x[[trait]] + 1 # positive numbers required
    }
    if(lda[trait,] == 0){
      x[[trait]] <- log(x[[trait]])
    } else {
      x[[trait]] <- (x[[trait]]^lda[trait,] - 1)/lda[trait,]
    }
  }
  
  return(x)
}

qualdat_BC <- bc1(qualdat, lda)

### becasue this function above does not work, I still write one by one to get the transformation data,
## do you have any idea for this one, It is really hard for me to do this one 
### Lindsay's code to get the transformed data
bc <- function(x, lda){ # function to transform data after lambda is determined
  while(min(x, na.rm = TRUE) <= 0){
    x <- x+ 1 # positive numbers required
  }
  if(lda == 0){
    x <- log(x)
  } else {
    x <- (x^lda - 1)/lda
  }
  return(x)
}
qualdat$HD_1 <- bc(qualdat$HD_1,lda=0.3)##### 0.3
qualdat$FD_1 <- bc(qualdat$FD_1,lda=-0.075)##### -0.075
qualdat$HD_50. <- bc(qualdat$HD_50.,lda=0.65)#### 0.65
qualdat$FD_50. <- bc(qualdat$FD_50.,lda=0.175)#### 0.175
qualdat$Yld_kg <- bc(qualdat$Yld_kg,lda=0)#### 0
qualdat$SDW_kg<- bc(qualdat$SDW_kg,lda=0.175)#### 0.175
qualdat$Cml_cm<- bc(qualdat$Cml_cm,lda=0.975)#### 0.975
qualdat$CmD_BI_mm<- bc(qualdat$CmD_BI_mm,lda=0.175)#### 0.175
qualdat$CmD_LI_mm<- bc(qualdat$CmD_LI_mm,lda=0.575)#### 0.575
qualdat$CmN.<- bc(qualdat$CmN.,lda=0.375)#### 0.375
qualdat$CmDW_g <- bc(qualdat$CmDW_g,lda=0.15)#### 0.15
qualdat$Bcirc_cm<- bc(qualdat$Bcirc_cm,lda=0.3)#### 0.3
qualdat$FNMain <- bc(qualdat$FNMain,lda=0.225)#### 0.225
qualdat$FNsmall<- bc(qualdat$FNsmall,lda=0.175)#### 0.175
qualdat$TFN. <- bc(qualdat$TFN.,lda=0.15) ##### 0.15
qualdat$Lg <- bc(qualdat$Lg,lda=-1) ####### 0.55
qualdat$FD <- bc(qualdat$FD,lda=0.55) ####### 0.55
qualdat$CCirc_cm <- bc(qualdat$CCirc_cm,lda=0)###### 0
qualdat$SRD <- bc(qualdat$SRD,lda=0.4)####### 0.4
qualdat$ADD <- bc(qualdat$ADD,lda=-1) ####### -1
qualdat$GS <- bc(qualdat$GS,lda=1) ####### 1
str(qualdat)

identical(qualdat, qualdat_BC) # this confirms results the same; can delete above code that does one at a time.
####write the data out 
write.csv(qualdat, file = "~/Documents/whole traits/traits1718normalited1.csv",row.names = T, na = ".")