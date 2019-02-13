####import the data####
library(readr)
qualdat  <- read.csv("~/Documents/whole traits/Copy of trait1718-4.csv" , na.strings = c("",".","NA")) # on Yongli's computer
# mywd <- "." # on Lindsay's computer
qualdat <- read.csv(file.path(mywd, "Copy of trait1718-4.csv"), na.strings = c("",".","NA"))
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
###change several variables format
qualdat$GS <- as.numeric(as.character(qualdat$GS))
qualdat$Entry=as.factor(qualdat$Entry)
qualdat$Rep=as.factor(qualdat$Rep)
qualdat$Year=as.factor(qualdat$Year)
####check the data formate
str(qualdat)
#####data summary
#####data summary
####data summary Possible functions used in sapply include mean, sd, var, min, max, median, range, and quantile
install.packages("pastecs")
library(pastecs)
options(scipen=100)
options(digits=2)
stat.desc(qualdat)
summary <-stat.desc(qualdat)
write.csv(summary,file="~/Documents/whole traits/summary.csv",row.names = T,
             eol = "\n", na = "NA")
#write.csv(qualdat1,file="~/Documents/whole traits/data1718updated.csv",row.names = T, na = ".")

####scatterplot
####scatterplot
####Scatterplot Matrices for flower traits ( can hep to check some outlies)
pdf(paste("Scatterplot", 1 ,".pdf",sep="")) 
pairs(~HD_1+FD_1+HD_50.+FD_50.,data=qualdat, 
      main="Simple Scatterplot Matrix for flowering traits")
dev.off()
####scater plots for other traits 
pdf(paste("Scatterplot", 2 ,".pdf",sep="")) 
pairs(~Yld_kg+SDW_kg+Cml_cm+CmD_BI_mm+CmD_LI_mm+CmN.+CmDW_g+Bcirc_cm+FNMain+FNsmall+TFN.+FD+CCirc_cm+GS+Lg,data=qualdat, 
      main="Simple Scatterplot Matrix for flowering traits")
dev.off()

#####BoxPlot
#####BoxPlot
####using function PlotboxFunc loop to get the box for each variable 
####out_start: the first variable for boxplot, 
####out_end: the last variable for boxplot
plotboxFunc <- function(out_start,out_end, mydata, na.rm = TRUE) {
  pdf(paste("Boxplot", 1 ,".pdf",sep="")) 
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

####qq plot and histogram 
####qq plot and histogram 
qualdat1 <- qualdat[7:27]
histqq <- function(mydata,na.rm=TRUE){
 pdf(paste("Q-Q normal plot and Histogram", 1 ,".pdf",sep="")) 
 par(mar=rep(2,4))
 par(mfrow=c(4,4))
 for(colName in names(mydata)){ 
 hist(mydata[,colName],100,col="lightblue",xlab=colName,main=paste0("Histogram of ",colName),cex.main=1) 
 qqnorm(mydata[,colName],main=paste0("Q-Q normal plot of ",colName),cex.main=1) 
 qqline(mydata[,colName])
 } 
dev.off()
}
histqq(qualdat1)

####Normality
####Normality
#### Boxcox function R codes from online to get lambda
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
###????question below
###I tried to add this lines "return(lambda.hat)", but it does not work, I copied the output into CSV 
###and then import it. I does not take me a lot of time, but it is better to return it as a csv or excel. 

### get the lambdal for the variable need do data transformation 
##quesiton: trying to write a loop but it does not work, do you have any idea?
## The loop was not constructed correctly, since there were not curly brackets
## after it.  I have fixed it, and also indexed by trait name to make things
## easier.
####This one from Lindsay's contribution is GREAT, Thank you so much. 

lda<- read.csv("~/Documents/whole traits/lambda1.csv",row.name=1)
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

####write the data out 
write.csv(qualdat_BC, file = "~/Documents/whole traits/traits1718normalited1.csv",row.names = T, na = ".")