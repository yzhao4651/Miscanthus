####import the data####
library(readr)
qualdat  <- read.csv("~/Documents/whole traits/Copy of trait1718-8.csv" , na.strings = c("",".","NA"))
qualdat  <- read.csv("data/Copy of trait1718-8.csv" , na.strings = c("",".","NA")) # on Yongli's computer
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
###select the data need for analysis 
qualdat <- qualdat[,c(3,5:6,4,7:27)]
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
###save this data set 
save(qualdat,file="qualdat.RData")
####check the data formate
str(qualdat)
#####data summary
#####data summary
####data summary Possible functions used in sapply include mean, sd, var, min, max, median, range, and quantile
install.packages("pastecs")
library(pastecs)
options(scipen=100)
options(digits=4)
stat.desc(qualdat)
summary <-stat.desc(qualdat)
write.csv(summary,file="data/summary.csv",row.names = T,
             eol = "\n", na = "NA")

####scatterplot
Scatterplot<- function(mydata){
  panel.cor <- function(x, y, digits=4, prefix="", cex.cor) 
  {
    usr <- par("usr"); on.exit(par(usr)) 
    par(usr = c(0, 1, 0, 1)) 
    r <- abs(cor(x, y, use="pairwise")) 
    txt <- format(c(r, 0.123456789), digits=digits)[1] 
    txt <- paste(prefix, txt, sep="") 
    if(missing(cex.cor)) cex <- 0.9/strwidth(txt) 
    
    test <- cor.test(x, y, use="pairwise")
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " ")) 
    text(0.5, 0.5, txt, cex = cex * r) 
    text(.5, .8, Signif, cex = cex, col=2)
  }
  pairs(mydata, lower.panel=panel.smooth, upper.panel=panel.cor)
}
source("Function/Scatterplot.R")
####flowering traits 
pdf(paste("Scatterplot of flowering traits", 1 ,".pdf",sep="")) 
Scatterplot(qualdat[,5:8])  
dev.off() 
####other traits with 2017 years
pdf(paste("Scatterplot of 2017 traits", 1 ,".pdf",sep="")) 
Scatterplot(qualdat[,c(9:20)])  
dev.off() 
####other traits with 2018 years
pdf(paste("Scatterplot of 2018 traits", 1 ,".pdf",sep="")) 
Scatterplot(qualdat[,c(9:16,20:25)])  
dev.off() 
####other traits with 2017, 2018 years
pdf(paste("Scatterplot of 1718 traits", 1 ,".pdf",sep="")) 
Scatterplot(qualdat[,c(9:18,20)])  
dev.off() 
#####BoxPlot
#####BoxPlot

source("Function/plotboxFunc.R")
pdf(paste("Boxplot of all traits", 1 ,".pdf",sep=""))
plotboxFunc(4,25,qualdat)
dev.off()
####qq plot and histogram 
####qq plot and histogram 
source("Function/histqq.R")
pdf(paste("Q-Q normal plot and Histogram", 1 ,".pdf",sep=""))
histqq(qualdat[4:25])
dev.off()
####Normality
####Normality
#### Boxcox function R codes from online to get lambda
source("Function/bcplot function.txt")
pdf(paste("lambda", 1 ,".pdf",sep=""))####question: only produce one image --> Did we solve this earlier?  seems to product all images now
par(mar=rep(2,4))
par(mfrow=c(4,4))
out_start=5
out_end=25
# set up empty data frame to contain lambda values
lda <- data.frame(row.names = colnames(qualdat)[out_start:out_end],
                  lambda = rep(NA_real_, out_end - out_start + 1))
for (i in out_start:out_end){ 
  trtname <- colnames(qualdat)[i]
  lda[trtname, 1] <- bcplot(na.omit(qualdat[i])) # put result of bcplot into data frame
} 
dev.off()
write.csv(lda,file="data/lda.csv")
###question: how i can get a csv file with column name with lambda?
### Unfortunately the bcplot function does not return any value, so it can't 
### be used programatically to send values to a vector that can be written to a
### file.  You could possible edit the function with a statement like
### `return(lambda.hat)`
###????question below
###I tried to add this lines "return(lambda.hat)", but it does not work, I copied the output into CSV 
###and then import it. I does not take me a lot of time, but it is better to return it as a csv or excel. 
### --> Lindsay's answer
### bcplot was returning lambda.hat, but you didn't update your script to save that value
### anywhere.  I have edited it.

### get the lambdal for the variable need do data transformation 
##quesiton: trying to write a loop but it does not work, do you have any idea?
## The loop was not constructed correctly, since there were not curly brackets
## after it.  I have fixed it, and also indexed by trait name to make things
## easier.
####This one from Lindsay's contribution is GREAT, Thank you so much. 

lda<- read.csv("data/lda.csv",row.name=1)
lda <- read.csv(file.path(mywd, "lambda1.csv"), row.name=1)

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
source("Function/bc1.R")
qualdat_BC <- bc1(qualdat, lda)

####write the data out 
write.csv(qualdat_BC, file = "data/traits1718normalited1.csv",row.names = T, na = ".")

