####import the data
####import the data
ranefvalueall<- read.csv("~/Documents/whole traits/ranefvalueall.csv",na.strings = c("",".","NA"))
###check data format
str(ranefvalueall)
###rename of the column name 
colnames(ranefvalueall)[colnames(ranefvalueall)=="X"] <- "Entry"
###change 0 value to NA
ranefvalueall[ranefvalueall== 0] <- NA
###check data format
str(ranefvalueall)
####qq plot and histogram ####approximately normally distribution 
####qq plot and histogram 
ranefvalueall1<- data.frame(ranefvalueall,row.names=1)
histqq <- function(mydata, na.rm = TRUE){
 pdf(paste("Q-Q normal plot and Histogram", 2 ,".pdf",sep="")) 
 par(mar=rep(2,4))
 par(mfrow=c(4,4))
 for(colName in names(mydata)){ 
  hist(mydata[,colName],100,col="lightblue",xlab=colName,main=paste0("Histogram of ",colName),cex.main=1) 
  qqnorm(mydata[,colName],main=paste0("Q-Q normal plot of ",colName),cex.main=1) 
  qqline(mydata[,colName])
 } 
dev.off()
}
histqq(ranefvalueall1)


####plot BLUP and LineBlup
####trying to write function to plot 
####Compare BLUP to line averages on a scatterplot
qualdat <- read.csv("~/Documents/whole traits/data1718updated1.csv",na.strings = c("",".","NA"),row.names=1)
####check the format of dataset
str(qualdat)
indx <- sapply(qualdat[,c(7:10,16,19:23,25:27)], is.integer)
qualdat[,c(7:10,16,19:23,25:27)][indx] <- lapply(qualdat[,c(7:10,16,19:23,25:27)][indx], function(x) as.numeric(as.character(x)))
###change several variables format
qualdat$GS <- as.numeric(as.character(qualdat$GS))
qualdat$Entry=as.factor(qualdat$Entry)
qualdat$Rep=as.factor(qualdat$Rep)
qualdat$Year=as.factor(qualdat$Year)
####check the format of dataset
str(qualdat)
blupline <- function(x,y,na.rm=T){
  pdf(paste("Linemean and LineBlup", 1 ,".pdf",sep="")) 
  par(mar = c(5,4,1,1))
  par(mfrow=c(3,3))
  for(colName in names(x)){
 lmean = tapply(y[,colName], y$Entry, na.rm=T, mean, data=y)
 plot(x[,colName], lmean[row.names(x)], col="blue", main=paste0("Var_",colName),cex.main=0.9,
      xlab="", ylab="" ) 
 mtext(text = paste0("LineBlup_", colName, sep = ""),side = 1, line = 2, cex=0.7)
 mtext(text = paste0("line mean_", colName, sep = ""), side = 2, line = 2,cex=0.7)### subsetted lmean to match LINEBLUP
 }
dev.off()
}

blupline(ranefvalueall1,qualdat)

#### PCA analysis with ppca for flowering trait
#### PCA analysis with ppca for flowering trait
##install.packages("BiocManager")
###source("https://bioconductor.org/biocLite.R") ###intall.packages("pcaMethods")
###biocLite("pcaMethods")
###library(pcaMethods)
BiocManager::install("pcaMethods") ###intall.packages("pcaMethods")
library(pcaMethods)
#citation("pcaMethods")
ranefvaluef<- ranefvalueall[!(is.na(ranefvalueall[2]))|!(is.na(ranefvalueall[3])) | !(is.na(ranefvalueall[4]))|!(is.na(ranefvalueall[5])),]
pc <- pca(ranefvaluef[2:5], nPcs=1, method="ppca",center = TRUE)
fblupimputed <- data.frame(completeObs(pc))
flocompletepc <- cbind(ranefvaluef[1],fblupimputed)
###PC 
fprin_comp <- prcomp(flocompletepc[2:5])
###PC screen plot
pdf(paste("PC screen plot", 1 ,".pdf",sep="")) 
par(mfrow=c(1,2))
plot(fprin_comp)
plot(fprin_comp,type="line", main=paste0("PC screen plot"), cex.main=0.9)
abline(h=1,lty=3, col="red")
dev.off()
###get pc1 for each individual
fprin1 <- cbind(ranefvaluef[1], fprin_comp$x[,1])
names(fprin1)
colnames(fprin1)[colnames(fprin1)=="fprin_comp$x[, 1]"] <- "fprin"
alltraitsblup <- plyr::join_all(list(ranefvalueall,fprin1), by='Entry')
write.csv(alltraitsblup, file = "~/Documents/whole traits/alltraitsblup.csv",row.names = T)
#str(fprin1)
#write.csv(fprin1, file = "~/Documents/whole traits/fprin1.csv",row.names = T)





