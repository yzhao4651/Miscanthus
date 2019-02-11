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






