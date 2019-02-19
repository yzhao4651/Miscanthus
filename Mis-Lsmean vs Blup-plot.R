####This codes for plotting the Lsmean anD BLUP of all traits
####import the data
ranefvalueall<- read.csv("data/ranefvalueall.csv",na.strings = c("",".","NA"))
###check data format
str(ranefvalueall)
###rename of the column name 
colnames(ranefvalueall)[colnames(ranefvalueall)=="X"] <- "Entry"
###check data format
str(ranefvalueall)
####firstly using Q-Q plot and histogram to check if the BLUP is normally distribution
####qq plot and histogram ####approximately normally distribution 
####qq plot and histogram 
ranefvalueall1<- data.frame(ranefvalueall,row.names=1)
source("histqq.R")
pdf(paste("Q-Q normal plot and Histogram", 2 ,".pdf",sep=""))
histqq(ranefvalueall1)

####plot BLUP VS LineBlup
###firstly load the data of qualdat
load("qualdat.RData")
####check the format of dataset
str(qualdat)
####trying to write function to plot 
####Compare BLUP to line averages on a scatterplot
blupline <- function(x,y,na.rm=T){
  pdf(paste("Linemean VSnLineBlup", 1 ,".pdf",sep="")) 
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






