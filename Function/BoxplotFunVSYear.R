####function PlotboxFunc loop to get the box for each variable 
####out_start: the first variable for boxplot, 
####out_end: the last variable for boxplot
#### writing the output pdf file like "pdf(paste("Boxplot of all traits", 1 ,".pdf",sep=""))  "
####in front of the function when citing this function. 
####In this way, the number 1 can be changed to the different number when using different variables 
#### in case the file will be overlap every time to use
BoxplotFunVSYear <- function(out_start,out_end, mydata, na.rm = TRUE) {
  par(mar = c(5,4,1,1))
  par(mfrow=c(3,3))
  for (i in out_start:out_end){
    boxplot(mydata[,i] ~ mydata$Year, col=(c("pink")), 
            xlab="", ylab="", main=paste("Boxplot of ", names(mydata[i]), sep = ""), cex.main=1)
    mtext(text = "Year",side = 1, line = 1.5, cex=0.7)
    mtext(text = paste("var_", names(mydata[i]), sep = ""), side = 2, line = 2,cex=0.7) #side 2 = left
  }
  dev.off()
}
