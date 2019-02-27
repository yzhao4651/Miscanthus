####function of Histogram and Q-Q plot for each variable with loop
#### writing the output pdf file "pdf(paste("Q-Q normal plot and Histogram", 1 ,".pdf",sep="")) "
####in front of the function when citing this function. 
####In this way, the number 1 can be changed to the different number when using different variables 
#### in case the file will be overlap every time to use
histqq <- function(mydata){
  par(mar=rep(2,4))
  par(mfrow=c(4,4))
  for(colName in names(mydata)){ 
    hist(mydata[,colName],100,col="lightblue",xlab=colName,main=paste0("Histogram of ",colName),cex.main=1) 
    qqnorm(mydata[,colName],main=paste0("Q-Q normal plot of ",colName),cex.main=1) 
    qqline(mydata[,colName])
  }
dev.off()
}
