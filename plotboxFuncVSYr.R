plotboxFuncVSYear <- function(out_start,out_end, mydata, na.rm = TRUE) {
  par(mar = c(5,4,1,1))
  par(mfrow=c(3,3))
  for (i in out_start:out_end){
    boxplot(mydata[,i] ~ mydata$Year, col=(c("pink")), 
            xlab="", ylab="", main=paste("Boxplot of ", names(mydata[i]), sep = ""), cex.main=1)
    mtext(text = "Rep",side = 1, line = 1.5, cex=0.7)
    mtext(text = paste("var_", names(mydata[i]), sep = ""), side = 2, line = 2,cex=0.7) #side 2 = left
  }
  dev.off()
}