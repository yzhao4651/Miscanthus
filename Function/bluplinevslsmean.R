####function for blup line VS lsmean line
### x is the data set used for calculating the lsmean
### y is the data set with BLUP of all traits 
### pdf(paste("Linemean VS LineBlup", 1 ,".pdf",sep="")) 
### write before the function when using this function 
bluplinevslsmean <- function(x,y,na.rm=T){
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