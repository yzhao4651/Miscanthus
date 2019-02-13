
#### PCA analysis with ppca for flowering trait
#### PCA analysis with ppca for flowering trait
##install.packages("BiocManager")
###source("https://bioconductor.org/biocLite.R") ###intall.packages("pcaMethods")
###biocLite("pcaMethods")
###library(pcaMethods)
BiocManager::install("pcaMethods") ###intall.packages("pcaMethods")
library(pcaMethods)
####import the data
#ranefvalueall<- read.csv("~/Documents/whole traits/ranefvalueall.csv",na.strings = c("",".","NA"))
ranefvalueall<- read.csv("ranefvalueall.csv",na.strings = c("",".","NA"))
###check data format
str(ranefvalueall)
###rename of the column name 
colnames(ranefvalueall)[colnames(ranefvalueall)=="X"] <- "Entry"

###check data format
str(ranefvalueall)
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
#write.csv(alltraitsblup, file = "~/Documents/whole traits/alltraitsblup.csv",row.names = T)
write.csv(alltraitsblup, file = "alltraitsblup.csv",row.names = T)
#str(fprin1)
#write.csv(fprin1, file = "~/Documents/whole traits/fprin1.csv",row.names = T)





