
#### PCA analysis with ppca for flowering trait
#### PCA analysis with ppca for flowering trait
##install.packages("BiocManager")
###source("https://bioconductor.org/biocLite.R") ###intall.packages("pcaMethods")
###biocLite("pcaMethods")
###library(pcaMethods)
####import the data
#ranefvalueall<- read.csv("~/Documents/whole traits/ranefvalueall2.csv",na.strings = c("",".","NA"))
ranefvalueall<- read.csv("data/ranefvalueall2.csv",na.strings = c("",".","NA"))
###check data format
str(ranefvalueall)
###rename of the column name 
colnames(ranefvalueall)[colnames(ranefvalueall)=="X"] <- "Entry"
###check data format
str(ranefvalueall)

###PC1 for flowering labels as days
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("pcaMethods") ###intall.packages("pcaMethods")
#citation("pcaMethods")
ranefvaluef<- ranefvalueall[!(is.na(ranefvalueall[2]))|!(is.na(ranefvalueall[3])) | !(is.na(ranefvalueall[4]))|!(is.na(ranefvalueall[5])),]
library(pcaMethods)
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
fprind <- cbind(ranefvaluef[1], fprin_comp$x[,1])
names(fprind)
colnames(fprind)[colnames(fprind)=="fprin_comp$x[, 1]"] <- "fprind"


###PC1 for flowering labels as week
###PC1 for flowering labels as week
ranefvaluef<- ranefvalueall[!(is.na(ranefvalueall[24]))|!(is.na(ranefvalueall[25])) | !(is.na(ranefvalueall[26]))|!(is.na(ranefvalueall[27])),]
library(pcaMethods)
pc <- pca(ranefvaluef[24:27], nPcs=1, method="ppca",center = TRUE)
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
fprinW <- cbind(ranefvaluef[1], fprin_comp$x[,1])
names(fprinW)
colnames(fprinW)[colnames(fprinW)=="fprin_comp$x[, 1]"] <- "fprinW"

###PC1 for flowering labels as two weeks 
###PC1 for flowering labels as two weeks 
#citation("pcaMethods")
ranefvaluef<- ranefvalueall[!(is.na(ranefvalueall[28]))|!(is.na(ranefvalueall[29])) | !(is.na(ranefvalueall[30]))|!(is.na(ranefvalueall[31])),]
library(pcaMethods)
pc <- pca(ranefvaluef[28:31], nPcs=1, method="ppca",center = TRUE)
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
fprinGW <- cbind(ranefvaluef[1], fprin_comp$x[,1])
names(fprinGW)
colnames(fprinGW)[colnames(fprinGW)=="fprin_comp$x[, 1]"] <- "fprinGW"

###PC1 for flowering labels as month 
###PC1 for flowering labels as month  
#citation("pcaMethods")
ranefvaluef<- ranefvalueall[!(is.na(ranefvalueall[32]))|!(is.na(ranefvalueall[33])) | !(is.na(ranefvalueall[34]))|!(is.na(ranefvalueall[35])),]
pc <- pca(ranefvaluef[32:35], nPcs=1, method="ppca",center = TRUE)
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
fprinM <- cbind(ranefvaluef[1], fprin_comp$x[,1])
names(fprinM)
colnames(fprinM)[colnames(fprinM)=="fprin_comp$x[, 1]"] <- "fprinM"
###combine all ofthe traits together
alltraitsblup <- plyr::join_all(list(ranefvalueall,fprind,fprinW,fprinGW,fprinM), by='Entry')
str(alltraitsblup)
#write.csv(alltraitsblup, file = "~/Documents/whole traits/alltraitsblup.csv",row.names = T)
write.csv(alltraitsblup, file = "data/alltraitsblup2.csv",row.names = T)
#str(fprin1)
#write.csv(fprin1, file = "~/Documents/whole traits/fprin1.csv",row.names = T)
str(ranefvalueall)
#citation("pcaMethods")

