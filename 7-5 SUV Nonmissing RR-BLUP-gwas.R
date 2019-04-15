####rrBLUP method
####import the genotype data with NON missing values SNPs 
###S.106.5163
myY <- read.csv("data/myYS.106.5163.csv")
str(myY)
myGD <- read.csv("data/myGDS.106.5163.csv", row.names=1)
myGM <- read.csv("data/myGMS.106.5163.csv")
myQ<- read.csv("data/myQS.106.5163.csv", row.names=1)
str(myQ)
###tansform the genotype data set into rrblup data set
myGDt <- t(myGD)-1
#myGDtdataframe <- data.frame(myGDt)
#head(myGDtdataframe)
### change the row name into the first column 
myGDtdataframe <- data.frame(myGDt)
Name <- rownames(myGDtdataframe)
rownames(myGDtdataframe) <- NULL
myGDtdataframe <- cbind(Name,myGDtdataframe)
###combine the genotype with the chromosome position (dataset: myGM)
genorrblup <- merge(myGM,myGDtdataframe,by="Name")
###change the ID of the name of phenotype to same to the ID of the genotype ( like change "_" to "." )
head(myY)
as.character(myY$Taxa)
myY$Taxa <- gsub("-",".", myY$Taxa)
## myY$Taxa <- make.names(mY$Taxa) # make.names(mY$Taxa) would also work
#head(myY)
#str(myY)
###question: how can incorprated the population structure? 
###question: how can incorprated the population structure? 
###run rrblup
###run rrblup
#install.packages("rrBLUP")
#install.packages("qqman")
### this one is without Kin matrix, 
library(rrBLUP)
library(qqman)
myGDk <- as.matrix(myGD - 1)
dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
#dimnames(myGDk)[[1]] <- make.names(dimnames(myGDk)[[1]])
###get kinship matrix
kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                 n.core=1,shrink=FALSE,return.imputed=FALSE)
####do GWAS analysis again with Kmatrix
gwasResults.S106 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                     min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.S106, file = "rrBLUPS106/rrBLUP_GWAS_results.S106.csv")
str(gwasResults.S106)

###this one for qq plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS106",
                      paste("S106 QQ plot of VAR_", names(gwasResults.S106[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.S106[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS106",
                      paste("S106 Manhattan plots of _", names(gwasResults.S106[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.S106$Name,
                       CHR = gwasResults.S106$Chromosome,
                       BP = gwasResults.S106$Position,
                       P = 10 ^ -gwasResults.S106[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}


###Culm127
myY <- read.csv("data/myYS.116.4685.csv")
str(myY)
myGD <- read.csv("data/myGDS.116.4685.csv", row.names=1)
myGM <- read.csv("data/myGMS.116.4685.csv")
myQ<- read.csv("data/myQS.116.4685.csv", row.names=1)
str(myQ)
###tansform the genotype data set into rrblup data set
myGDt <- t(myGD)-1
#myGDtdataframe <- data.frame(myGDt)
#head(myGDtdataframe)
### change the row name into the first column 
myGDtdataframe <- data.frame(myGDt)
Name <- rownames(myGDtdataframe)
rownames(myGDtdataframe) <- NULL
myGDtdataframe <- cbind(Name,myGDtdataframe)
###combine the genotype with the chromosome position (dataset: myGM)
genorrblup <- merge(myGM,myGDtdataframe,by="Name")
###change the ID of the name of phenotype to same to the ID of the genotype ( like change "_" to "." )
head(myY)
as.character(myY$Taxa)
myY$Taxa <- gsub("-",".", myY$Taxa)
## myY$Taxa <- make.names(mY$Taxa) # make.names(mY$Taxa) would also work
#head(myY)
#str(myY)
###question: how can incorprated the population structure? 
###question: how can incorprated the population structure? 
###run rrblup
###run rrblup
#install.packages("rrBLUP")
#install.packages("qqman")
### this one is without Kin matrix, 
library(rrBLUP)
library(qqman)
myGDk <- as.matrix(myGD - 1)
dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
#dimnames(myGDk)[[1]] <- make.names(dimnames(myGDk)[[1]])
###get kinship matrix
kmatrix <- A.mat(myGDk,min.MAF=0.05, max.missing=NULL,impute.method="mean",tol=0.02,
                 n.core=1,shrink=FALSE,return.imputed=FALSE)
####do GWAS analysis again with Kmatrix
gwasResults.S116 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                         min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.S116, file = "rrBLUPS116/rrBLUP_GWAS_results.S116.csv")
str(gwasResults.S116)

###this one for qq plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS116",
                      paste("S116 QQ plot of VAR_", names(gwasResults.S116[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.S116[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS116",
                      paste("S116 Manhattan plots of _", names(gwasResults.S116[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.S116$Name,
                       CHR = gwasResults.S116$Chromosome,
                       BP = gwasResults.S116$Position,
                       P = 10 ^ -gwasResults.S116[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}



###S.125.4182
myY <- read.csv("data/myYS.125.4182.csv")
str(myY)
myGD <- read.csv("data/myGDS.125.4182.csv", row.names=1)
myGM <- read.csv("data/myGMS.125.4182.csv")
myQ<- read.csv("data/myQS.125.4182.csv", row.names=1)
str(myQ)
###tansform the genotype data set into rrblup data set
myGDt <- t(myGD)-1
#myGDtdataframe <- data.frame(myGDt)
#head(myGDtdataframe)
### change the row name into the first column 
myGDtdataframe <- data.frame(myGDt)
Name <- rownames(myGDtdataframe)
rownames(myGDtdataframe) <- NULL
myGDtdataframe <- cbind(Name,myGDtdataframe)
###combine the genotype with the chromosome position (dataset: myGM)
genorrblup <- merge(myGM,myGDtdataframe,by="Name")
###change the ID of the name of phenotype to same to the ID of the genotype ( like change "_" to "." )
head(myY)
as.character(myY$Taxa)
myY$Taxa <- gsub("-",".", myY$Taxa)
## myY$Taxa <- make.names(mY$Taxa) # make.names(mY$Taxa) would also work
#head(myY)
#str(myY)
###question: how can incorprated the population structure? 
###question: how can incorprated the population structure? 
###run rrblup
###run rrblup
#install.packages("rrBLUP")
#install.packages("qqman")
### this one is without Kin matrix, 
library(rrBLUP)
library(qqman)
myGDk <- as.matrix(myGD - 1)
dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
#dimnames(myGDk)[[1]] <- make.names(dimnames(myGDk)[[1]])
###get kinship matrix
kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                 n.core=1,shrink=FALSE,return.imputed=FALSE)
####do GWAS analysis again with Kmatrix
gwasResults.S125 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                         min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.S125, file = "rrBLUPS125/rrBLUP_GWAS_results.S125.csv")
str(gwasResults.S125)

###this one for qq plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS125",
                      paste("S125 QQ plot of VAR_", names(gwasResults.S125[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.S125[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS125",
                      paste("S125 Manhattan plots of _", names(gwasResults.S125[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.S125$Name,
                       CHR = gwasResults.S125$Chromosome,
                       BP = gwasResults.S125$Position,
                       P = 10 ^ -gwasResults.S125[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}



###S.135.3628
myY <- read.csv("data/myYS.135.3628.csv")
str(myY)
myGD <- read.csv("data/myGDS.135.3628.csv", row.names=1)
myGM <- read.csv("data/myGMS.135.3628.csv")
myQ<- read.csv("data/myQS.135.3628.csv", row.names=1)
str(myQ)
###tansform the genotype data set into rrblup data set
myGDt <- t(myGD)-1
#myGDtdataframe <- data.frame(myGDt)
#head(myGDtdataframe)
### change the row name into the first column 
myGDtdataframe <- data.frame(myGDt)
Name <- rownames(myGDtdataframe)
rownames(myGDtdataframe) <- NULL
myGDtdataframe <- cbind(Name,myGDtdataframe)
###combine the genotype with the chromosome position (dataset: myGM)
genorrblup <- merge(myGM,myGDtdataframe,by="Name")
###change the ID of the name of phenotype to same to the ID of the genotype ( like change "_" to "." )
head(myY)
as.character(myY$Taxa)
myY$Taxa <- gsub("-",".", myY$Taxa)
## myY$Taxa <- make.names(mY$Taxa) # make.names(mY$Taxa) would also work
#head(myY)
#str(myY)
###question: how can incorprated the population structure? 
###question: how can incorprated the population structure? 
###run rrblup
###run rrblup
#install.packages("rrBLUP")
#install.packages("qqman")
### this one is without Kin matrix, 
library(rrBLUP)
library(qqman)
myGDk <- as.matrix(myGD - 1)
dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
#dimnames(myGDk)[[1]] <- make.names(dimnames(myGDk)[[1]])
###get kinship matrix
kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                 n.core=1,shrink=FALSE,return.imputed=FALSE)
####do GWAS analysis again with Kmatrix
gwasResults.S135 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                         min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.S135, file = "rrBLUPS135/rrBLUP_GWAS_results.S135.csv")
str(gwasResults.S135)

###this one for qq plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS135",
                      paste("S135 QQ plot of VAR_", names(gwasResults.S135[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.S135[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS135",
                      paste("S135 Manhattan plots of _", names(gwasResults.S135[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.S135$Name,
                       CHR = gwasResults.S135$Chromosome,
                       BP = gwasResults.S135$Position,
                       P = 10 ^ -gwasResults.S135[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}

###S.145.2787
myY <- read.csv("data/myYS.145.2787.csv")
str(myY)
myGD <- read.csv("data/myGDS.145.2787.csv", row.names=1)
myGM <- read.csv("data/myGMS.145.2787.csv")
myQ<- read.csv("data/myQS.145.2787.csv", row.names=1)
str(myQ)
###tansform the genotype data set into rrblup data set
myGDt <- t(myGD)-1
#myGDtdataframe <- data.frame(myGDt)
#head(myGDtdataframe)
### change the row name into the first column 
myGDtdataframe <- data.frame(myGDt)
Name <- rownames(myGDtdataframe)
rownames(myGDtdataframe) <- NULL
myGDtdataframe <- cbind(Name,myGDtdataframe)
###combine the genotype with the chromosome position (dataset: myGM)
genorrblup <- merge(myGM,myGDtdataframe,by="Name")
###change the ID of the name of phenotype to same to the ID of the genotype ( like change "_" to "." )
head(myY)
as.character(myY$Taxa)
myY$Taxa <- gsub("-",".", myY$Taxa)
## myY$Taxa <- make.names(mY$Taxa) # make.names(mY$Taxa) would also work
#head(myY)
#str(myY)
###question: how can incorprated the population structure? 
###question: how can incorprated the population structure? 
###run rrblup
###run rrblup
#install.packages("rrBLUP")
#install.packages("qqman")
### this one is without Kin matrix, 
library(rrBLUP)
library(qqman)
myGDk <- as.matrix(myGD - 1)
dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
#dimnames(myGDk)[[1]] <- make.names(dimnames(myGDk)[[1]])
###get kinship matrix
kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                 n.core=1,shrink=FALSE,return.imputed=FALSE)
####do GWAS analysis again with Kmatrix
gwasResults.S145 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                         min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.S145, file = "rrBLUPS145/rrBLUP_GWAS_results.S145.csv")
str(gwasResults.S145)

###this one for qq plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS145",
                      paste("S145 QQ plot of VAR_", names(gwasResults.S145[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.S145[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS145",
                      paste("S145 Manhattan plots of _", names(gwasResults.S145[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.S145$Name,
                       CHR = gwasResults.S145$Chromosome,
                       BP = gwasResults.S145$Position,
                       P = 10 ^ -gwasResults.S145[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}


