####rrBLUP method
####import the genotype data with NON missing values SNPs 
###flowering 119
myY <- read.csv("data/myYf.119.2537.csv")
myGD <- read.csv("data/myGDf.119.2537.csv", row.names=1)
myGM <- read.csv("data/myGMf.119.2537.csv")
myQ<- read.csv("data/myQf.119.2537.csv", row.names=1)
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
gwasResults.flo119 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                     min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.flo119, file = "rrBLUPF119/rrBLUP_GWAS_results.flo119.csv")
str(gwasResults.flo119)

###this one for qq plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPF119",
                      paste("flo119 QQ plot of VAR_", names(gwasResults.flo119[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.flo119[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPF119",
                      paste("flo119 Manhattan plots of _", names(gwasResults.flo119[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.flo119$Name,
                       CHR = gwasResults.flo119$Chromosome,
                       BP = gwasResults.flo119$Position,
                       P = 10 ^ -gwasResults.flo119[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}


###flowering 11
####import the genotype data with NON missing values SNPs 
myY <- read.csv("data/myYf.106.4185.csv")
myGD <- read.csv("data/myGDf.106.4185.csv",row.names=1)
myGM <- read.csv("data/myGMf.106.4185.csv")
myQ<- read.csv("data/myQf.106.4185.csv",row.names=1)
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
gwasResults.flo106 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                           min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.flo106, file = "rrBLUPF106/rrBLUP_GWAS_results.flo106.csv")

###this one for qq plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPF106",
                      paste("flo106 QQ plot of VAR_", names(gwasResults.flo106[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.flo106[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPF106",
                      paste("flo106 Manhattan plots of _", names(gwasResults.flo106[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.flo106$Name,
                       CHR = gwasResults.flo106$Chromosome,
                       BP = gwasResults.flo106$Position,
                       P = 10 ^ -gwasResults.flo106[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}
###f.116.3077
###f.116.3077
####import the genotype data with NON missing values SNPs 
myY <- read.csv("data/myYf.116.3077.csv")
myGD <- read.csv("data/myGDf.116.3077.csv",row.names=1)
myGM <- read.csv("data/myGMf.116.3077.csv")
myQ<- read.csv("data/myQf.116.3077.csv",row.names=1)
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
gwasResults.flo106 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                           min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.flo106, file = "rrBLUPF106/rrBLUP_GWAS_results.flo106.csv")
str(gwasResults.flo106)

###this one for qq plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPF106",
                      paste("Flo 106 QQ plot of VAR_", names(gwasResults.flo106[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.flo106[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPF106",
                      paste("Flo 106 Manhattan plots of _", names(gwasResults.flo106[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.flo106$Name,
                       CHR = gwasResults.flo106$Chromosome,
                       BP = gwasResults.flo106$Position,
                       P = 10 ^ -gwasResults.flo106[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}


###f.122.2272
###f.122.2272
####import the genotype data with NON missing values SNPs 
myY <- read.csv("data/myYf.122.2272.csv")
myGD <- read.csv("data/myGDf.122.2272.csv",row.names=1)
myGM <- read.csv("data/myGMf.122.2272.csv")
myQ<- read.csv("data/myQf.122.2272.csv",row.names=1)
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
gwasResults.flo122 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                            min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.flo122, file = "rrBLUPF122/rrBLUP_GWAS_results.flo122.csv")

###this one for qq plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPF122",
                      paste("flo122 QQ plot of VAR_", names(gwasResults.flo122[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.flo122[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPF122",
                      paste("flo122 Manhattan plots of _", names(gwasResults.flo122[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.flo122$Name,
                       CHR = gwasResults.flo122$Chromosome,
                       BP = gwasResults.flo122$Position,
                       P = 10 ^ -gwasResults.flo122[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}
