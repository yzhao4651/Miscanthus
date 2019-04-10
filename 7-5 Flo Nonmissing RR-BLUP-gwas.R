####rrBLUP method
####import the genotype data with NON missing values SNPs 
myY <- read.csv("data/myYflo.106.4178.csv")
myGD <- read.csv("data/myGDflo.106.4178.csv",row.names=1)
myGM <- read.csv("data/myGMflo.106.4178.csv")
myQ<- read.csv("data/myQflo.106.4178.csv",row.names=1)
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
write.csv(gwasResults.flo106, file = "rrBLUP/rrBLUP_GWAS_results.flo106.csv")
str(gwasResults.flo106)

###this one for qq plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUP",
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
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUP",
                      paste("flo106 Manhattan plots of _", names(gwasResults.flo106[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.flo106$Name,
                       CHR = gwasResults.flo106$Chromosome,
                       BP = gwasResults.flo106$Position,
                       P = 10 ^ -gwasResults.flo106[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}



####import the genotype data with NON missing values SNPs 
myY <- read.csv("data/myYflo.115.3154.csv")
myGD <- read.csv("data/myGDflo.115.3154.csv",row.names=1)
myGM <- read.csv("data/myGMflo.115.3154.csv")
myQ<- read.csv("data/myQflo.115.3154.csv",row.names=1)
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
gwasResults.flo115 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                           min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.flo115, file = "rrBLUP/rrBLUP_GWAS_results.flo115.csv")

###this one for qq plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUP",
                      paste("flo115 QQ plot of VAR_", names(gwasResults.flo115[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.flo115[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=26
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUP",
                      paste("flo115 Manhattan plots of _", names(gwasResults.flo115[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.flo115$Name,
                       CHR = gwasResults.flo115$Chromosome,
                       BP = gwasResults.flo115$Position,
                       P = 10 ^ -gwasResults.flo115[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}
###Culm 106
###Culm 106
####import the genotype data with NON missing values SNPs 
myY <- read.csv("data/myYCulm.106.4220.csv")
myGD <- read.csv("data/myGDCulm.106.4220.csv",row.names=1)
myGM <- read.csv("data/myGMCulm.106.4220.csv")
myQ<- read.csv("data/myQCulm.106.4220.csv",row.names=1)
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
gwasResults.Clum106 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                           min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.Clum106, file = "rrBLUP/rrBLUP_GWAS_results.Clum106.csv")
str(gwasResults.Clum106)

###this one for qq plots
out_start=4
out_end=17
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUP",
                      paste("Clum106 QQ plot of VAR_", names(gwasResults.Clum106[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.Clum106[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=17
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUP",
                      paste("Clum106 Manhattan plots of _", names(gwasResults.Clum106[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.Clum106$Name,
                       CHR = gwasResults.Clum106$Chromosome,
                       BP = gwasResults.Clum106$Position,
                       P = 10 ^ -gwasResults.Clum106[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}


###Culm
###Culm 123
####import the genotype data with NON missing values SNPs 
myY <- read.csv("data/myYCulm.123.2670.csv")
myGD <- read.csv("data/myGDCulm.123.2670.csv",row.names=1)
myGM <- read.csv("data/myGMCulm.123.2670.csv")
myQ<- read.csv("data/myQCulm.123.2670.csv",row.names=1)
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
gwasResults.Clum123 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                            min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.Clum123, file = "rrBLUP/rrBLUP_GWAS_results.Clum123.csv")

###this one for qq plots
out_start=4
out_end=17
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUP",
                      paste("Clum123 QQ plot of VAR_", names(gwasResults.Clum123[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.Clum123[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=17
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUP",
                      paste("Clum123 Manhattan plots of _", names(gwasResults.Clum123[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.Clum123$Name,
                       CHR = gwasResults.Clum123$Chromosome,
                       BP = gwasResults.Clum123$Position,
                       P = 10 ^ -gwasResults.Clum123[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}
