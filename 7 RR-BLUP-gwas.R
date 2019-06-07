####rrBLUP method
####import the genotype data
myY <- read.csv("data/myYimputedSNP19.csv")
myGD <- read.csv("data/myGDimputedSNP19.csv",row.names=1)
myGM <- read.csv("data/myGMimputedSNP19.csv")
myQ<- read.csv("data/myQimputedSNP19.csv",row.names=1)
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
gwasResults2 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                     min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults2, file = "rrBLUP/rrBLUP_GWAS_results2.csv")

###this one for qq plots
out_start=4
out_end=41
for (i in out_start:out_end){
  mypath <- file.path("~/Documents/R-corde for miscanthus project/Miscanthus/rrBLUP",
                      paste("QQ plot of VAR_", names(gwasResults2[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults2[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=41
for (i in out_start:out_end){
  mypath <- file.path("~/Documents/R-corde for miscanthus project/Miscanthus/rrBLUP",
                      paste("Manhattan plots of _", names(gwasResults2[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults2$Name,
                       CHR = gwasResults2$Chromosome,
                       BP = gwasResults2$Position,
                       P = 10 ^ -gwasResults2[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}


###this one for OWA
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYimputedSNP19OWA3.csv")
str(myY)
myGD <- read.csv("data/myGDimputedSNP19OWA3.csv", row.names=1)
str(myGD)
myGM <- read.csv("data/myGMimputedSNP19OWA3.csv")
str(myGM)
myQ<- read.csv("data/myQimputedSNP19OWA3.csv",row.names=1)
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
gwasResultsowa <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                       min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResultsowa, file = "rrBLUP/rrBLUP_GWAS_resultsowa3.csv")

###this one for qq plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUP",
                      paste("QQ plot of VAR3YR_", names(gwasResultsowa[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResultsowa[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUP",
                      paste("Manhattan plots of 3YR_", names(gwasResultsowa[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResultsowa$Name,
                       CHR = gwasResultsowa$Chromosome,
                       BP = gwasResultsowa$Position,
                       P = 10 ^ -gwasResultsowa[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}

