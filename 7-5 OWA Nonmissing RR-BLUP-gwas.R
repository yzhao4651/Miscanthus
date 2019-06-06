####rrBLUP method
####import the genotype data with NON missing values SNPs 
###S.106.5163
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYO.106.4834.csv")
str(myY)
myGD <- read.csv("data/myGDO.106.4834.csv", row.names=1)
myGM <- read.csv("data/myGMO.106.4834.csv")
myQ<- read.csv("data/myQO.106.4834.csv", row.names=1)
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
write.csv(gwasResults.S106, file = "rrBLUPS106/rrBLUP_GWAS_results.O106.csv")
str(gwasResults.S106)

###this one for qq plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS106",
                      paste("O106 QQ plot of VAR_", names(gwasResults.S106[i]), ".jpeg", sep = ""))
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
                      paste("O106 Manhattan plots of _", names(gwasResults.S106[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.S106$Name,
                       CHR = gwasResults.S106$Chromosome,
                       BP = gwasResults.S106$Position,
                       P = 10 ^ -gwasResults.S106[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}


gwasResults.O106 <- read.csv("rrBLUPS106/rrBLUP_GWAS_results.O106.csv",row.names=1)
str(gwasResults.O106)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("rrBLUP.",colnames(data[i]))] <- 10^-data[i]
  }
  return(data)
}
gwasResults.O106 <- adj_P_function(gwasResults.O106, 4, 4)

write.csv(gwasResults.O106, file = "rrBLUPS106/rrBLUP_GWAS_results.O106.csv")


###myYO.116.4232.csv
myY <- read.csv("data/myYO.116.4232.csv")
str(myY)
myGD <- read.csv("data/myGDO.116.4232.csv", row.names=1)
myGM <- read.csv("data/myGMO.116.4232.csv")
myQ<- read.csv("data/myQO.116.4232.csv", row.names=1)
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
write.csv(gwasResults.S116, file = "rrBLUPS116/rrBLUP_GWAS_results.O116.csv")
str(gwasResults.S116)

###this one for qq plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS116",
                      paste("O116 QQ plot of VAR_", names(gwasResults.S116[i]), ".jpeg", sep = ""))
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
                      paste("O116 Manhattan plots of _", names(gwasResults.S116[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.S116$Name,
                       CHR = gwasResults.S116$Chromosome,
                       BP = gwasResults.S116$Position,
                       P = 10 ^ -gwasResults.S116[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}

gwasResults.O106 <- read.csv("rrBLUPS116/rrBLUP_GWAS_results.O116.csv",row.names=1)
str(gwasResults.O106)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("rrBLUP.",colnames(data[i]))] <- 10^-data[i]
  }
  return(data)
}
gwasResults.O106 <- adj_P_function(gwasResults.O106, 4, 4)

write.csv(gwasResults.O106, file = "rrBLUPS116/rrBLUP_GWAS_results.O116.csv")


###O.126.3659
myY <- read.csv("data/myYO.126.3659.csv")
str(myY)
myGD <- read.csv("data/myGDO.126.3659.csv", row.names=1)
myGM <- read.csv("data/myGMO.126.3659.csv")
myQ<- read.csv("data/myQO.126.3659.csv", row.names=1)
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
write.csv(gwasResults.S125, file = "rrBLUPS125/rrBLUP_GWAS_results.O125.csv")
str(gwasResults.S125)

###this one for qq plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS125",
                      paste("O125 QQ plot of VAR_", names(gwasResults.S125[i]), ".jpeg", sep = ""))
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
                      paste("O125 Manhattan plots of _", names(gwasResults.S125[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.S125$Name,
                       CHR = gwasResults.S125$Chromosome,
                       BP = gwasResults.S125$Position,
                       P = 10 ^ -gwasResults.S125[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}

gwasResults.O106 <- read.csv("rrBLUPS125/rrBLUP_GWAS_results.O125.csv",row.names=1)
str(gwasResults.O106)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("rrBLUP.",colnames(data[i]))] <- 10^-data[i]
  }
  return(data)
}
gwasResults.O106 <- adj_P_function(gwasResults.O106, 4, 4)

write.csv(gwasResults.O106, file = "rrBLUPS125/rrBLUP_GWAS_results.O126.csv")



###O.135.3628
myY <- read.csv("data/myYO.135.2855.csv")
str(myY)
myGD <- read.csv("data/myGDO.135.2855.csv", row.names=1)
myGM <- read.csv("data/myGMO.135.2855.csv")
myQ<- read.csv("data/myQO.135.2855.csv", row.names=1)
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
write.csv(gwasResults.S135, file = "rrBLUPS135/rrBLUP_GWAS_results.O135.csv")
str(gwasResults.S135)

###this one for qq plots
out_start=4
out_end=4
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPS135",
                      paste("O135 QQ plot of VAR_", names(gwasResults.S135[i]), ".jpeg", sep = ""))
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
                      paste("O135 Manhattan plots of _", names(gwasResults.S135[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.S135$Name,
                       CHR = gwasResults.S135$Chromosome,
                       BP = gwasResults.S135$Position,
                       P = 10 ^ -gwasResults.S135[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}


gwasResults.O106 <- read.csv("rrBLUPS135/rrBLUP_GWAS_results.O135.csv",row.names=1)
str(gwasResults.O106)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("rrBLUP.",colnames(data[i]))] <- 10^-data[i]
  }
  return(data)
}
gwasResults.O106 <- adj_P_function(gwasResults.O106, 4, 4)

write.csv(gwasResults.O106, file = "rrBLUPS135/rrBLUP_GWAS_results.O135.csv")





