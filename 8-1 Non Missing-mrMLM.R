############# this is prepare the dataset for mrMLMM software
############# this is prepare the dataset for mrMLMM software
allsnp2sub <- read.csv("data/allsnp2sub.csv")
##import the SNPs files
##Suv
##S.106.5163
myGMmrMLMM <- read.csv("data/myGMS.106.5163.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDS.106.5163.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMS.106.5163.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYS.106.5163.csv")
str(myYmrMlMM)
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMS.106.5163.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMS.106.5163.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMS.106.5163.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/S106")

####S.116.4685
####S.116.4685
myGMmrMLMM <- read.csv("data/myGMS.116.4685.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDS.116.4685.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMS.116.4685.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYS.116.4685.csv")
str(myYmrMlMM)
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMS.116.4685.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMS.116.4685.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMS.116.4685.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/S116")

####S.125.4182
####S.125.4182
myGMmrMLMM <- read.csv("data/myGMS.125.4182.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDS.125.4182.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMS.125.4182.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYS.125.4182.csv")
str(myYmrMlMM)
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMS.125.4182.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMS.125.4182.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMS.125.4182.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/S125")


###S.135.3628
###S.135.3628
myGMmrMLMM <- read.csv("data/myGMS.135.3628.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDS.135.3628.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMS.135.3628.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYS.135.3628.csv")
str(myYmrMlMM)
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMS.135.3628.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMS.135.3628.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMS.135.3628.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/S135")

####S.145.2787
####S.145.2787
myGMmrMLMM <- read.csv("data/myGMS.145.2787.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDS.145.2787.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMS.145.2787.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYS.145.2787.csv")
str(myYmrMlMM)
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMS.145.2787.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMS.145.2787.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMS.145.2787.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/S145")

###Trait for Flowering traits
###Trait for Flowering traits
allsnp2sub <- read.csv("data/allsnp2sub.csv")
##import the SNPs files
##import the SNPs files
myGMmrMLMM <- read.csv("data/myGMf.119.2537.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDf.119.2537.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMf.119.2537.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYf.119.2537.csv")
str(myYmrMlMM)

myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMf.119.2537.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMf.119.2537.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMf.119.2537.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:24,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/F119")

####f.106.4185
####f.106.4185
myGMmrMLMM <- read.csv("data/myGMf.106.4185.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDf.106.4185.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMf.106.4185.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYf.106.4185.csv")
str(myYmrMlMM)
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMf.106.4185.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMf.106.4185.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMf.106.4185.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:24,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/F106")

####f.116.3077
####f.116.3077
myGMmrMLMM <- read.csv("data/myGMf.116.3077.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDf.116.3077.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMf.116.3077.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYf.116.3077.csv")
str(myYmrMlMM)
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMf.116.3077.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMf.116.3077.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMf.116.3077.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:24,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/F116")


###f.122.2272
###f.122.2272
myGMmrMLMM <- read.csv("data/myGMf.122.2272.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDf.122.2272.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMf.122.2272.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYf.122.2272.csv")
str(myYmrMlMM)
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMf.122.2272.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMf.122.2272.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMf.122.2272.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:24,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/F122")


###Clum data sets
####C.130.2001
####C.130.2001
myGMmrMLMM <- read.csv("data/myGMC.130.2001.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDC.130.2001.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMC.130.2001.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYC.130.2001.csv")
str(myYmrMlMM)
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMC.130.2001.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMC.130.2001.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMC.130.2001.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:15,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/C130")

####C.127.2216
####C.127.2216
myGMmrMLMM <- read.csv("data/myGMC.127.2216.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDC.127.2216.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMC.127.2216.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYC.127.2216.csv")
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMC.127.2216.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMC.127.2216.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMC.127.2216.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:15,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/C127")

####C.125.2562
####C.125.2562
myGMmrMLMM <- read.csv("data/myGMC.125.2562.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDC.125.2562.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMC.125.2562.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYC.125.2562.csv")
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMC.125.2562.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMC.125.2562.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMC.125.2562.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:15,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/C125")

####C.119.2988
####C.119.2988
myGMmrMLMM <- read.csv("data/myGMC.119.2988.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDC.119.2988.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMC.119.2988.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYC.119.2988.csv")
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMC.119.2988.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMC.119.2988.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMC.119.2988.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:15,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/C119")


####C.106.4202
####C.106.4202
myGMmrMLMM <- read.csv("data/myGMC.106.4202.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDC.106.4202.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMC.106.4202.csv", row.names = FALSE, na = "NA")
str(subgenomrMLMM)

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYC.106.4202.csv")
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMMC.106.4202.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMC.106.4202.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMC.106.4202.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:15,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/C106")