############# this is prepare the dataset for mrMLMM software
############# this is prepare the dataset for mrMLMM software
allsnp2sub <- read.csv("data/allsnp2sub.csv")
###step 4 import the myGM dataset for the GAPIT to select matched Genotype from step 3
##Suv
myGMmrMLMM <- read.csv("data/myGMS.106.5163.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDS.106.5163.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))
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
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMM.106.5163.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYS.106.5163.csv")
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
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMM.106.5163.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMS.106.5163.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/S106")

myGMmrMLMM <- read.csv("data/myGMS.106.5163.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDS.106.5163.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))
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
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMM.106.5163.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYS.106.5163.csv")
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
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMM.106.5163.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMS.106.5163.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/S106")
myGMmrMLMM <- read.csv("data/myGMS.106.5163.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDS.106.5163.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))
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
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMM.106.5163.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYS.106.5163.csv")
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
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMM.106.5163.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMS.106.5163.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/S106")
myGMmrMLMM <- read.csv("data/myGMS.106.5163.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDS.106.5163.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))
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
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMM.106.5163.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYS.106.5163.csv")
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
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMM.106.5163.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMS.106.5163.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/S106")
myGMmrMLMM <- read.csv("data/myGMS.106.5163.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDS.106.5163.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))
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
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMM.106.5163.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYS.106.5163.csv")
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
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMM.106.5163.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMS.106.5163.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/S106")




