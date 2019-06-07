
###for Culm data with 106.4834.csv
###for Culm data with 106.4834.csv
####import all of the genotype data
myY <- read.csv("data/myYO.106.4834.csv")
str(myY)

myY$Taxa <- make.names(myY$Taxa)
colnames(myY)[which(names(myY) == "Taxa")] <- "<phenotype>"
source("Function/Changecolname.R")
myY <- colRename(myY)
str(myY)
write.csv(myY, file = "mrMLMM2/mrmyYO.106.4834.csv", row.names = FALSE, na = "NA")
myGD <- read.csv("data/myGDO.106.4834.csv")
str(myGD)
myGD <- data.frame(myGD,row.names = 1)
myGDtran <- data.frame(t(myGD))
str(myGDtran)
source("Function/FirstColumn.R")
myGDtran <- FirstColumn(myGDtran)
str(myGDtran)
colnames(myGDtran)[which(names(myGDtran) == "Taxa")] <- "rs"
myGM <- read.csv("data/myGMO.106.4834.csv")
str(myGM)
colnames(myGM)[which(names(myGM) == "Name")] <- "rs"
str(myGM)
myQ<- read.csv("data/myQO.106.4834.csv")
allsnp1 <- read.csv("data/allsnp1.csv")
str(allsnp1)
mrMLMMgeno <-allsnp1[match(colnames(myGD), allsnp1$rs, nomatch=0),]
mrMLMMgeno <- mrMLMMgeno[,c(2,5)]
str(mrMLMMgeno)
mrMLMMgeno2 <- plyr::join_all(list(data.frame(myGM),mrMLMMgeno,myGDtran),by="rs")
str(mrMLMMgeno2)
###change the name in order to fit the software requirment
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "rs")] <- "rs#"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Chromosome")] <- "chrom"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Position")] <- "pos"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "genotype.for.code.1")] <- "genotype for code 1"
str(mrMLMMgeno2)
write.csv(mrMLMMgeno2, file = "mrMLMM2/genomrmyYO.106.4834.csv", row.names = FALSE, na = "NA")

install.packages("mrMLM")
install.packages("tibble")
install.packages("zoo")
install.packages("lpSolve")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\genomrmyYO.106.4834.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mrmyYO.106.4834.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/OWA106")

###for Culm data with 123.2670
###for Culm data with 123.2670
####import all of the genotype data
myY <- read.csv("data/myYO.116.4232.csv")
myY <- myY[1:2]
myY$Taxa <- make.names(myY$Taxa)
colnames(myY)[which(names(myY) == "Taxa")] <- "<phenotype>"
source("Function/Changecolname.R")
myY <- colRename(myY)
str(myY)
write.csv(myY, file = "mrMLMM2/mrmyYO.116.4232.csv", row.names = FALSE, na = "NA")
myGD <- read.csv("data/myGDO.116.4232.csv")
str(myGD)
myGD <- data.frame(myGD,row.names = 1)
myGDtran <- data.frame(t(myGD))
str(myGDtran)
source("Function/FirstColumn.R")
myGDtran <- FirstColumn(myGDtran)
str(myGDtran)
colnames(myGDtran)[which(names(myGDtran) == "Taxa")] <- "rs"
myGM <- read.csv("data/myGMO.116.4232.csv")
str(myGM)
colnames(myGM)[which(names(myGM) == "Name")] <- "rs"
str(myGM)
myQ<- read.csv("data/myQO.116.4232.csv")
allsnp1 <- read.csv("data/allsnp1.csv")
str(allsnp1)
mrMLMMgeno <-allsnp1[match(colnames(myGD), allsnp1$rs, nomatch=0),]
mrMLMMgeno <- mrMLMMgeno[,c(2,5)]
str(mrMLMMgeno)
mrMLMMgeno2 <- plyr::join_all(list(data.frame(myGM),mrMLMMgeno,myGDtran),by="rs")
str(mrMLMMgeno2)
###change the name in order to fit the software requirment
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "rs")] <- "rs#"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Chromosome")] <- "chrom"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Position")] <- "pos"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "genotype.for.code.1")] <- "genotype for code 1"
str(mrMLMMgeno2)
write.csv(mrMLMMgeno2, file = "mrMLMM2/genomrmyYO.116.4232.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\genomrmyYO.116.4232.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mrmyYO.116.4232.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/OWA116")

###for Culm data with 136.3387
###for Culm data with 136.3387
####import all of the genotype data
myY <- read.csv("data/myYO.126.3659.csv")
myY <- myY[1:2]
myY$Taxa <- make.names(myY$Taxa)
colnames(myY)[which(names(myY) == "Taxa")] <- "<phenotype>"
source("Function/Changecolname.R")
myY <- colRename(myY)
str(myY)
write.csv(myY, file = "mrMLMM2/mrmyYO.126.3659.csv", row.names = FALSE, na = "NA")
myGD <- read.csv("data/myGDO.126.3659.csv")
str(myGD)
myGD <- data.frame(myGD,row.names = 1)
myGDtran <- data.frame(t(myGD))
str(myGDtran)
source("Function/FirstColumn.R")
myGDtran <- FirstColumn(myGDtran)
str(myGDtran)
colnames(myGDtran)[which(names(myGDtran) == "Taxa")] <- "rs"
myGM <- read.csv("data/myGMO.135.3628.csv")
str(myGM)
colnames(myGM)[which(names(myGM) == "Name")] <- "rs"
str(myGM)
myQ<- read.csv("data/myQO.126.3659.csv")
allsnp1 <- read.csv("data/allsnp1.csv")
str(allsnp1)
mrMLMMgeno <-allsnp1[match(colnames(myGD), allsnp1$rs, nomatch=0),]
mrMLMMgeno <- mrMLMMgeno[,c(2,5)]
str(mrMLMMgeno)
mrMLMMgeno2 <- plyr::join_all(list(data.frame(myGM),mrMLMMgeno,myGDtran),by="rs")
str(mrMLMMgeno2)
###change the name in order to fit the software requirment
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "rs")] <- "rs#"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Chromosome")] <- "chrom"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Position")] <- "pos"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "genotype.for.code.1")] <- "genotype for code 1"
str(mrMLMMgeno2)
write.csv(mrMLMMgeno2, file = "mrMLMM2/genomrmyYO.126.3659.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\genomrmyYO.126.3659.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mrmyYO.126.3659.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/OWA126")

###for Culm data with 146.2728
###for Culm data with 146.2728
####import all of the genotype data
myY <- read.csv("data/myYO.135.2855.csv")
myY <- myY[1:2]
myY$Taxa <- make.names(myY$Taxa)
colnames(myY)[which(names(myY) == "Taxa")] <- "<phenotype>"
source("Function/Changecolname.R")
myY <- colRename(myY)
str(myY)
write.csv(myY, file = "mrMLMM2/mrmyYO.135.2855.csv", row.names = FALSE, na = "NA")
myGD <- read.csv("data/myGDO.135.2855.csv")
str(myGD)
myGD <- data.frame(myGD,row.names = 1)
myGDtran <- data.frame(t(myGD))
str(myGDtran)
source("Function/FirstColumn.R")
myGDtran <- FirstColumn(myGDtran)
str(myGDtran)
colnames(myGDtran)[which(names(myGDtran) == "Taxa")] <- "rs"
myGM <- read.csv("data/myGMO.135.2855.csv")
str(myGM)
colnames(myGM)[which(names(myGM) == "Name")] <- "rs"
str(myGM)
myQ<- read.csv("data/myQO.135.2855.csv")
allsnp1 <- read.csv("data/allsnp1.csv")
str(allsnp1)
mrMLMMgeno <-allsnp1[match(colnames(myGD), allsnp1$rs, nomatch=0),]
mrMLMMgeno <- mrMLMMgeno[,c(2,5)]
str(mrMLMMgeno)
mrMLMMgeno2 <- plyr::join_all(list(data.frame(myGM),mrMLMMgeno,myGDtran),by="rs")
str(mrMLMMgeno2)
###change the name in order to fit the software requirment
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "rs")] <- "rs#"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Chromosome")] <- "chrom"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "Position")] <- "pos"
colnames(mrMLMMgeno2)[which(names(mrMLMMgeno2) == "genotype.for.code.1")] <- "genotype for code 1"
str(mrMLMMgeno2)
write.csv(mrMLMMgeno2, file = "mrMLMM2/genomrmyYO.135.2855.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\genomrmyYO.135.2855.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mrmyYO.135.2855.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/OWA135")

