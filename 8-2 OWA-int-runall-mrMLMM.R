###change some number in the genotype to int
subgenomrMLMM <- read.csv("data/subgenomrMLMM.csv")
subgenomrMLMM1 <- data.frame(subgenomrMLMM[,c(1,5:157)],row.names = 1)
###seperate two parts
subgenomrMLMMtran <- data.frame(t(subgenomrMLMM1))
subgenomrMLMMtran[subgenomrMLMMtran < 1 & subgenomrMLMMtran > 0.5] <- 1
subgenomrMLMMtran[subgenomrMLMMtran < 0.5 & subgenomrMLMMtran > 0] <- 0
subgenomrMLMMtran[subgenomrMLMMtran < -0.5 & subgenomrMLMMtran > -1] <- -1
subgenomrMLMMtran[subgenomrMLMMtran < 0 & subgenomrMLMMtran > -0.5] <- 0
###Suv traits
###SUV traits
mymrMLMMSuv <- read.csv("mrMLMM2/myYmrMlMMSuv.csv")
str(mymrMLMMSuv)
###get the same
genomymrMLMMSuv <-subgenomrMLMMtran[match(mymrMLMMSuv$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
genomymrMLMMSuv <- data.frame(genomymrMLMMSuv)+1
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMSuv <- Select.MAF(genomymrMLMMSuv)
genomymrMLMMSuv1 <-subgenomrMLMM[1:4][match(genomymrMLMMSuv$rn,subgenomrMLMM$rs., nomatch=0),]
colnames(genomymrMLMMSuv1)[which(names(genomymrMLMMSuv1) == "rs.")] <- "rn"
genomymrMLMMSuv <- plyr::join_all(list(genomymrMLMMSuv1,genomymrMLMMSuv),by="rn")
str(genomymrMLMMSuv)
###change the name in order to fit the software requirment
colnames(genomymrMLMMSuv)[which(names(genomymrMLMMSuv) == "rn")] <- "rs#"
colnames(genomymrMLMMSuv)[which(names(genomymrMLMMSuv) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMSuv, file = "mrMLMM2/subgenomrMLMMSuv3.csv", row.names = FALSE, na = "NA")

####including all of the methods: 

library("mrMLM")
mrMLM(fileGen="mrMLMM2\\subgenomrMLMMSuv3.csv",
      filePhe="mrMLMM2\\myYmrMlMMSuv.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/ResultSuvall3")

###using the one with imputed number in SNPs datasets
subgenomrMLMM <- read.csv("data/subgenomrMLMM.csv")
str(subgenomrMLMM)
subgenomrMLMM1 <- data.frame(subgenomrMLMM[,c(1,5:157)],row.names = 1)
subgenomrMLMMtran <- data.frame(t(subgenomrMLMM1))
str(subgenomrMLMMtran)
#mymrMLMMflo <- read.csv("data/myYmrMLMMflo.csv")# this one has no any missing values 

mymrMLMMSuv <- read.csv("mrMLMM2/myYmrMlMMSuv.csv")
str(mymrMLMMSuv)
###get the same
genomymrMLMMSuv <-subgenomrMLMMtran[match(mymrMLMMSuv$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
genomymrMLMMSuv <- data.frame(genomymrMLMMSuv)+1
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMSuv <- Select.MAF(genomymrMLMMSuv)
genomymrMLMMSuv1 <-subgenomrMLMM[1:4][match(genomymrMLMMSuv$rn,subgenomrMLMM$rs., nomatch=0),]
colnames(genomymrMLMMSuv1)[which(names(genomymrMLMMSuv1) == "rs.")] <- "rn"
genomymrMLMMSuv <- plyr::join_all(list(genomymrMLMMSuv1,genomymrMLMMSuv),by="rn")
str(genomymrMLMMSuv)
###change the name in order to fit the software requirment
colnames(genomymrMLMMSuv)[which(names(genomymrMLMMSuv) == "rn")] <- "rs#"
colnames(genomymrMLMMSuv)[which(names(genomymrMLMMSuv) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMSuv, file = "mrMLMM2/subgenomrMLMMSuv.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="mrMLMM2\\subgenomrMLMMSuv.csv",
      filePhe="mrMLMM2\\myYmrMlMMSuv.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/ResultSuvall1")

###this one for OWA3
###using the one with imputed number in SNPs datasets
subgenomrMLMM <- read.csv("data/subgenomrMLMM.csv")
str(subgenomrMLMM)
subgenomrMLMM1 <- data.frame(subgenomrMLMM[,c(1,5:157)],row.names = 1)
subgenomrMLMMtran <- data.frame(t(subgenomrMLMM1))
str(subgenomrMLMMtran)
#mymrMLMMflo <- read.csv("data/myYmrMLMMflo.csv")# this one has no any missing values 
mymrMLMMSuv <- read.csv("data/myYimputedSNP19OWA3.csv")
mymrMLMMSuv$Taxa <- make.names(mymrMLMMSuv$Taxa)

#mymrMLMMSuv <- read.csv("mrMLMM2/myYmrMlMMSuv.csv")
str(mymrMLMMSuv)
###get the same
genomymrMLMMSuv <-subgenomrMLMMtran[match(mymrMLMMSuv$Taxa,rownames(subgenomrMLMMtran), nomatch=0),]
genomymrMLMMSuv <- data.frame(genomymrMLMMSuv)+1
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMSuv <- Select.MAF(genomymrMLMMSuv)
genomymrMLMMSuv1 <-subgenomrMLMM[1:4][match(genomymrMLMMSuv$rn,subgenomrMLMM$rs., nomatch=0),]
colnames(genomymrMLMMSuv1)[which(names(genomymrMLMMSuv1) == "rs.")] <- "rn"
genomymrMLMMSuv <- plyr::join_all(list(genomymrMLMMSuv1,genomymrMLMMSuv),by="rn")
str(genomymrMLMMSuv)
###change the name in order to fit the software requirment
colnames(genomymrMLMMSuv)[which(names(genomymrMLMMSuv) == "rn")] <- "rs#"
colnames(genomymrMLMMSuv)[which(names(genomymrMLMMSuv) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMSuv, file = "mrMLMM2/subgenomrMLMMOWA3.csv", row.names = FALSE, na = "NA")
str(genomymrMLMMSuv)

colnames(mymrMLMMSuv)[which(names(mymrMLMMSuv) == "Taxa")] <- "<phenotype>"
colnames(mymrMLMMSuv)[which(names(mymrMLMMSuv) == "OWA")] <- "Trait1"
write.csv(mymrMLMMSuv, file="mrMLMM2/myYmrMlMMOWA3.csv",row.names = FALSE,na = "NA")
library("mrMLM")
mrMLM(fileGen="mrMLMM2\\subgenomrMLMMOWA3.csv",
      filePhe="mrMLMM2\\myYmrMlMMOWA3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/ResultOWA3all1")

