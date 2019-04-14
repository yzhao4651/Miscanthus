###change some number in the genotype to int
subgenomrMLMM <- read.csv("data/subgenomrMLMM.csv")
subgenomrMLMM1 <- data.frame(subgenomrMLMM[,c(1,5:157)],row.names = 1)
###seperate two parts
subgenomrMLMMtran <- data.frame(t(subgenomrMLMM1))
subgenomrMLMMtran[subgenomrMLMMtran < 1 & subgenomrMLMMtran > 0.5] <- 1
subgenomrMLMMtran[subgenomrMLMMtran < 0.5 & subgenomrMLMMtran > 0] <- 0
subgenomrMLMMtran[subgenomrMLMMtran < -0.5 & subgenomrMLMMtran > -1] <- -1
subgenomrMLMMtran[subgenomrMLMMtran < 0 & subgenomrMLMMtran > -0.5] <- 0
###flowering traits
###flowering traits
mymrMLMMflo <- read.csv("mrMLMM2/myYmrMlMMflo3.csv")
str(mymrMLMMflo)
###get the same
genomymrMLMMflo <-subgenomrMLMMtran[match(mymrMLMMflo$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
genomymrMLMMflo <- data.frame(genomymrMLMMflo)+1
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMflo <- Select.MAF(genomymrMLMMflo)
genomymrMLMMflo1 <-subgenomrMLMM[1:4][match(genomymrMLMMflo$rn,subgenomrMLMM$rs., nomatch=0),]
colnames(genomymrMLMMflo1)[which(names(genomymrMLMMflo1) == "rs.")] <- "rn"
genomymrMLMMflo <- plyr::join_all(list(genomymrMLMMflo1,genomymrMLMMflo),by="rn")
str(genomymrMLMMflo)
###change the name in order to fit the software requirment
colnames(genomymrMLMMflo)[which(names(genomymrMLMMflo) == "rn")] <- "rs#"
colnames(genomymrMLMMflo)[which(names(genomymrMLMMflo) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMflo, file = "mrMLMM2/subgenomrMLMMflo3.csv", row.names = FALSE, na = "NA")

####including all of the methods: 

library("mrMLM")
mrMLM(fileGen="mrMLMM2\\subgenomrMLMMflo3.csv",
      filePhe="mrMLMM2\\myYmrMlMMflo3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/Resultfloall3")

###using the one with imputed number in SNPs datasets
subgenomrMLMM <- read.csv("data/subgenomrMLMM.csv")
str(subgenomrMLMM)
subgenomrMLMM1 <- data.frame(subgenomrMLMM[,c(1,5:157)],row.names = 1)
subgenomrMLMMtran <- data.frame(t(subgenomrMLMM1))
str(subgenomrMLMMtran)
#mymrMLMMflo <- read.csv("data/myYmrMLMMflo.csv")# this one has no any missing values 
mymrMLMMflo <- read.csv("mrMLMM2/myYmrMLMMflo3.csv")# this one has two more column and then has missing values 
str(mymrMLMMflo)
###get the same
genomymrMLMMflo <-subgenomrMLMMtran[match(mymrMLMMflo$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
genomymrMLMMflo <- data.frame(genomymrMLMMflo)+1
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMflo <- Select.MAF(genomymrMLMMflo)
genomymrMLMMflo1 <-subgenomrMLMM[1:4][match(genomymrMLMMflo$rn,subgenomrMLMM$rs., nomatch=0),]
colnames(genomymrMLMMflo1)[which(names(genomymrMLMMflo1) == "rs.")] <- "rn"
genomymrMLMMflo <- plyr::join_all(list(genomymrMLMMflo1,genomymrMLMMflo),by="rn")
str(genomymrMLMMflo)
###change the name in order to fit the software requirment
colnames(genomymrMLMMflo)[which(names(genomymrMLMMflo) == "rn")] <- "rs#"
colnames(genomymrMLMMflo)[which(names(genomymrMLMMflo) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMflo, file = "mrMLMM2/subgenomrMLMMflo.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="mrMLMM2\\subgenomrMLMMflo.csv",
      filePhe="mrMLMM2\\myYmrMlMMflo3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/Resultfloall1")
