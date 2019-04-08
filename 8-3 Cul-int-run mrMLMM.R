####import all of the genotype data 
####import all of the genotype data
subgenomrMLMM <- read.csv("data/subgenomrMLMM.csv")
###seperate two parts
subgenomrMLMM1 <- data.frame(subgenomrMLMM[,c(1,5:157)],row.names = 1)
###change the some number to the int
subgenomrMLMMtran <- data.frame(t(subgenomrMLMM1))
subgenomrMLMMtran[subgenomrMLMMtran < 1 & subgenomrMLMMtran > 0.5] <- 1
subgenomrMLMMtran[subgenomrMLMMtran < 0.5 & subgenomrMLMMtran > 0] <- 0
subgenomrMLMMtran[subgenomrMLMMtran < -0.5 & subgenomrMLMMtran > -1] <- -1
subgenomrMLMMtran[subgenomrMLMMtran < 0 & subgenomrMLMMtran > -0.5] <- 0

####import all of the phenotyp data
###Clum trait
mymrMLMMculm <- read.csv("mrMLMM2/myYmrMlMMculm3.csv")
str(mymrMLMMculm)
###get the same
genomymrMLMMculm <-subgenomrMLMMtran[match(mymrMLMMculm$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
genomymrMLMMculm <- data.frame(genomymrMLMMculm)+1
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMculm <- Select.MAF(genomymrMLMMculm)
##get the mateched SNP 
genomymrMLMMculm1 <-subgenomrMLMM[1:4][match(genomymrMLMMculm$rn,subgenomrMLMM$rs., nomatch=0),]
colnames(genomymrMLMMculm1)[which(names(genomymrMLMMculm1) == "rs.")] <- "rn"
genomymrMLMMculm <- plyr::join_all(list(genomymrMLMMculm1,genomymrMLMMculm),by="rn")
str(genomymrMLMMculm)
###change the name in order to fit the software requirment
colnames(genomymrMLMMculm)[which(names(genomymrMLMMculm) == "rn")] <- "rs#"
colnames(genomymrMLMMculm)[which(names(genomymrMLMMculm) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMculm, file = "mrMLMM2/subgenomrMLMMculm3.csv", row.names = FALSE, na = "NA")

####including all of the methods: 

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/Resultculmall3")

###this one has the SNPS with imputed number not int
###this one has the SNPS with imputed number not int
####import all of the genotype data
subgenomrMLMM <- read.csv("data/subgenomrMLMM.csv")
subgenomrMLMM1 <- data.frame(subgenomrMLMM[,c(1,5:157)],row.names = 1)
###seperate two parts
subgenomrMLMMtran <- data.frame(t(subgenomrMLMM1))
####import all of the phenotyp data
###Clum trait
mymrMLMMculm <- read.csv("mrMLMM2/myYmrMLMMculm3.csv")

###get the same
genomymrMLMMculm <-subgenomrMLMMtran[match(mymrMLMMculm$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMculm <- Select.MAF(genomymrMLMMculm)
##get the mateched SNP 
genomymrMLMMculm <-subgenomrMLMM[match(genomymrMLMMculm$rn,subgenomrMLMM$rn., nomatch=0),]
###change the name in order to fit the software requirment
colnames(genomymrMLMMculm)[which(names(genomymrMLMMculm) == "rn")] <- "rs#"
colnames(genomymrMLMMculm)[which(names(genomymrMLMMculm) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMculm, file = "data/subgenomrMLMMculm.csv", row.names = FALSE, na = "NA")
subgenomrMLMMculm <- read.csv("data/subgenomrMLMMculm.csv")
str(subgenomrMLMMculm)

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/Resultculmall1")