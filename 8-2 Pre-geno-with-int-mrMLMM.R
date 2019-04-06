####import all of the genotype data
subgenomrMLMM <- read.csv("data/subgenomrMLMM.csv")
str(subgenomrMLMM)
subgenomrMLMM1 <- data.frame(subgenomrMLMM[,c(1,5:157)],row.names = 1)
str(subgenomrMLMM1)
###seperate two parts
subgenomrMLMMtran <- data.frame(t(subgenomrMLMM1))
subgenomrMLMMtran[subgenomrMLMMtran < 1 & subgenomrMLMMtran > 0.5] <- 1
subgenomrMLMMtran[subgenomrMLMMtran < 0.5 & subgenomrMLMMtran > 0] <- 0
subgenomrMLMMtran[subgenomrMLMMtran < -0.5 & subgenomrMLMMtran > -1] <- -1
subgenomrMLMMtran[subgenomrMLMMtran < 0 & subgenomrMLMMtran > -0.5] <- 0

rownames(subgenomrMLMMtran)
str(subgenomrMLMMtran)
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
str(genomymrMLMMculm)

##get the mateched SNP 
genomymrMLMMculm1 <-subgenomrMLMM[1:4][match(genomymrMLMMculm$rn,subgenomrMLMM$rs., nomatch=0),]

colnames(genomymrMLMMculm1)[which(names(genomymrMLMMculm1) == "rs.")] <- "rn"

genomymrMLMMculm <- plyr::join_all(list(genomymrMLMMculm1,genomymrMLMMculm),by="rn")

str(genomymrMLMMculm)
###change the name in order to fit the software requirment
colnames(genomymrMLMMculm)[which(names(genomymrMLMMculm) == "rn")] <- "rs#"
colnames(genomymrMLMMculm)[which(names(genomymrMLMMculm) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMculm, file = "mrMLMM2/subgenomrMLMMculm2.csv", row.names = FALSE, na = "NA")

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


###Survive
mymrMLMMSuv <- read.csv("mrMLMM2/myYmrMlMMSuv2.csv")
str(mymrMLMMSuv)
###get the same
genomymrMLMMSuv <-subgenomrMLMMtran[match(mymrMLMMSuv$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
genomymrMLMMSuv  <- data.frame(genomymrMLMMSuv)+1
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMSuv <- Select.MAF(genomymrMLMMSuv)
##get the mateched SNP 
genomymrMLMMSuv1 <-subgenomrMLMM[1:4][match(genomymrMLMMSuv$rn,subgenomrMLMM$rs., nomatch=0),]
str(genomymrMLMMSuv1)
colnames(genomymrMLMMSuv1)[which(names(genomymrMLMMSuv1) == "rs.")] <- "rn"
genomymrMLMMSuv <- plyr::join_all(list(genomymrMLMMSuv1,genomymrMLMMSuv),by="rn")
str(genomymrMLMMSuv)
###change the name in order to fit the software requirment
colnames(genomymrMLMMSuv)[which(names(genomymrMLMMSuv) == "rn")] <- "rs#"
colnames(genomymrMLMMSuv)[which(names(genomymrMLMMSuv) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMSuv, file = "mrMLMM2/subgenomrMLMMSuv2.csv", row.names = FALSE, na = "NA")

install.packages("rlang")

library(mrMLM.GUI) 
mrMLM.GUI()

