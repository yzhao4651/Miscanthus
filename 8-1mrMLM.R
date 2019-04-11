
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

###flowering traits

mymrMLMMflo <- read.csv("data/myYmrMLMMflo.csv")
str(mymrMLMMflo)
###get the same
genomymrMLMMflo <-subgenomrMLMMtran[match(mymrMLMMflo$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMflo <- Select.MAF(genomymrMLMMflo)
##get the mateched SNP 
genomymrMLMMflo <-subgenomrMLMM[match(genomymrMLMMflo$rn,subgenomrMLMM$rn., nomatch=0),]
###change the name in order to fit the software requirment
colnames(genomymrMLMMflo)[which(names(genomymrMLMMflo) == "rn")] <- "rs#"
colnames(genomymrMLMMflo)[which(names(genomymrMLMMflo) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMflo, file = "data/subgenomrMLMMflo.csv", row.names = FALSE, na = "NA")

###flowering traits

mymrMLMMflo <- read.csv("data/myYmrMLMMflo.csv")
str(mymrMLMMflo)
###get the same
genomymrMLMMflo <-subgenomrMLMMtran[match(mymrMLMMflo$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMflo <- Select.MAF(genomymrMLMMflo)
##get the mateched SNP 
genomymrMLMMflo <-subgenomrMLMM[match(genomymrMLMMflo$rn,subgenomrMLMM$rn., nomatch=0),]
###change the name in order to fit the software requirment
colnames(genomymrMLMMflo)[which(names(genomymrMLMMflo) == "rn")] <- "rs#"
colnames(genomymrMLMMflo)[which(names(genomymrMLMMflo) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMflo, file = "data/subgenomrMLMMflo.csv", row.names = FALSE, na = "NA")

###flowering number traits

mymrMLMMflonum <- read.csv("data/myYmrMLMMflonum.csv")
str(mymrMLMMflonum)
###get the same
genomymrMLMMflonum <-subgenomrMLMMtran[match(mymrMLMMflonum$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMflonum <- Select.MAF(genomymrMLMMflonum)
##get the mateched SNP 
genomymrMLMMflonum <-subgenomrMLMM[match(genomymrMLMMflonum$rn,subgenomrMLMM$rn., nomatch=0),]
###change the name in order to fit the software requirment
colnames(genomymrMLMMflonum)[which(names(genomymrMLMMflonum) == "rn")] <- "rs#"
colnames(genomymrMLMMflonum)[which(names(genomymrMLMMflonum) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMflonum, file = "data/subgenomrMLMMflonum.csv", row.names = FALSE, na = "NA")

###Survive
mymrMLMMSuv <- read.csv("data/myYmrMLMMSuv.csv")
str(mymrMLMMSuv)
###get the same
genomymrMLMMSuv <-subgenomrMLMMtran[match(mymrMLMMSuv$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
###select the MAF>0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMSuv <- Select.MAF(genomymrMLMMSuv)
##get the mateched SNP 
genomymrMLMMSuv <-subgenomrMLMM[match(genomymrMLMMSuv$rn,subgenomrMLMM$rn., nomatch=0),]
###change the name in order to fit the software requirment
colnames(genomymrMLMMSuv)[which(names(genomymrMLMMSuv) == "rn")] <- "rs#"
colnames(genomymrMLMMSuv)[which(names(genomymrMLMMSuv) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMSuv, file = "data/subgenomrMLMMSuv.csv", row.names = FALSE, na = "NA")





###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
###check the name of Taxa
myYmrMlMM$Taxa 
###get the column names of phenotype
names(myYmrMlMM)
###change the name of column
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<Phenotype>"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Surv")] <- "Trait1"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "HD_1")] <- "Trait2"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "FD_1")] <- "Trait3"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "HD_50.")] <- "Trait4"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "FD_50.")] <- "Trait5"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Yld_kg")] <- "Trait6"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "SDW_kg")] <- "Trait7"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "CmDW_g")] <- "Trait8"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Cml_cm")] <- "Trait9"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "CmD_BI_mm")] <- "Trait10"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "CmD_LI_mm")] <- "Trait11"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "CmN.")] <- "Trait12"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Bcirc_cm")] <- "Trait13"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Lg")] <- "Trait14"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "GS")] <- "Trait15"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "FD")] <- "Trait16"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "TFN.")] <- "Trait17"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "FNMain")] <- "Trait18"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "FNsmall")] <- "Trait19"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "CCirc_cm")] <- "Trait20"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "SRD")] <- "Trait21"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "ADD")] <- "Trait22"
#colnames(myYmrMlMM)[which(names(myYmrMlMM) == "fprin")] <- "Trait23"

### finally i got this function online, and it will save a lot of time to change the column name
## but I also keep the code above in order to find the original name for the traits
###using this function to change the name of this 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
###check the column name if changed 
names(myYmrMlMM)
###write out the dataset
write.csv(myYmrMlMM, file = "mrMLMM/myYmrMlMM.csv", row.names = FALSE, na = "NA")


install.packages("rlang")

library(mrMLM.GUI) 
mrMLM.GUI()

