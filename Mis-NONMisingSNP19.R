############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY

myYmrMlMM.clum <- read.csv("data/myYmrMlMMclum.csv",na.strings = c("",".","NA"))


############# step2 
############# step2
############download the SNP data that Lindsay gave to me ################### 
############loading the genotype SNP data with missing value called datacomb3################
load("~/Documents/MiacnathusSNPinformation /160322filteredSNPs.RData")
########## inspect the genotype data
str(datacomb3) # 594 genotypes, 46,177 markers
dimnames(datacomb3)[[1]] # this is how the genotypes are named in the SNP data

write.csv(dimnames(datacomb3)[[1]], file = "~/Documents/MiacnathusSNPinformation /d3taxa1.csv",row.names = FALSE)
#############changing the rowname to the first columne##############################
datacomb3 <- data.frame(datacomb3)
Taxa <- rownames(datacomb3)
rownames(datacomb3) <- NULL
datacomb3 <- cbind(Taxa,datacomb3)

##############step3 select matched SNP 
##############step3 select matched SNP 
### checking Taxa if match in SNP and in phenotype data 
allblup$Taxa %in% datacomb3$Taxa
### submet the SNP data set match with the phenotype with total missing value
SNP <- datacomb3[match(allblup$Taxa, datacomb3$Taxa,nomatch=0),]

###change colomn to row for datacomb2
SNP <- data.frame(SNP,row.names=1)


##############step4 select suitable SNP and individual for the next analysis 
##############step4 select suitable SNP and individual for the next analysis 

########1)remove row with missing value more than 0.6, 
### Remove the row with 50% missing value (more than 50% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.50),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]
###SNPnmcol <- SNP[, -which(colMeans(is.na(SNP)) > 0)]


############select the MAF
#####select the MAF with SNP from the last 
SNPmeans <- colMeans(SNPnmcol, na.rm = TRUE)
#range(SNPmeans)
#levels(as.factor(SNPmeans))
# determine for which SNPs 0 is the major allele
ZeroIsMajor <- SNPmeans < 1
# set up an empty vector to indicate which SNPs to keep
SNPsToKeep <- logical(length(ZeroIsMajor))
# define a minimum number of accessions with minor allele
minWithMinor <- 3
# figure out which to keep
SNPsToKeep[ZeroIsMajor] <- colSums(SNPnmcol[, ZeroIsMajor] == 1 | SNPnmcol[, ZeroIsMajor] == 2, na.rm = TRUE) >= minWithMinor
SNPsToKeep[!ZeroIsMajor] <- colSums(SNPnmcol[, !ZeroIsMajor] == 1 | SNPnmcol[, !ZeroIsMajor] == 0, na.rm = TRUE) >= minWithMinor
# subset genotype matrix 
###using this mymat as new genotype for the next step. 
mymat <- SNPnmcol[, SNPsToKeep]


############step5 myGD
############step5 myGD
mymatg <-mymat+1 
mymatg <- data.frame(mymatg)
Taxa <- rownames(mymatg)
rownames(mymatg) <- NULL
mymatg <- cbind(Taxa,mymatg)
myGDorderg <-mymatg[order(mymatg$Taxa),]

#### write out the phenotype with correct TAXA for GAPIT analysis 
write.csv(myGDorderg, file = "~/Documents/whole traits/myGD127SNP19.csv",row.names = FALSE)

##############step6 myY
##############step6 myY
myYorderg <- allblup[match(myGDorderg$Taxa,allblup$Taxa, nomatch=0),]
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myYorderg)
#### write out the phenotype with correct TAXA
write.csv(myYorderg, file = "~/Documents/whole traits/myY127SNP19.csv", row.names = FALSE, na = "")

##############step7 myQ
##############step7 myQ
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("~/Documents/MiacnathusSNPinformation /myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####get the name of collumn 
names(myQ)
####change the Sample_name to Taxa
myQ$Taxa <- myQ$Sample_name
####check if the name match from both data
myYorderg$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ1 <- myQ[match(myYorderg$Taxa,myQ$Taxa, nomatch=0),]
myQorder <- myQ1[order(myQ1$Taxa),]
myQorder <- myQorder[,c(10,2:9)]
str(myQorder)
####write out the dataset
write.csv(myQorder, file = "~/Documents/whole traits/myQ127SNP19.csv",row.names = FALSE)

#############step8 myGM
#############step8 myGM
########## GM (Genetic mapping dataset)###############################
#load GM dataset get myGM and myQ
load("~/Documents/MiacnathusSNPinformation /161025forGAPIT.RData")
### select the GM dataset with the same ID with GD
#myGDorderg <- read.csv("~/Documents/whole traits/myGD127SNP19.csv")

n <- data.frame(myGDorderg,row.names = 1)
a <- data.frame(t(n))
myGM2 <- myGM[match(dimnames(a)[[1]],myGM$Name, nomatch=0),]
####output GM data with csv format
write.csv(myGM2, "~/Documents/whole traits/myGM127SNP19.csv",row.names = FALSE)

#### this one for rrBLUP
#### this one for rrBLUP
#### this one for rrBLUP
####rrBLUP analysis
####write out for the rrBLUP analysis 
#####changing the rowname to the first columne
mymat <- data.frame(mymat)
Taxa <- rownames(mymat)
rownames(mymat) <- NULL
mymat <- cbind(Taxa,mymat)
myGDorder <- mymat[order(mymat$Taxa),]
str(myGDorder)
#### write out the phenotype with correct TAXA
write.csv(myGDorder, file = "~/Documents/whole traits/myGD127SNPrrblup19.csv",row.names = FALSE)

