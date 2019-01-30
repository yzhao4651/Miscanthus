############step 1 emerge all the phenotype data 
############step 1 emerge all the phenotype data 
### data 1 with entry and accession that matched the SNP 
Altable <- read.csv("~/Documents/MiacnathusSNPinformation /Altable19.csv", header = TRUE, na="")
Altable <- data.frame(Altable)
str(Altable)
###data 2 with flowering trait with pc1
ftpc1 <- read.csv("~/Documents/whole traits/fprin1.csv",stringsAsFactors = FALSE, header = TRUE)
ftpc1$Entry <- ftpc1$X1
ftpc1$flopc1 <- ftpc1$X2
ftpc1 <- ftpc1[,c(4:5)]
str(ftpc1)
###data3 with all of other traits 
alltraitsblup <- read.csv("~/Documents/whole traits/alltraitsblup.csv",stringsAsFactors = FALSE, header = TRUE)
is.na(alltraitsblup) <- !alltraitsblup
alltraitsblup$Entry <- alltraitsblup$X
str(alltraitsblup)
####combine all data tegather  
allblup <- plyr::join_all(list(alltraitsblup,ftpc1,Altable[1:2]), by='Entry')
allblup$Taxa <- allblup$Accession
allblup <- allblup[,c(16,2:12,14)]
allblup <- subset(allblup,!(is.na(allblup$Taxa)))
str(allblup)
write.csv(ftdaysblup, file = "~/Documents/whole traits/allblup.csv",row.names = T)

############# step2 
############# step2
############download the SNP data that Lindsay gave to me ################### 
#### import the SNP values with already imputed missing value datacomb2
load("~/Documents/MiacnathusSNPinformation /160324EMimputedSNP_Msi.RData")
### this data set in value part, change the value into data.frame
##check the name of that value:
names(myA.EM.Msi)
myA.EM.Msi$imputed[myA.EM.Msi$imputed > 1] <- 1
myA.EM.Msi$imputed[myA.EM.Msi$imputed < -1] <- -1
####change the data to data.frame
datacomb2 <- data.frame(myA.EM.Msi[["imputed"]])
write.csv(dimnames(datacomb2)[[1]], file = "~/Documents/MiacnathusSNPinformation /d2taxa.csv", row.names = FALSE)
#############changing the rowname to the first columne##############################
datacomb2 <- data.frame(datacomb2)
Taxa <- rownames(datacomb2)
rownames(datacomb2) <- NULL
datacomb2 <- cbind(Taxa,datacomb2)

##############step3 select matched SNP 
##############step3 select matched SNP 
### checking Taxa if match in SNP and in phenotype data 
allblup$Taxa %in% datacomb2$Taxa
### submet the SNP data set match with the phenotype with total missing value
SNP <- datacomb2[match(allblup$Taxa, datacomb2$Taxa,nomatch=0),]

#######omit the missing values 
SNP <- na.omit(SNP)
###change colomn to row for datacomb2
SNP <- data.frame(SNP,row.names=1)

############step4 select the MAF 
############step4 select the MAF
#####select the MAF with SNP from the last 
SNPmeans <- colMeans(SNP, na.rm = TRUE)
#range(SNPmeans)
#levels(as.factor(SNPmeans))
# determine for which SNPs 0 is the major allele
ZeroIsMajor <- SNPmeans < 1
# set up an empty vector to indicate which SNPs to keep
SNPsToKeep <- logical(length(ZeroIsMajor))
# define a minimum number of accessions with minor allele
minWithMinor <- 3
# figure out which to keep
SNPsToKeep[ZeroIsMajor] <- colSums(SNP[, ZeroIsMajor] == 1 | SNP[, ZeroIsMajor] == 2, na.rm = TRUE) >= minWithMinor
SNPsToKeep[!ZeroIsMajor] <- colSums(SNP[, !ZeroIsMajor] == 1 | SNP[, !ZeroIsMajor] == 0, na.rm = TRUE) >= minWithMinor
# subset genotype matrix 
###using this mymat as new genotype for the next step. 
mymat <- SNP[, SNPsToKeep]

############step5 myGD
############step5 myGD
mymatg <-mymat+1 
mymatg <- data.frame(mymatg)
Taxa <- rownames(mymatg)
rownames(mymatg) <- NULL
mymatg <- cbind(Taxa,mymatg)
myGDorderg <-mymatg[order(mymatg$Taxa),]

#### write out the phenotype with correct TAXA for GAPIT analysis 
write.csv(myGDorderg, file = "~/Documents/whole traits/myGDimputedSNP19.csv",row.names = FALSE)

##############step6 myY
##############step6 myY
myYorderg <- allblup[match(myGDorderg$Taxa,allblup$Taxa, nomatch=0),]
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myYorderg)
#### write out the phenotype with correct TAXA
write.csv(myYorder, file = "~/Documents/whole traits/myYimputedSNP19.csv", row.names = FALSE, na = "")

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
write.csv(myQorder, file = "~/Documents/whole traits/myQimputedSNP19.csv",row.names = FALSE)

#############step8 myGM
#############step8 myGM
########## GM (Genetic mapping dataset)###############################
#load GM dataset get myGM and myQ
load("~/Documents/MiacnathusSNPinformation /161025forGAPIT.RData")
### select the GM dataset with the same ID with GD
#myGDorderg <- read.csv("~/Documents/whole traits/myGDordfrrblup19.csv")

n <- data.frame(myGDorderg,row.names = 1)
a <- data.frame(t(n))
myGM2 <- myGM[match(dimnames(a)[[1]],myGM$Name, nomatch=0),]
####output GM data with csv format
write.csv(myGM2, "~/Documents/whole traits/myGMimputedSNP19.csv",row.names = FALSE)

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
write.csv(myGDorder, file = "~/Documents/whole traits/myGDimputedSNPrrblup19.csv",row.names = FALSE)

