
### import the data set for GWAS anaysis 
allblup <- read.csv("data/alltraits.csv",row.names = 1)

############# step2 Download the genotype data
############# step2 Download the genotype data
############download the SNP data that Lindsay gave to me ################### 
#### import the SNP values with already imputed missing value datacomb2
###trying to load in to GitHub, But this one is big, can not load from here. 
load("~/Documents/MiacnathusSNPinformation /160324EMimputedSNP_Msi.RData")
# load("~/DOE Msi study/yield manuscript/phenotypic and GWAS analysis/160324EMimputedSNP_Msi.RData") # Lindsay's version
### this data set in value part, change the value into data.frame
##check the name of that value:
names(myA.EM.Msi)
myA.EM.Msi$imputed[myA.EM.Msi$imputed > 1] <- 1
myA.EM.Msi$imputed[myA.EM.Msi$imputed < -1] <- -1
####change the data to data.frame
datacomb2 <- data.frame(myA.EM.Msi[["imputed"]] + 1) # convert to be 0, 1, 2 instead of -1, 0, 1
write.csv(dimnames(datacomb2)[[1]], file = "data/datacomb2.csv", row.names = FALSE)
#############changing the rowname to the first columne##############################
datacomb2 <- data.frame(datacomb2)
Taxa <- rownames(datacomb2)
rownames(datacomb2) <- NULL
datacomb2 <- cbind(Taxa,datacomb2)

##############step3 select Genotype data that matched phenotype 
##############step3 select Genotype data that matched phenotype  
### checking Taxa if match in SNP and in phenotype data 
allblup$Taxa %in% datacomb2$Taxa
allblup$Taxa[!allblup$Taxa %in% datacomb2$Taxa] # Msa, Mxg, a few Msi not genotyped
### submet the SNP data set match with the phenotype with total missing value
SNP <- datacomb2[match(allblup$Taxa, datacomb2$Taxa,nomatch=0),]
#######omit the missing values 
SNP <- na.omit(SNP)
###change colomn to row for datacomb2
SNP <- data.frame(SNP,row.names=1)


############step4 select the MAF larger than 1%
############step4 select the MAF larger than 1% 
##I think i asked this question before, but still get confused, I asked here again.
##Question: I still confused this one also about the MAF. Do you know where has formular for calculating MAF, 

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
# figure out which to keep based on how many individuals have minor allele
SNPsToKeep[ZeroIsMajor] <- colSums(SNP[, ZeroIsMajor] == 1 | SNP[, ZeroIsMajor] == 2, na.rm = TRUE) >= minWithMinor
SNPsToKeep[!ZeroIsMajor] <- colSums(SNP[, !ZeroIsMajor] == 1 | SNP[, !ZeroIsMajor] == 0, na.rm = TRUE) >= minWithMinor
# this is how I would filter based on MAF
MAF <- SNPmeans/2
MAF[!ZeroIsMajor] <- 1 - MAF[!ZeroIsMajor] # flip the frequencies if allele 2 is the common one
hist(MAF) # most near zero, max 0.5
SNPsToKeep <- SNPsToKeep & MAF > 0.01 # update SNPsToKeep; 3 ind. with minor allele AND MAF > 0.01. 35,279 SNPs.
# subset genotype matrix 
###using this mymat as new genotype for the next step. 
mymat <- SNP[, SNPsToKeep]

############step5 myGD: Genotype data for GAPIT and FarmCPU
############step5 myGD: Genotype data for GAPIT and FarmCPU
mymatg <-mymat
mymatg <- data.frame(mymatg)
Taxa <- rownames(mymatg)
rownames(mymatg) <- NULL
mymatg <- cbind(Taxa,mymatg)
myGDorderg <-mymatg[order(mymatg$Taxa),]
#### write out the phenotype with correct TAXA for GAPIT analysis 
write.csv(myGDorderg, file = "data/myGDimputedSNP19.csv",row.names = FALSE)

# free up RAM by removing large unused objects
rm(myA.EM.Msi, mymat, mymatg, SNP)

##############step6 myY: Phenotype data for GAPIT and FarmCPU 
##############step6 myY: Phenotype data for GAPIT and FarmCPU 
myYorderg <- allblup[match(myGDorderg$Taxa,allblup$Taxa, nomatch=0),]
myYorderg <- droplevels(myYorderg)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myYorderg)
#### write out the phenotype with correct TAXA
write.csv(myYorderg, file = "data/myYimputedSNP19.csv", row.names = FALSE, na = "")

##############step7 myQ: population structure
##############step7 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
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
write.csv(myQorder, file = "data/myQimputedSNP19.csv",row.names = FALSE)

#############step8 myGM:Genetic mapping information
#############step8 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)###############################
#load GM dataset get myGM and myQ
load("~/Documents/MiacnathusSNPinformation /161025forGAPIT.RData")
#load("~/DOE Msi study/yield manuscript/phenotypic and GWAS analysis/161025forGAPIT.RData") # on Lindsay's computer
### select the GM dataset with the same ID with GD
#myGDorderg <- read.csv("~/Documents/whole traits/myGDordfrrblup19.csv")
n <- data.frame(myGDorderg,row.names = 1)
a <- data.frame(t(n))
myGM2 <- myGM[match(dimnames(a)[[1]],myGM$Name, nomatch=0),]
# put the SNPs back in order by chromosome and position
identical(as.character(myGM2$Name), colnames(myGDorderg)[-1])
snporder <- order(myGM2$Chromosome, myGM2$Position)
myGM2 <- myGM2[snporder,]
myGDorderg <- myGDorderg[,c(1, snporder + 1)]
identical(as.character(myGM2$Name), colnames(myGDorderg)[-1])
####output GM data with csv format
write.csv(myGM2, "data/myGMimputedSNP19.csv",row.names = FALSE)

