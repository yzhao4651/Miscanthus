############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
allblup <- read.csv("data/myYimputedSNP19.csv",na.strings = c("",".","NA"))
allblupCulm <- allblup[,c(1,3:13,17:19)]
str(allblupCulm)
###delete all the rows with missing values 
allblupCulm<- allblupCulm[-which(rowSums(is.na(allblupCulm)) ==14),]
allblupCulm <- droplevels(allblupCulm)

str(allblupCulm)
############# step2 
############# step2
############download the SNP data that Lindsay gave to me ################### 
############loading the genotype SNP data with missing value called datacomb3####
load("~/Documents/MiacnathusSNPinformation /160322filteredSNPs.RData")
load("C:/Users/Admin/Desktop/New folder/miscanthus study-1/Misthcanthus CCA data analysis/160322filteredSNPs.RData")
########## inspect the genotype data
#str(datacomb3) # 594 genotypes, 46,177 markers
#dimnames(datacomb3)[[1]] # this is how the genotypes are named in the SNP data
#############changing the rowname to the first columne###########################
source("Function/subsetNONmissingSNP.GAPIT.R")
SNP <- subsetNONmisingSNP.GAPIT(datacomb3,allblupCulm)

##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 35, 
### Remove the row with 35% missing value (more than 35% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.42),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.123.2670<- Select.MAF(SNPnmcol)


#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.123.2670))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.123.2670 <- getmyGD(myGM,myGD.123.2670)
myGM.123.2670 <- getmyGM(myGM,myGD.123.2670)
###check if they are shared the same name in the same order
identical(as.character(myGM.123.2670$Name), colnames(myGD.123.2670))
write.csv(myGM.123.2670, file = "data/myGMCulm.123.2670.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.123.2670 <- FirstColumn(myGD.123.2670)
myGD.123.2670$Taxa <- make.names(myGD.123.2670$Taxa)
write.csv(myGD.123.2670, file = "data/myGDCulm.123.2670.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
allblupCulm$Taxa <- make.names(allblupCulm$Taxa)
myGD.123.2670$Taxa %in% allblupCulm$Taxa
myY.123.2670 <- allblupCulm[match(myGD.123.2670$Taxa,allblupCulm$Taxa, nomatch=0),]
myY.123.2670 <- droplevels(myY.123.2670)
######### order the Taxa

#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.123.2670)
#### write out the phenotype with correct TAXA
write.csv(myY.123.2670, file = "data/myYCulm.123.2670.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myQ$Taxa <- make.names(myQ$Taxa)
myY.123.2670$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.123.2670 <- myQ[match(myY.123.2670$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.123.2670, file = "data/myQCulm.123.2670.csv",row.names = FALSE)




##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 35, 
### Remove the row with 35% missing value (more than 35% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.33),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]
source("Function/SelectMAF-GAPIT.R")
myGD.113.3489<- Select.MAF(SNPnmcol)
#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.113.3489))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.113.3489<- getmyGD(myGM,myGD.113.3489)
myGM.113.3489 <- getmyGM(myGM,myGD.113.3489)
###check if they are shared the same name in the same order
identical(as.character(myGM.113.3489$Name), colnames(myGD.113.3489))
write.csv(myGM.113.3489, file = "data/myGMCulm.113.3489.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.113.3489 <- FirstColumn(myGD.113.3489)
myGD.113.3489$Taxa <- make.names(myGD.113.3489$Taxa)
write.csv(myGD.113.3489, file = "data/myGDCulm.113.3489.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
allblupCulm$Taxa <- make.names(allblupCulm$Taxa)
myY.113.3489 <- allblupCulm[match(myGD.113.3489$Taxa,allblupCulm$Taxa, nomatch=0),]
myY.113.3489 <- droplevels(myY.113.3489)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.113.3489)
#### write out the phenotype with correct TAXA
write.csv(myY.113.3489, file = "data/myYCulm.113.3489.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myQ$Taxa <- make.names(myQ$Taxa)
myY.113.3489$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.113.3489 <- myQ[match(myY.113.3489$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.113.3489, file = "data/myQCulm.113.3489.csv",row.names = FALSE)


##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 35, 
### Remove the row with 35% missing value (more than 35% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.31),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]
source("Function/SelectMAF-GAPIT.R")
myGD.106.4220<- Select.MAF(SNPnmcol)
#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.106.4220))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.106.4220<- getmyGD(myGM,myGD.106.4220)
myGM.106.4220 <- getmyGM(myGM,myGD.106.4220)
###check if they are shared the same name in the same order
identical(as.character(myGM.106.4220$Name), colnames(myGD.106.4220))
write.csv(myGM.106.4220, file = "data/myGMCulm.106.4220.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.106.4220 <- FirstColumn(myGD.106.4220)
myGD.106.4220$Taxa <- make.names(myGD.106.4220$Taxa)
write.csv(myGD.106.4220, file = "data/myGDCulm.106.4220.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
allblupCulm$Taxa <- make.names(allblupCulm$Taxa)
myY.106.4220 <- allblupCulm[match(myGD.106.4220$Taxa,allblupCulm$Taxa, nomatch=0),]
myY.106.4220 <- droplevels(myY.106.4220)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.106.4220)
#### write out the phenotype with correct TAXA
write.csv(myY.106.4220, file = "data/myYCulm.106.4220.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myQ$Taxa <- make.names(myQ$Taxa)
myY.106.4220$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.106.4220 <- myQ[match(myY.106.4220$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.106.4220, file = "data/myQCulm.106.4220.csv",row.names = FALSE)



