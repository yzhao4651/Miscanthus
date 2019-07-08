############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
allblup <- read.csv("data/alltraits.csv", na.strings = c("",".","NA"),row.names = 1)
str(allblup)
###import the SRD trait
allblup.1 <- read.csv("data/alltraitsSRD.csv", na.strings = c("",".","NA"),row.names = 1)
allblupSRD <- allblup.1[,c(1,3)]
str(allblupSRD)
allblup.2 <- read.csv("data/alltraitsOWA.csv", na.strings = c("",".","NA"),row.names = 1)
allblupOWA <- allblup.2[,c(1,3)]
str(allblupOWA)
##combine both 
alltraits <- plyr::join_all(list(allblup,allblupOWA,allblupSRD),by="Taxa")
str(alltraits)
alltraitsorder <- alltraits[,c(1,40,3:17,41,19:39)]
str(alltraitsorder)
write.csv(alltraitsorder,file="data/alltraits1.csv",row.names = FALSE)

###import the whole data again 
allblup <- read.csv("data/alltraits1.csv", na.strings = c("",".","NA"))
str(allblup)
allblupClu<- allblup[,c(1, 3:13,17:19)]
str(allblupClu)
allblupClu<- allblupClu[-which(rowSums(is.na(allblupClu[1:15])) ==14),]
str(allblupClu)
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
SNP <- subsetNONmisingSNP.GAPIT(datacomb3,allblupClu)

##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 30, 
### Remove the row with 30% missing value (more than 30% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.3),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.106.4202 <- Select.MAF(SNPnmcol)

#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)

### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.106.4202))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.106.4202 <- getmyGD(myGM,myGD.106.4202)

myGM.106.4202 <- getmyGM(myGM,myGD.106.4202)

###check if they are shared the same name in the same order
identical(as.character(myGM.106.4202$Name), colnames(myGD.106.4202))
write.csv(myGM.106.4202, file = "data/myGMC.106.4202.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.106.4202 <- FirstColumn(myGD.106.4202)

write.csv(myGD.106.4202, file = "data/myGDC.106.4202.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.106.4202 <- allblupClu[match(myGD.106.4202$Taxa,allblupClu$Taxa, nomatch=0),]
myY.106.4202 <- droplevels(myY.106.4202)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])

#### write out the phenotype with correct TAXA
write.csv(myY.106.4202, file = "data/myYC.106.4202.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.106.4202$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.106.4202 <- myQ[match(myY.106.4202$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.106.4202, file = "data/myQC.106.4202.csv",row.names = FALSE)

###0.351
###0.351
###0.351
##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 351, 
### Remove the row with 351% missing value (more than 351% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.336),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.116.3293 <- Select.MAF(SNPnmcol)

#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)

### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.116.3293))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.116.3293 <- getmyGD(myGM,myGD.116.3293)

myGM.116.3293 <- getmyGM(myGM,myGD.116.3293)

###check if they are shared the same name in the same order
identical(as.character(myGM.116.3293$Name), colnames(myGD.116.3293))
write.csv(myGM.116.3293, file = "data/myGMC.116.3293.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.116.3293 <- FirstColumn(myGD.116.3293)

write.csv(myGD.116.3293, file = "data/myGDC.116.3293.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.116.3293 <- allblupClu[match(myGD.116.3293$Taxa,allblupClu$Taxa, nomatch=0),]
myY.116.3293 <- droplevels(myY.116.3293)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])

#### write out the phenotype with correct TAXA
write.csv(myY.116.3293, file = "data/myYC.116.3293.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.116.3293$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.116.3293 <- myQ[match(myY.116.3293$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.116.3293, file = "data/myQC.116.3293.csv",row.names = FALSE)

###0.42
###0.42
###0.42
### Remove the row with 42% missing value (more than 42% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.42),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.125.2562 <- Select.MAF(SNPnmcol)

#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)

### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.125.2562))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.125.2562 <- getmyGD(myGM,myGD.125.2562)
myGM.125.2562 <- getmyGM(myGM,myGD.125.2562)

###check if they are shared the same name in the same order
identical(as.character(myGM.125.2562$Name), colnames(myGD.125.2562))
write.csv(myGM.125.2562, file = "data/myGMC.125.2562.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.125.2562 <- FirstColumn(myGD.125.2562)
write.csv(myGD.125.2562, file = "data/myGDC.125.2562.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.125.2562 <- allblupClu[match(myGD.125.2562$Taxa,allblupClu$Taxa, nomatch=0),]
myY.125.2562 <- droplevels(myY.125.2562)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])

#### write out the phenotype with correct TAXA
write.csv(myY.125.2562, file = "data/myYC.125.2562.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.125.2562$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.125.2562 <- myQ[match(myY.125.2562$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.125.2562, file = "data/myQC.125.2562.csv",row.names = FALSE)


###0.63
###0.63
###0.63
### Remove the row with 63% missing value (more than 63% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.63),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]
############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.130.2001 <- Select.MAF(SNPnmcol)
#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.130.2001))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.130.2001 <- getmyGD(myGM,myGD.130.2001)
myGM.130.2001 <- getmyGM(myGM,myGD.130.2001)
###check if they are shared the same name in the same order
identical(as.character(myGM.130.2001$Name), colnames(myGD.130.2001))
write.csv(myGM.130.2001, file = "data/myGMC.130.2001.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.130.2001 <- FirstColumn(myGD.130.2001)
write.csv(myGD.130.2001, file = "data/myGDC.130.2001.csv",row.names = FALSE)
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.130.2001 <- allblupClu[match(myGD.130.2001$Taxa,allblupClu$Taxa, nomatch=0),]
myY.130.2001 <- droplevels(myY.130.2001)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
#### write out the phenotype with correct TAXA
write.csv(myY.130.2001, file = "data/myYC.130.2001.csv", row.names = FALSE, na = "")
##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.130.2001$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.130.2001 <- myQ[match(myY.130.2001$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.130.2001, file = "data/myQC.130.2001.csv",row.names = FALSE)






