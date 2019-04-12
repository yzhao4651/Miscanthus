############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
allblup <- read.csv("data/myYimputedSNP19.csv",na.strings = c("",".","NA"))
allblup$Taxa <- make.names(allblup$Taxa)
allblupflo <- allblup[,c(1,14:16,20:39)]
str(allblupflo)
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
SNP <- subsetNONmisingSNP.GAPIT(datacomb3,allblupflo)

##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 35, 
### Remove the row with 35% missing value (more than 35% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.34),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.142.2942<- Select.MAF(SNPnmcol)


#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)

### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.142.2942))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.142.2942 <- getmyGD(myGM,myGD.142.2942)
myGM.142.2942 <- getmyGM(myGM,myGD.142.2942)
###check if they are shared the same name in the same order
identical(as.character(myGM.142.2942$Name), colnames(myGD.142.2942))
write.csv(myGM.142.2942, file = "data/myGMflo.142.2942.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.142.2942 <- FirstColumn(myGD.142.2942)
myGD.142.2942$Taxa <- make.names(myGD.142.2942$Taxa)
write.csv(myGD.142.2942, file = "data/myGDflo.142.2942.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 

myY.142.2942 <- allblupflo[match(myGD.142.2942$Taxa,allblupflo$Taxa, nomatch=0),]
myY.142.2942 <- droplevels(myY.142.2942)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.142.2942)
#### write out the phenotype with correct TAXA
write.csv(myY.142.2942, file = "data/myYflo.142.2942.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myQ$Taxa <- make.names(myQ$Taxa)
myY.142.2942$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.142.2942 <- myQ[match(myY.142.2942$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.142.2942, file = "data/myQflo.142.2942.csv",row.names = FALSE)

##########change 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 35, 
### Remove the row with 35% missing value (more than 35% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.317),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.132.3814<- Select.MAF(SNPnmcol)


#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)

### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.132.3814))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.132.3814 <- getmyGD(myGM,myGD.132.3814)
myGM.132.3814 <- getmyGM(myGM,myGD.132.3814)
###check if they are shared the same name in the same order
identical(as.character(myGM.132.3814$Name), colnames(myGD.132.3814))
write.csv(myGM.132.3814, file = "data/myGMflo.132.3814.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.132.3814 <- FirstColumn(myGD.132.3814)
myGD.132.3814$Taxa <- make.names(myGD.132.3814$Taxa)
write.csv(myGD.132.3814, file = "data/myGDflo.132.3814.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 

myY.132.3814 <- allblupflo[match(myGD.132.3814$Taxa,allblupflo$Taxa, nomatch=0),]
myY.132.3814 <- droplevels(myY.132.3814)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.132.3814)
#### write out the phenotype with correct TAXA
write.csv(myY.132.3814, file = "data/myYflo.132.3814.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myQ$Taxa <- make.names(myQ$Taxa)
myY.132.3814$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.132.3814 <- myQ[match(myY.132.3814$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.132.3814, file = "data/myQflo.132.3814.csv",row.names = FALSE)


##########change 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 35, 
### Remove the row with 35% missing value (more than 35% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.26),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.112.4876<- Select.MAF(SNPnmcol)


#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.112.4876))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.112.4876 <- getmyGD(myGM,myGD.112.4876)
myGM.112.4876 <- getmyGM(myGM,myGD.112.4876)
###check if they are shared the same name in the same order
identical(as.character(myGM.112.4876$Name), colnames(myGD.112.4876))
write.csv(myGM.112.4876, file = "data/myGMflo.112.4876.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.112.4876 <- FirstColumn(myGD.112.4876)
myGD.112.4876$Taxa <- make.names(myGD.112.4876$Taxa)
write.csv(myGD.112.4876, file = "data/myGDflo.112.4876.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 

myY.112.4876 <- allblupflo[match(myGD.112.4876$Taxa,allblupflo$Taxa, nomatch=0),]
myY.112.4876 <- droplevels(myY.112.4876)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.112.4876)
#### write out the phenotype with correct TAXA
write.csv(myY.112.4876, file = "data/myYflo.112.4876.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myQ$Taxa <- make.names(myQ$Taxa)
myY.112.4876$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.112.4876 <- myQ[match(myY.112.4876$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.112.4876, file = "data/myQflo.112.4876.csv",row.names = FALSE)

##########change 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 35, 
### Remove the row with 35% missing value (more than 35% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.25),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.103.5305<- Select.MAF(SNPnmcol)


#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.103.5305))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.103.5305 <- getmyGD(myGM,myGD.103.5305)
myGM.103.5305 <- getmyGM(myGM,myGD.103.5305)
###check if they are shared the same name in the same order
identical(as.character(myGM.103.5305$Name), colnames(myGD.103.5305))
write.csv(myGM.103.5305, file = "data/myGMflo.112.4876.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.103.5305 <- FirstColumn(myGD.103.5305)
myGD.103.5305$Taxa <- make.names(myGD.103.5305$Taxa)
write.csv(myGD.103.5305, file = "data/myGDflo.103.5305.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 

myY.103.5305 <- allblupflo[match(myGD.103.5305$Taxa,allblupflo$Taxa, nomatch=0),]
myY.103.5305 <- droplevels(myY.103.5305)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.103.5305)
#### write out the phenotype with correct TAXA
write.csv(myY.103.5305, file = "data/myYflo.103.5305.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myQ$Taxa <- make.names(myQ$Taxa)
myY.103.5305$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.103.5305 <- myQ[match(myY.103.5305$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.103.5305, file = "data/myQflo.103.5305.csv",row.names = FALSE)


