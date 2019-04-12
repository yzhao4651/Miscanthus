############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
allblup <- read.csv("data/myYimputedSNP19.csv",na.strings = c("",".","NA"))
allblupflo <- allblup[,c(1,14:16,20:39)]
str(allblupflo)
###delete all the rows with missing values 
allblupflo<- allblupflo[-which(rowSums(is.na(allblupflo)) ==23),]

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
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.4215),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]
############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.115.3154<- Select.MAF(SNPnmcol)


#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.115.3154))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.115.3154 <- getmyGD(myGM,myGD.115.3154)
myGM.115.3154 <- getmyGM(myGM,myGD.115.3154)
###check if they are shared the same name in the same order
identical(as.character(myGM.115.3154$Name), colnames(myGD.115.3154))
write.csv(myGM.115.3154, file = "data/myGMflo.142.2942.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.115.3154 <- FirstColumn(myGD.115.3154)
myGD.115.3154$Taxa <- make.names(myGD.115.3154$Taxa)
write.csv(myGD.115.3154, file = "data/myGDflo.115.3154.csv",row.names = FALSE)

##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
allblupflo$Taxa <- make.names(allblupflo$Taxa)
myY.115.3154 <- allblupflo[match(myGD.115.3154$Taxa,allblupflo$Taxa, nomatch=0),]
myY.115.3154 <- droplevels(myY.115.3154)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.115.3154)
#### write out the phenotype with correct TAXA
write.csv(myY.115.3154, file = "data/myYflo.115.3154.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myQ$Taxa <- make.names(myQ$Taxa)
myY.115.3154$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.115.3154 <- myQ[match(myY.115.3154$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.115.3154, file = "data/myQflo.115.3154.csv",row.names = FALSE)


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
myGD.106.4178<- Select.MAF(SNPnmcol)

#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.106.4178))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.106.4178 <- getmyGD(myGM,myGD.106.4178)
myGM.106.4178 <- getmyGM(myGM,myGD.106.4178)
###check if they are shared the same name in the same order
identical(as.character(myGM.106.4178$Name), colnames(myGD.106.4178))
write.csv(myGM.106.4178, file = "data/myGMflo.106.4178.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.106.4178 <- FirstColumn(myGD.106.4178)
myGD.106.4178$Taxa <- make.names(myGD.106.4178$Taxa)
write.csv(myGD.106.4178, file = "data/myGDflo.106.4178.csv",row.names = FALSE)

##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.106.4178 <- allblupflo[match(myGD.106.4178$Taxa,allblupflo$Taxa, nomatch=0),]
myY.106.4178 <- droplevels(myY.106.4178)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.106.4178)
#### write out the phenotype with correct TAXA
write.csv(myY.106.4178, file = "data/myYflo.106.4178.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myQ$Taxa <- make.names(myQ$Taxa)
myY.106.4178$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.106.4178 <- myQ[match(myY.106.4178$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.106.4178, file = "data/myQflo.106.4178.csv",row.names = FALSE)


