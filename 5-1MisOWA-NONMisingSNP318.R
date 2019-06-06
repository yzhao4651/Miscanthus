############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
allblup <- read.csv("data/alltraitsOWA.csv", na.strings = c("",".","NA"),row.names = 1)
allblup <- allblup[,c(1,3)]
str(allblup)
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
SNP <- subsetNONmisingSNP.GAPIT(datacomb3,allblup)

##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 38, 
### Remove the row with 38% missing value (more than 38% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.318),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.126.3659 <- Select.MAF(SNPnmcol)

str(myGD.126.3659)
row.names(myGD.126.3659)

#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
myGM$Name
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.126.3659))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.126.3659 <- getmyGD(myGM,myGD.126.3659)
str(myGD.126.3659)
myGM.126.3659 <- getmyGM(myGM,myGD.126.3659)
myGM.126.3659$Name
str(myGM.126.3659)
###check if they are shared the same name in the same order
identical(as.character(myGM.126.3659$Name), colnames(myGD.126.3659))
write.csv(myGM.126.3659, file = "data/myGMO.135.3628.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.126.3659 <- FirstColumn(myGD.126.3659)
str(myGD.126.3659)
write.csv(myGD.126.3659, file = "data/myGDO.126.3659.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.126.3659 <- allblup[match(myGD.126.3659$Taxa,allblup$Taxa, nomatch=0),]
myY.126.3659 <- droplevels(myY.126.3659)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.126.3659)
#### write out the phenotype with correct TAXA
write.csv(myY.126.3659, file = "data/myYO.126.3659.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.126.3659$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.126.3659 <- myQ[match(myY.126.3659$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.126.3659, file = "data/myQO.126.3659.csv",row.names = FALSE)








