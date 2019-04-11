############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
allblup <- read.csv("data/alltraits.csv", na.strings = c("",".","NA"),row.names = 1)
str(allblup)
allblupflo <- allblup[,c(1, 3:13,17:19)]
allblupflo<- allblupflo[-which(rowSums(is.na(allblupflo[1:15])) == 14),]
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
########1)remove row with missing value more than 38, 
### Remove the row with 38% missing value (more than 38% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.55),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.127.2216 <- Select.MAF(SNPnmcol)

#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
myGM$Name
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.127.2216))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.127.2216 <- getmyGD(myGM,myGD.127.2216)
str(myGD.127.2216)
myGM.127.2216 <- getmyGM(myGM,myGD.127.2216)
myGM.127.2216$Name
str(myGM.127.2216)
###check if they are shared the same name in the same order
identical(as.character(myGM.127.2216$Name), colnames(myGD.127.2216))
write.csv(myGM.127.2216, file = "data/myGMf.127.2216.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.127.2216 <- FirstColumn(myGD.127.2216)
str(myGD.127.2216)
write.csv(myGD.127.2216, file = "data/myGDf.127.2216.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.127.2216 <- allblup[match(myGD.127.2216$Taxa,allblup$Taxa, nomatch=0),]
myY.127.2216 <- droplevels(myY.127.2216)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.127.2216)
#### write out the phenotype with correct TAXA
write.csv(myY.127.2216, file = "data/myYf.127.2216.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.127.2216$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.127.2216 <- myQ[match(myY.127.2216$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.127.2216, file = "data/myQf.127.2216.csv",row.names = FALSE)








