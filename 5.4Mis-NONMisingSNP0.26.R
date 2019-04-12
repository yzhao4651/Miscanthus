############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
allblup <- read.csv("data/myYimputedSNP19.csv",na.strings = c("",".","NA"))

############# step2 
############# step2
############download the SNP data that Lindsay gave to me ################### 
############loading the genotype SNP data with missing value called datacomb3####
load("~/Documents/MiacnathusSNPinformation /160322filteredSNPs.RData")
########## inspect the genotype data
str(datacomb3) # 594 genotypes, 46,177 markers
dimnames(datacomb3)[[1]] # this is how the genotypes are named in the SNP data
#############changing the rowname to the first columne###########################
source("Function/subsetNONmissingSNP.GAPIT.R")
SNP <- subsetNONmisingSNP.GAPIT(datacomb3,allblup)

##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 26, 
### Remove the row with 26% missing value (more than 26% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.26),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.112.4877 <- Select.MAF(SNPnmcol)


#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)
### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.112.4877)[-1])
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.112.4877 <- getmyGD(myGM,myGD.112.4877)
myGM.112.4877 <- getmyGM(myGM,myGD.112.4877)
###check if they are shared the same name in the same order
identical(as.character(myGM.112.4877$Name), colnames(myGD.112.4877)[-1])
write.csv(myGM.112.4877, file = "data/myGM.112.4877.csv",row.names = FALSE)
write.csv(myGD.112.4877, file = "data/myGD.112.4877.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.112.4877 <- allblup[match(myGD.112.4877$Taxa,allblup$Taxa, nomatch=0),]
myY.112.4877 <- droplevels(myY.112.4877)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])
str(myY.112.4877)
#### write out the phenotype with correct TAXA
write.csv(myY.112.4877, file = "data/myY.112.4877.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.112.4877$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.112.4877 <- myQ[match(myY.112.4877$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.112.4877, file = "data/myQ.112.4877.csv",row.names = FALSE)








