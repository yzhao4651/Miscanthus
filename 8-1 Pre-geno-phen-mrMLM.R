############# this is prepare the dataset for mrMLMM software
############# this is prepare the dataset for mrMLMM software
############step 1: import all SNP genotype with imputed value
############step 1: import all SNP genotype with imputed value
###trying to load in to GitHub, But this one is big, can not load from here. 
load("~/Documents/MiacnathusSNPinformation /160324EMimputedSNP_Msi.RData")
#load("C:/Users/Admin/Desktop/New folder/miscanthus study-1/Misthcanthus CCA data analysis/160324EMimputedSNP_Msi.RData")
# load("~/DOE Msi study/yield manuscript/phenotypic and GWAS analysis/160324EMimputedSNP_Msi.RData") # Lindsay's version
### this data set in value part, change the value into data.frame
##check the name of that value:
names(myA.EM.Msi)
myA.EM.Msi$imputed[myA.EM.Msi$imputed > 1] <- 1
myA.EM.Msi$imputed[myA.EM.Msi$imputed < -1] <- -1
####change the data to data.frame
datacomb2 <- data.frame(myA.EM.Msi[["imputed"]]) 
###change the column to the row, which make the row is SNPs and column is individuals 
datacomb2tran <- data.frame(t(datacomb2))
### change row name to the first column 
rn <- rownames(datacomb2tran)
rownames(datacomb2tran) <- NULL
datacomb2tran <- cbind(rn,datacomb2tran)
###check the data
str(datacomb2tran)

#### step2: import the allSNP including genotype.for.code.1 to merge with the all SNP from above step  
#### step2: import the allSNP including genotype.for.code.1 to merge with the all SNP from above step 
###import the all marker information ( this one already get the alleles)
#install.packages("readxl")
library(readxl)
allsnp2 <- data.frame(read_excel("data/allsnp.xlsx"))# I got all SNPs from the webs you sent to me and manually 
                                                     ##changed some parts and then get "genotype for code 1"
###change the data all snp to data frame
colnames(allsnp2)[colnames(allsnp2)=="Original.marker.name"] <- "rn"
###check if the name in the datacomb3 is same to the name in the allsnps
all(allsnp2[,2] %in% as.character(datacomb2tran[,1]))
#str(allsnp)
####find the datacomb2 data match the allsnpsub, pay attention the datacomb2tran should put the front of the %in%.
datacomb2transub <- subset(datacomb2tran, (as.character(datacomb2tran[,1]) %in% allsnp2[,2]))
###select several column from all of the allsnps
allsnp2sub <- data.frame(allsnp2[,c(1:2,3:7)])
str(allsnp2sub)

#### merge allsnpsub and datacomb2transub
datacomb2tranmerge <- merge(allsnp2sub,datacomb2transub,by="rn")
##checking the data
str(datacomb2tranmerge)
### get the genotype with the genotype.for. code.1 for using later on 
datacomb2tranmerge2 <- datacomb2tranmerge[,c(1,5,8:582)]
###check the data again
str(datacomb2tranmerge2)

###step 3 load GM data sets to select matched Genotype from step 2
###step 3 load GM data sets to select matched Genotype from step 2
###importing the myGM from data set
myGM <- read.csv("data/myGM.csv",row.names=1)
str(myGM)
all(myGM$Name %in% datacomb2tranmerge2$rn) # This comes out FALSE
###get all the match with the datacomb2tranmerge
myGM2 <- myGM[match(datacomb2tranmerge2$rn, myGM$Name, nomatch=0),] 
###change the column name of the myGM3 to rn#
colnames(myGM2)[which(names(myGM2) == "Name")] <- "rn"
###merger myGM and datacomb2tranmerge2, after they have the same name
datacomb2tranmerge4 <- merge(myGM2, datacomb2tranmerge2, by="rn")
snporder <- order(datacomb2tranmerge4$Chromosome, datacomb2tranmerge4$Position)
datacomb2tranmerge4 <- datacomb2tranmerge4[snporder,]
###check 
str(datacomb2tranmerge3)

### write out the datasets with 
allsnp2sub <- data.frame(allsnp2[,c(2,5)])
str(allsnp2sub)
write.csv(allsnp2sub, file = "data/allsnp2sub.csv", row.names = FALSE, na = "NA")
allsnp2sub <- read.csv("data/allsnp2sub.csv")

###step 4 import the myGM dataset for the GAPIT to select matched Genotype from step 3
###step 4 import the myGM dataset for the GAPIT to select matched Genotype from step 3
####import the GM data set for GAPIT to select the genotype dataset
myGMmrMLMM <- read.csv("data/myGMimputedSNP19.csv")
str(myGMmrMLMM)
###check they both if get the same name
all(myGMmrMLMM$Name %in% datacomb2tranmerge4$rn) # comes out to FALSE
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno<- datacomb2tranmerge4[match(myGMmrMLMM$Name, datacomb2tranmerge4$rn,nomatch=0),]
str(subgeno)

###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
###import the genotype and phonetype:

myYmrMlMM<- read.csv("data/myYimputedSNP19.csv")
#str(myYmrMlMM)
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
str(myYmrMlMM)
###subset the phonetype which match from myY
### change the first column to the row name in order to do transform 
subgeno1 <- data.frame(subgeno,row.names = 1)
str(subgeno1)
### transform the data: and then the row is individual and column is SNPs
subgenotran <- data.frame(t(subgeno1[4:577]))
str(subgenotran)
###check if the name of phenotype is the same to the name of genotype
all(myYmrMlMM$Taxa %in% dimnames(subgenotran)[[1]]) # TRUE
### select the the individuals matched the individual in phenotype dataset 
subgenotranonly<- subgenotran[match(myYmrMlMM$Taxa,dimnames(subgenotran)[[1]],nomatch=0),]
subgenotranonly <- data.frame(subgenotranonly)+1
source("Function/SelectMAF-mrMLMM.R")
subgenotranonly <- Select.MAF(subgenotranonly)

###combine all (emger part of subgeno and subgenotranonly)
##subgenotranonly$rn <- gsub("\\.", "-",subgenotranonly$rn)
#str(subgeno[,1:4])
#str(subgenotranonly)
subgenomrMLMM <- plyr::join(data.frame(subgeno[,1:4]),data.frame(subgenotranonly),by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "data/subgenomrMLMM.csv", row.names = FALSE, na = "NA")
#write.csv(subgenomrMLMM, file = "C:/Users/Admin/Desktop/Miscanthus/Miscanthus/subgenomrMLMM.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###step6 chagne the name of the phenotype and also seperate them with less missing values

###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myYmrMlMM<- read.csv("data/myYimputedSNP19.csv")
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
myYmrMlMM$Taxa 
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"

##depending on the missing value and then separate the data to three data set: 

myYmrMlMMSuv <- myYmrMlMM[1:2]
myYmrMlMMflo <- myYmrMlMM[,c(1,14:16,20:39)]
myYmrMlMMculm <- myYmrMlMM[,c(1,3:13,17:19)]

###get the column names of phenotype in order to check the later one 
#names(myYmrMlMM)
##[1] "Taxa"      "Surv"      "CmDW_g"    "Cml_cm"    "CmD_BI_mm" "CmD_LI_mm" "CmN."      "Bcirc_cm"  "Yld_kg"    "SDW_kg"    "CCirc_cm" 
##[12] "Lg"        "GS"        "TFN."      "FNMain"    "FNsmall"   "FD"        "SRD"       "ADD"       "HD_1"      "FD_1"      "HD_50."   
##[23] "FD_50."    "HW_1"      "FW_1"      "HW_50"     "FW_50"     "GHW_1"     "GFW_1"     "GHW_50"    "GFW_50"    "HM_1"      "FM_1"     
##[34] "HM_50"     "FM_50"     "fprind"    "fprinW"    "fprinGW"   "fprinM"
#names(myYmrMlMMculm)
#[1] "Taxa"      "CmDW_g"    "Cml_cm"    "CmD_BI_mm" "CmD_LI_mm" "CmN."      "Bcirc_cm"  "Yld_kg"    "SDW_kg"    "CCirc_cm"  "Lg"       
#[12] "GS"        "FD"        "SRD"       "ADD"  
#names(myYmrMlMMflo)
#[1] "Taxa"    "TFN."    "FNMain"  "FNsmall" "HD_1"    "FD_1"    "HD_50."  "FD_50."  "HW_1"    "FW_1"    "HW_50"   "FW_50"   "GHW_1"   "GFW_1"  
#[15] "GHW_50"  "GFW_50"  "HM_1"    "FM_1"    "HM_50"   "FM_50"   "fprind"  "fprinW"  "fprinGW" "fprinM" 
###using this function to change the name of traits names
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
###check the column name if changed 


names(myYmrMlMM)
###write out the dataset
write.csv(myYmrMlMM, file = "mrMLMM2/myYmrMlMM.csv", row.names = FALSE, na = "NA")
myYmrMlMMSuv <- myYmrMlMM[1:2]
write.csv(myYmrMlMMSuv, file = "mrMLMM2/myYmrMlMMSuv.csv", row.names = FALSE, na = "NA")
###Flowering traits 
myYmrMlMMflo <- myYmrMlMM[,c(1,14:16,20:39)]
###omit the rews with all missing values 
myYmrMlMMflo<- myYmrMlMMflo[-which(rowSums(is.na(myYmrMlMMflo)) ==23),]
###write out the dataset
write.csv(myYmrMlMMflo, file = "mrMLMM2/myYmrMlMMflo3.csv", row.names = FALSE, na = "NA")

###culm traits
myYmrMlMMculm <- myYmrMlMM[,c(1,3:13,17:19)]
myYmrMlMMculm<- myYmrMlMMculm[-which(rowSums(is.na(myYmrMlMMculm)) ==14),]
###write out the dataset
write.csv(myYmrMlMMculm, file = "mrMLMM2/myYmrMlMMculm3.csv", row.names = FALSE, na = "NA")
