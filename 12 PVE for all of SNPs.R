setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")

### calculating the PVE for each of SNPs 
##I am trying to use loop to get all of the results one time
extension <- "csv"
fileNames <- Sys.glob(paste("allresults/*_2.", extension, sep = ""))#with _2, because some other files without_1, _2,_3 are not SNPs results
###actually, there is one file called full_data, which including all of the SNPs for all traits. 
### Maybe it will be better to use full_data to do this loop? 
myY <- read.csv("data/myYimputedSNP19.csv")
myGD <- read.csv("data/myGDimputedSNP19.csv")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  for(j in 1:nlevels(sample$Method)){
    sample_M <- sample[which(sample$Method==j),]
    ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% sample_M$SNP])
    ###check if the new data set match with names from FarmCPU
    ####change the first column to the row name:
    ftdGDSNPs2 <- ftdGDSNPs[,-1]
    rownames(ftdGDSNPs2) <- ftdGDSNPs[,1]
    #rownames(kinshipftdGDSNPs7) <- kinshipftdGDSNPs6[,1]
    ##### Phenotype is ftdphenotype
    #####genotype is ftdGDSNPs
    colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
    flowerSNPs <- merge(myY, ftdGDSNPs,by="Taxa")
    #### method according to the Xiaolei Liu
    lm1 <-lm(flowerSNPs[[levels(sample$Trait.name)]] ~ flowerSNPs[,(ncol(myY)+1): ncol(flowerSNPs)], 
             flowerSNPs) ##i do not think i did right for indicating the multiple variables here. 
    ###i tried below one, but each method will give different number of SNPs, do you have any idea to 
    ###indicate each of SNPs from each method
    lm1 <-lm(flowerSNPs$ADD ~ colnames(ftdGDSNPs2)[1]+colnames(ftdGDSNPs2)[2]+colnames(ftdGDSNPs2)[3]+
               colnames(ftdGDSNPs2)[3], flowerSNPs)
    af <- anova(lm1)
    afss <- af$"Sum Sq"
    PctExp <- afss/sum(afss)*100
    mzList[[i]] <-  data.frame(PctExp)
  }
}
resultPVE = do.call("rbind", mzList)
