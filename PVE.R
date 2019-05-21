setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
###import the data fprinW 
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinW.GWAS.Results.csv")
str(fprinW)
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
str(fprinW)
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
###import the total SNP with hit gene function 
allsnpblast <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/allsnpblast.csv")
str(allsnpblast)
colnames(allsnpblast)[colnames(allsnpblast)=="Original.marker.name"] <- "SNP"
fprinW$SNP %in% allsnpblast$SNP
##join both data sets
SNPall <- plyr::join_all(list(fprinW,allsnpblast), by="SNP")
str(SNPall)
SNPallpPC_w <- subset(SNPall,SNPall$P.value < 0.05)
write.csv(SNPallpPC_w , "Flowerresult/SNPallpPC_w.csv", row.names = FALSE)


###import the data fprind
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprind.GWAS.Results.csv")
str(fprinW)
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
str(fprinW)

###import the total SNP with hit gene function 
allsnpblast <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/allsnpblast.csv")
str(allsnpblast)
colnames(allsnpblast)[colnames(allsnpblast)=="Original.marker.name"] <- "SNP"
fprinW$SNP %in% allsnpblast$SNP
##join both data sets
SNPall <- plyr::join_all(list(fprinW,allsnpblast), by="SNP")
str(SNPall)
SNPallpPC_D <- subset(SNPall,SNPall$P.value < 0.05)
write.csv(SNPallpPC_D, "Flowerresult/SNPallpPC_D.csv", row.names = FALSE)

###import the data fprinGW
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinGW.GWAS.Results.csv")
str(fprinW)
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
str(fprinW)
###import the total SNP with hit gene function 
allsnpblast <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/allsnpblast.csv")
str(allsnpblast)
colnames(allsnpblast)[colnames(allsnpblast)=="Original.marker.name"] <- "SNP"
fprinW$SNP %in% allsnpblast$SNP
##join both data sets
SNPall <- plyr::join_all(list(fprinW,allsnpblast), by="SNP")
str(SNPall)
SNPallpPC_GW <- subset(SNPall,SNPall$P.value < 0.05)
write.csv(SNPallpPC_GW, "Flowerresult/SNPallpPC_GW.csv", row.names = FALSE)


###import the data fprinM
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinM.GWAS.Results.csv")
str(fprinW)
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
str(fprinW)
###import the total SNP with hit gene function 
allsnpblast <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/allsnpblast.csv")
str(allsnpblast)
colnames(allsnpblast)[colnames(allsnpblast)=="Original.marker.name"] <- "SNP"
fprinW$SNP %in% allsnpblast$SNP
##join both data sets
SNPall <- plyr::join_all(list(fprinW,allsnpblast), by="SNP")
str(SNPall)
SNPallpPC_M <- subset(SNPall,SNPall$P.value < 0.05)
write.csv(SNPallpPC_M, "Flowerresult/SNPallpPC_M.csv", row.names = FALSE)

###import the data FD_1
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FD_1.GWAS.Results.csv")
str(fprinW)
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
str(fprinW)
###import the total SNP with hit gene function 
allsnpblast <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/allsnpblast.csv")
str(allsnpblast)
colnames(allsnpblast)[colnames(allsnpblast)=="Original.marker.name"] <- "SNP"
fprinW$SNP %in% allsnpblast$SNP
##join both data sets
SNPall <- plyr::join_all(list(fprinW,allsnpblast), by="SNP")
str(SNPall)
SNPallpFD_1 <- subset(SNPall,SNPall$P.value < 0.05)
write.csv(SNPallpFD_1 , "Flowerresult/SNPallpFD_1 .csv", row.names = FALSE)


###import the data FW_1
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FW_1.GWAS.Results.csv")
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
str(fprinW)
###import the total SNP with hit gene function 
allsnpblast <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/allsnpblast.csv")
str(allsnpblast)
colnames(allsnpblast)[colnames(allsnpblast)=="Original.marker.name"] <- "SNP"
fprinW$SNP %in% allsnpblast$SNP
##join both data sets
SNPall <- plyr::join_all(list(fprinW,allsnpblast), by="SNP")
str(SNPall)
SNPallpFW_1 <- subset(SNPall,SNPall$P.value < 0.05)
write.csv(SNPallpFW_1, "Flowerresult/SNPallpFW_1.csv", row.names = FALSE)



###import the data GFW_1
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.GFW_1.GWAS.Results.csv")
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
str(fprinW)
###import the total SNP with hit gene function 
allsnpblast <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/allsnpblast.csv")
str(allsnpblast)
colnames(allsnpblast)[colnames(allsnpblast)=="Original.marker.name"] <- "SNP"
fprinW$SNP %in% allsnpblast$SNP
##join both data sets
SNPall <- plyr::join_all(list(fprinW,allsnpblast), by="SNP")
str(SNPall)
SNPallpGFW_1 <- subset(SNPall,SNPall$P.value < 0.05)
write.csv(SNPallpGFW_1, "Flowerresult/SNPallpGFW_1.csv", row.names = FALSE)


###import the data FM_1
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FM_1.GWAS.Results.csv")
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
str(fprinW)
###import the total SNP with hit gene function 
allsnpblast <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/allsnpblast.csv")
str(allsnpblast)
colnames(allsnpblast)[colnames(allsnpblast)=="Original.marker.name"] <- "SNP"
fprinW$SNP %in% allsnpblast$SNP
##join both data sets
SNPall <- plyr::join_all(list(fprinW,allsnpblast), by="SNP")
str(SNPall)
SNPallpFM_1 <- subset(SNPall,SNPall$P.value < 0.05)
write.csv(SNPallpFM_1, "Flowerresult/SNPallpFM_1.csv", row.names = FALSE)

setwd("Flowerresult/")
extension <- "csv"
fileNames <- Sys.glob(paste("SNPallp*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mzList[[i]] = data.frame(sample, filename = rep(fileNames[i]))
}
resultFarmCPU1 = do.call("rbind", mzList)
str(resultFarmCPU)
resultFarmCPU1 <- resultFarmCPU1[,c(20,1:5,7:19)]
write.csv(resultFarmCPU1,"Flowerresult/resultFarmCPU1.csv",row.names = FALSE)


####fprinW
##get the the first colume: the name of the SNPs
### row change to colume 
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinW.GWAS.Results.csv")
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
fprinW <- subset(fprinW,fprinW$P.value < 0.05)
fprinW <- data.frame(t(fprinW))
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
##import the phenotype data set
#myGD <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myGDimputedSNP19.csv")
myY <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myYimputedSNP19.csv")
### change the first row to the column name 
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
### select the certain colume according to the fprinW
ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% names(fprinW)])
###check if the new data set match with names from FarmCPU
####change the first column to the row name:
ftdGDSNPs2 <- ftdGDSNPs[,-1]
rownames(ftdGDSNPs2) <- ftdGDSNPs[,1]
#rownames(kinshipftdGDSNPs7) <- kinshipftdGDSNPs6[,1]
##### Phenotype is ftdphenotype
#####genotype is ftdGDSNPs
colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
flowerSNPs=merge(myY, ftdGDSNPs,by="Taxa")
str(ftdGDSNPs)
#### method according to the Xiaolei Liu
lm1 <-lm(fprinW ~ PstI.TP202284+PstI.TP1027360+NsiI.TP632002+NsiI.TP728788, flowerSNPs)
af <- anova(lm1)
afss <- af$"Sum Sq"
PVEPC_W <- data.frame(cbind(row.names(af), PctExp=afss/sum(afss)*100))
write.csv(PVEPC_W,"Flowerresult/PVEPC_W.csv", row.names = FALSE)


##get the the first colume: the name of the SNPs
### row change to colume 
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprind.GWAS.Results.csv")
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
fprinW <- subset(fprinW,fprinW$P.value < 0.05)
fprinW <- data.frame(t(fprinW))
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
##import the phenotype data set
#myGD <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myGDimputedSNP19.csv")
myY <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myYimputedSNP19.csv")
### change the first row to the column name 
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
### select the certain colume according to the fprinW
ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% names(fprinW)])
###check if the new data set match with names from FarmCPU
####change the first column to the row name:
ftdGDSNPs2 <- ftdGDSNPs[,-1]
rownames(ftdGDSNPs2) <- ftdGDSNPs[,1]
#rownames(kinshipftdGDSNPs7) <- kinshipftdGDSNPs6[,1]
##### Phenotype is ftdphenotype
#####genotype is ftdGDSNPs
colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
flowerSNPs=merge(myY, ftdGDSNPs,by="Taxa")
str(ftdGDSNPs)
#### method according to the Xiaolei Liu
lm1 <-lm(fprind ~ NsiI.TP632002+NsiI.TP728788, flowerSNPs)
af <- anova(lm1)
afss <- af$"Sum Sq"
PVEPC_D <- data.frame(cbind(row.names(af),PctExp=afss/sum(afss)*100))
write.csv(PVEPC_D,"Flowerresult/PVEPC_D.csv", row.names = FALSE)

##get the the first colume: the name of the SNPs
### row change to colume 
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinGW.GWAS.Results.csv")
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
fprinW <- subset(fprinW,fprinW$P.value < 0.05)
fprinW <- data.frame(t(fprinW))
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
##import the phenotype data set
#myGD <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myGDimputedSNP19.csv")
myY <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myYimputedSNP19.csv")
### change the first row to the column name 
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
### select the certain colume according to the fprinW
ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% names(fprinW)])
###check if the new data set match with names from FarmCPU
####change the first column to the row name:
ftdGDSNPs2 <- ftdGDSNPs[,-1]
rownames(ftdGDSNPs2) <- ftdGDSNPs[,1]
#rownames(kinshipftdGDSNPs7) <- kinshipftdGDSNPs6[,1]
##### Phenotype is ftdphenotype
#####genotype is ftdGDSNPs
colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
flowerSNPs=merge(myY, ftdGDSNPs,by="Taxa")
str(ftdGDSNPs)
#### method according to the Xiaolei Liu
#### method according to the Xiaolei Liu
lm1 <-lm(fprinGW ~ PstI.TP202284+NsiI.TP42348+NsiI.TP222743+NsiI.TP632002+NsiI.TP728788, flowerSNPs)
af <- anova(lm1)
afss <- af$"Sum Sq"
PVEPC_GW <- data.frame(cbind(row.names(af),PctExp=afss/sum(afss)*100))
write.csv(PVEPC_GW,"Flowerresult/PVEPC_GW.csv", row.names = FALSE)

##get the the first colume: the name of the SNPs
### row change to colume 
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinM.GWAS.Results.csv")
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
fprinW <- subset(fprinW,fprinW$P.value < 0.05)
fprinW <- data.frame(t(fprinW))
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
##import the phenotype data set
#myGD <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myGDimputedSNP19.csv")
myY <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myYimputedSNP19.csv")
### change the first row to the column name 
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
### select the certain colume according to the fprinW
ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% names(fprinW)])
###check if the new data set match with names from FarmCPU
####change the first column to the row name:
ftdGDSNPs2 <- ftdGDSNPs[,-1]
rownames(ftdGDSNPs2) <- ftdGDSNPs[,1]
#rownames(kinshipftdGDSNPs7) <- kinshipftdGDSNPs6[,1]
##### Phenotype is ftdphenotype
#####genotype is ftdGDSNPs
colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
flowerSNPs=merge(myY, ftdGDSNPs,by="Taxa")
str(ftdGDSNPs)
#### method according to the Xiaolei Liu
#### method according to the Xiaolei Liu
lm1 <-lm(fprinM ~ PstI.TP114049+PstI.TP477722+NsiI.TP178513, flowerSNPs)
af <- anova(lm1)
afss <- af$"Sum Sq"
PVEPC_M <- data.frame(cbind(row.names(af),PctExp=afss/sum(afss)*100))
write.csv(PVEPC_M,"Flowerresult/PVEPC_M.csv", row.names = FALSE)

###
###
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FD_1.GWAS.Results.csv")
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
fprinW <- subset(fprinW,fprinW$P.value < 0.05)
fprinW <- data.frame(t(fprinW))
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
##import the phenotype data set
#myGD <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myGDimputedSNP19.csv")
myY <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myYimputedSNP19.csv")
### change the first row to the column name 
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
### select the certain colume according to the fprinW
ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% names(fprinW)])
###check if the new data set match with names from FarmCPU
####change the first column to the row name:
ftdGDSNPs2 <- ftdGDSNPs[,-1]
rownames(ftdGDSNPs2) <- ftdGDSNPs[,1]
#rownames(kinshipftdGDSNPs7) <- kinshipftdGDSNPs6[,1]
##### Phenotype is ftdphenotype
#####genotype is ftdGDSNPs
colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
flowerSNPs=merge(myY, ftdGDSNPs,by="Taxa")
str(ftdGDSNPs)
#### method according to the Xiaolei Liu
#### method according to the Xiaolei Liu
lm1 <-lm(FD_1 ~ PstI.TP114049+PstI.TP189070+PstI.TP477722+NsiI.TP321836+NsiI.TP789014, flowerSNPs)
af <- anova(lm1)
afss <- af$"Sum Sq"
PVEFD_1<- data.frame(cbind(row.names(af),PctExp=afss/sum(afss)*100))
write.csv(PVEFD_1,"Flowerresult/PVEFD_1.csv", row.names = FALSE)

fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FW_1.GWAS.Results.csv")
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
fprinW <- subset(fprinW,fprinW$P.value < 0.05)
fprinW <- data.frame(t(fprinW))
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
##import the phenotype data set
#myGD <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myGDimputedSNP19.csv")
myY <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myYimputedSNP19.csv")
### change the first row to the column name 
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
### select the certain colume according to the fprinW
ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% names(fprinW)])
###check if the new data set match with names from FarmCPU
####change the first column to the row name:
ftdGDSNPs2 <- ftdGDSNPs[,-1]
rownames(ftdGDSNPs2) <- ftdGDSNPs[,1]
#rownames(kinshipftdGDSNPs7) <- kinshipftdGDSNPs6[,1]
##### Phenotype is ftdphenotype
#####genotype is ftdGDSNPs
colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
flowerSNPs=merge(myY, ftdGDSNPs,by="Taxa")
str(ftdGDSNPs)
#### method according to the Xiaolei Liu
#### method according to the Xiaolei Liu
lm1 <-lm(FW_1 ~ PstI.TP114049+PstI.TP189070+PstI.TP477722+NsiI.TP321836+NsiI.TP789014, flowerSNPs)
af <- anova(lm1)
afss <- af$"Sum Sq"
PVEFW_1 <- data.frame(cbind(row.names(af),PctExp=afss/sum(afss)*100))
str(PVEFW_1)
write.csv(PVEFW_1,"Flowerresult/PVEFW_1.csv", row.names = FALSE)

fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.GFW_1.GWAS.Results.csv")
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
fprinW <- subset(fprinW,fprinW$P.value < 0.05)
fprinW <- data.frame(t(fprinW))
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
##import the phenotype data set
#myGD <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myGDimputedSNP19.csv")
myY <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myYimputedSNP19.csv")
### change the first row to the column name 
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
### select the certain colume according to the fprinW
ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% names(fprinW)])
###check if the new data set match with names from FarmCPU
####change the first column to the row name:
ftdGDSNPs2 <- ftdGDSNPs[,-1]
rownames(ftdGDSNPs2) <- ftdGDSNPs[,1]
#rownames(kinshipftdGDSNPs7) <- kinshipftdGDSNPs6[,1]
##### Phenotype is ftdphenotype
#####genotype is ftdGDSNPs
colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
flowerSNPs=merge(myY, ftdGDSNPs,by="Taxa")
str(ftdGDSNPs)
#### method according to the Xiaolei Liu
#### method according to the Xiaolei Liu
lm1 <-lm(GFW_1 ~ PstI.TP819930+NsiI.TP632002+NsiI.TP728788, flowerSNPs)
af <- anova(lm1)
afss <- af$"Sum Sq"
PVEGFW_1 <- data.frame(cbind(row.names(af),PctExp=afss/sum(afss)*100))
write.csv(PVEGFW_1,"Flowerresult/PVEGFW_1.csv", row.names = FALSE)

###FM_1
fprinW <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FM_1.GWAS.Results.csv")
fprinW <- fprinW[,c(1:4,7)]
fprinW <- fprinW[order(fprinW$P.value),]
fprinW$P.value <- p.adjust(fprinW$P.value, method ="bonferroni")
fprinW <- subset(fprinW,fprinW$P.value < 0.05)
fprinW <- data.frame(t(fprinW))
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
##import the phenotype data set
#myGD <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myGDimputedSNP19.csv")
myY <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/data/myYimputedSNP19.csv")
### change the first row to the column name 
colnames(fprinW) = as.character(unlist(fprinW[1, ])) # the first row will be the header
### select the certain colume according to the fprinW
ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% names(fprinW)])
###check if the new data set match with names from FarmCPU
####change the first column to the row name:
ftdGDSNPs2 <- ftdGDSNPs[,-1]
rownames(ftdGDSNPs2) <- ftdGDSNPs[,1]
#rownames(kinshipftdGDSNPs7) <- kinshipftdGDSNPs6[,1]
##### Phenotype is ftdphenotype
#####genotype is ftdGDSNPs
colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
flowerSNPs=merge(myY, ftdGDSNPs,by="Taxa")
str(ftdGDSNPs)
#### method according to the Xiaolei Liu
#### method according to the Xiaolei Liu
lm1 <-lm(FM_1 ~ PstI.TP319062+NsiI.TP178513+NsiI.TP579718, flowerSNPs)
af <- anova(lm1)
afss <- af$"Sum Sq"
PVEFM_1 <- data.frame(cbind(row.names(af),PctExp=afss/sum(afss)*100))
write.csv(PVEFM_1,"Flowerresult/PVEFM_1.csv", row.names = FALSE)


setwd("Flowerresult/")
extension <- "csv"
fileNames <- Sys.glob(paste("PVE*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mzList[[i]] = data.frame(sample, filename = rep(fileNames[i]))
}
FarmCPUPVE = do.call("rbind", mzList)
str(FarmCPUPVE)
#FarmCPUPVE <- FarmCPUPVE[,c(20,1:5,7:19)]
write.csv(FarmCPUPVE,"FarmCPUPVE1.csv",row.names = FALSE)