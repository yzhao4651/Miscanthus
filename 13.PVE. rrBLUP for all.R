### this one for MLMSUPER
myYSRD <- read.csv("data/myYimputedSNP19SRD.csv")
myYOWA <- read.csv("data/myYimputedSNP19OWA3.csv")
myY1 <- read.csv("data/myYimputedSNP19.csv")
myY1 <- myY1[,-c(2,18)]
myY <- Reduce(function(x, y) merge(x, y, all.x=TRUE), list(myY1,myYOWA,myYSRD))
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
flo.1 <- read.csv("Result.7.25/rrBLUPalltrait.PCA.csv")
str(flo.1)
levels(flo.1$Ind.SNP)
#colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
levels(flo.1$Trait.name)
flo <- flo.1[flo.1$Ind.SNP=="143+36088im",]

flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(flo$Trait.name)
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.7.15.116im")
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.7.25.116im")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/rrBLUPPCA.36088im")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}


PVE.fun <- function(sample,myY,myGD){
  if (nrow(sample)==1){
    ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% sample$SNP])
    colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
    flowerSNPs <- merge(myY, ftdGDSNPs,by="Taxa")
    a <- ncol(myY)+1
    lm1 <-lm(flowerSNPs[[levels(sample$Trait.name)]] ~ flowerSNPs[[a]])
    af <- anova(lm1)
    afss <- af$"Sum Sq"
    PctExp <- afss/sum(afss)*100
    PVE3 <- data.frame(cbind(colnames(ftdGDSNPs[2]),PctExp[[1]]))
    colnames(PVE3) <- c("SNP","PVE")
  } else {
    ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% sample$SNP])
    colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
    flowerSNPs <- merge(myY, ftdGDSNPs,by="Taxa")
    number=1
    a <- ncol(myY)+1
    b <- a+nrow(sample)-1
    out_variable = colnames(flowerSNPs[a:b])
    outcome <- matrix(NA, nrow=ncol(flowerSNPs[a:b]),
                      ncol = 1, dimnames = list(out_variable, "PctExp"))
    for(j in a:b){
      lm1 <-lm(flowerSNPs[[levels(sample$Trait.name)]] ~ flowerSNPs[[j]])
      af <- anova(lm1)
      afss <- af$"Sum Sq"
      PctExp <- afss/sum(afss)*100
      outcome[number] <-PctExp[[1]]
      number=number+1
      PVE3 <- data.frame(outcome)
      colnames(PVE3) <-  c("PVE")
      SNP <- rownames(PVE3)
      rownames(PVE3) <- NULL
      PVE3 <- cbind(SNP,PVE3)
    }
  }
  return(PVE3)
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPPCA.36088im/*.", extension, sep = ""))

#myY <- read.csv("data/myYimputedSNP19.csv")
myGD <- read.csv("data/myGDimputedSNP19.csv")

mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
rrBLUPPCAim <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(rrBLUPPCAim$filename)

rrBLUPPCAim$filename <- gsub('rrBLUPPCA.36088im/Trait.', '', rrBLUPPCAim$filename)
rrBLUPPCAim$filename <- gsub('.csv', '', rrBLUPPCAim$filename)

levels(as.factor(rrBLUPPCAim$filename))

rrBLUPPCAim <- unite(rrBLUPPCAim,"SNP.T",SNP,filename,sep="/")
flo <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

rrBLUPPCAim.flo <- merge(flo,rrBLUPPCAim,by="SNP.T")
rrBLUPPCAim.flo.S <- separate(rrBLUPPCAim.flo, SNP.T, c("SNP","Trait.name"),sep="/")


###for flowering 106
###for flowering 106
###for flowering 106
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYf.106.4185.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDf.106.4185.csv")
flo.1 <- read.csv("Result.7.25/rrBLUPalltrait.PCA.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="F106+4185",]
flo <- droplevels(flo)
levels(flo.$Trait.name)
flo <- flo[!(flo$Trait.name=="Surv"),]
levels(flo$Trait.name)

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/rrBLUPPCA.F106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPPCA.F106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
rrBLUPPCA.F106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(rrBLUPPCA.F106$filename)

rrBLUPPCA.F106$filename <- gsub('rrBLUPPCA.F106/Trait.', '', rrBLUPPCA.F106$filename)
rrBLUPPCA.F106$filename <- gsub('.csv', '', rrBLUPPCA.F106$filename)

levels(as.factor(rrBLUPPCA.F106$filename))

rrBLUPPCA.F106 <- unite(rrBLUPPCA.F106,"SNP.T",SNP,filename,sep="/")
flo <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

rrBLUPPCA.F106.flo <- merge(flo,rrBLUPPCA.F106,by="SNP.T")
rrBLUPPCA.F106.flo.S <- separate(rrBLUPPCA.F106.flo, SNP.T, c("SNP","Trait.name"),sep="/")


###for Culm 106
###for Culm 106
###for Culm 106
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYC.106.4202.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDC.106.4202.csv")
flo.1 <- read.csv("Result.7.25/rrBLUPalltrait.PCA.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="C106+4202",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/rrBLUPPCA.C106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPPCA.C106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
rrBLUPPCA.C106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(rrBLUPPCA.C106$filename)

rrBLUPPCA.C106$filename <- gsub('rrBLUPPCA.C106/Trait.', '', rrBLUPPCA.C106$filename)
rrBLUPPCA.C106$filename <- gsub('.csv', '', rrBLUPPCA.C106$filename)

levels(as.factor(rrBLUPPCA.C106$filename))

rrBLUPPCA.C106 <- unite(rrBLUPPCA.C106,"SNP.T",SNP,filename,sep="/")
flo <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

rrBLUPPCA.C106.flo <- merge(flo,rrBLUPPCA.C106,by="SNP.T")
rrBLUPPCA.C106.flo.S <- separate(rrBLUPPCA.C106.flo, SNP.T, c("SNP","Trait.name"),sep="/")


###for Culm 116
###for Culm 116
###for Culm 116
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYC.116.3293.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDC.116.3293.csv")
flo.1 <- read.csv("Result.7.25/rrBLUPalltrait.PCA.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="C116+3293",]
str(flo)
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/rrBLUPPCA.C116")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}



setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPPCA.C116/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  str(sample)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
rrBLUPPCA.C116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(rrBLUPPCA.C116$filename)

rrBLUPPCA.C116$filename <- gsub('rrBLUPPCA.C116/Trait.', '', rrBLUPPCA.C116$filename)
rrBLUPPCA.C116$filename <- gsub('.csv', '', rrBLUPPCA.C116$filename)

levels(as.factor(rrBLUPPCA.C116$filename))

rrBLUPPCA.C116 <- unite(rrBLUPPCA.C116,"SNP.T",SNP,filename,sep="/")
flo <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

rrBLUPPCA.C116.flo <- merge(flo,rrBLUPPCA.C116,by="SNP.T")
rrBLUPPCA.C116.flo.S <- separate(rrBLUPPCA.C116.flo, SNP.T, c("SNP","Trait.name"),sep="/")

###for Culm 125
###for Culm 125
###for Culm 125
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYC.125.2562.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDC.125.2562.csv")
flo.1 <- read.csv("Result.7.25/rrBLUPalltrait.PCA.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="C125+2562",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/rrBLUPPCA.C125")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPPCA.C125/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
rrBLUPPCA.C125 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(rrBLUPPCA.C125$filename)

rrBLUPPCA.C125$filename <- gsub('rrBLUPPCA.C125/Trait.', '', rrBLUPPCA.C125$filename)
rrBLUPPCA.C125$filename <- gsub('.csv', '', rrBLUPPCA.C125$filename)

levels(as.factor(rrBLUPPCA.C125$filename))

rrBLUPPCA.C125 <- unite(rrBLUPPCA.C125,"SNP.T",SNP,filename,sep="/")
flo <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

rrBLUPPCA.C125.flo <- merge(flo,rrBLUPPCA.C125,by="SNP.T")
rrBLUPPCA.C125.flo.S <- separate(rrBLUPPCA.C125.flo, SNP.T, c("SNP","Trait.name"),sep="/")

###combine all of them together 
rrBLUPPCA.PVE.a <- do.call("rbind",list(rrBLUPPCA.C106.flo.S,rrBLUPPCA.C116.flo.S,rrBLUPPCA.C125.flo.S,rrBLUPPCAim.flo.S))
write.csv(rrBLUPPCA.PVE.a,file="Result.7.25/rrBLUPPCA.PVE.a.csv")



