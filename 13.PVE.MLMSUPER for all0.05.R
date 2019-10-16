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
flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
str(flo.1)

flo <- flo.1[flo.1$Ind.SNP=="+36088im" & flo.1$Method=="MLM+SUPER",]
str(flo)
names(flo)
levels()
flo <- droplevels(flo)
levels(flo$Ind.SNP)
levels(as.factor(flo$Method))
levels(as.factor(flo$Trait.name))
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="fprin"),]
levels(as.factor(flo$Trait.name))
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.7.15.116im")
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.8.3.116im")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.36088im")
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
fileNames <- Sys.glob(paste("MLMSUPER0.05.36088im/*.", extension, sep = ""))

#myY <- read.csv("data/myYimputedSNP19.csv")
myGD <- read.csv("data/myGDimputedSNP19.csv")

mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMCUPERim <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMCUPERim$filename)

MLMCUPERim$filename <- gsub('MLMSUPER0.05.36088im/Trait.', '', MLMCUPERim$filename)
MLMCUPERim$filename <- gsub('.csv', '', MLMCUPERim$filename)

levels(as.factor(MLMCUPERim$filename))

MLMCUPERim.1 <- unite(MLMCUPERim,"SNP.T",SNP,filename,sep="/")
flo.2<- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMCUPERim.flo <- merge(flo.2,MLMCUPERim.1,by="SNP.T")
MLMCUPERim.flo.S <- separate(MLMCUPERim.flo, SNP.T, c("SNP","Trait.name"),sep="/")
levels(as.factor(MLMCUPERim.flo.S$Trait.name))


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

flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
str(flo.1)

flo <- flo.1[flo.1$Ind.SNP=="F106+4185" & flo.1$Method=="MLM+SUPER",]

flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.F106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.F106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMSUPER0.05.F106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMSUPER0.05.F106$filename)

MLMSUPER0.05.F106$filename <- gsub('MLMSUPER0.05.F106/Trait.', '', MLMSUPER0.05.F106$filename)
MLMSUPER0.05.F106$filename <- gsub('.csv', '', MLMSUPER0.05.F106$filename)

levels(as.factor(MLMSUPER0.05.F106$filename))

MLMSUPER0.05.F106.1 <- unite(MLMSUPER0.05.F106,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMSUPER0.05.F106.flo <- merge(flo.2,MLMSUPER0.05.F106.1,by="SNP.T")
MLMSUPER0.05.F106.flo.S <- separate(MLMSUPER0.05.F106.flo, SNP.T, c("SNP","Trait.name"),sep="/")

###for flowering 116
###for flowering 116
###for flowering 116
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYf.116.3077.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDf.116.3077.csv")

flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
str(flo.1)

flo <- flo.1[flo.1$Ind.SNP=="F116+3077" & flo.1$Method=="MLM+SUPER",]

flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.F116")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.F116/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMSUPER0.05.F116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMSUPER0.05.F116$filename)

MLMSUPER0.05.F116$filename <- gsub('MLMSUPER0.05.F116/Trait.', '', MLMSUPER0.05.F116$filename)
MLMSUPER0.05.F116$filename <- gsub('.csv', '', MLMSUPER0.05.F116$filename)

levels(as.factor(MLMSUPER0.05.F116$filename))

MLMSUPER0.05.F116.1 <- unite(MLMSUPER0.05.F116,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMSUPER0.05.F116.flo <- merge(flo.2,MLMSUPER0.05.F116.1,by="SNP.T")
MLMSUPER0.05.F116.flo.S <- separate(MLMSUPER0.05.F116.flo, SNP.T, c("SNP","Trait.name"),sep="/")




###for flowering 122
###for flowering 122
###for flowering 122
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYf.122.2272.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDf.122.2272.csv")
flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
str(flo.1)
flo <- flo.1[flo.1$Ind.SNP=="F122+2272" & flo.1$Method=="MLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.F122")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.F122/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMSUPER0.05.F122 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMSUPER0.05.F122$filename)

MLMSUPER0.05.F122$filename <- gsub('MLMSUPER0.05.F122/Trait.', '', MLMSUPER0.05.F122$filename)
MLMSUPER0.05.F122$filename <- gsub('.csv', '', MLMSUPER0.05.F122$filename)

levels(as.factor(MLMSUPER0.05.F122$filename))

MLMSUPER0.05.F122.1 <- unite(MLMSUPER0.05.F122,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMSUPER0.05.F122.flo <- merge(flo.2,MLMSUPER0.05.F122.1,by="SNP.T")
MLMSUPER0.05.F122.flo.S <- separate(MLMSUPER0.05.F122.flo, SNP.T, c("SNP","Trait.name"),sep="/")


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
flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="C106+4202" & flo.1$Method=="MLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.C106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.C106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMSUPER0.05.C106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMSUPER0.05.C106$filename)

MLMSUPER0.05.C106$filename <- gsub('MLMSUPER0.05.C106/Trait.', '', MLMSUPER0.05.C106$filename)
MLMSUPER0.05.C106$filename <- gsub('.csv', '', MLMSUPER0.05.C106$filename)

levels(as.factor(MLMSUPER0.05.C106$filename))

MLMSUPER0.05.C106.1 <- unite(MLMSUPER0.05.C106,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMSUPER0.05.C106.flo <- merge(flo.2,MLMSUPER0.05.C106.1,by="SNP.T")
MLMSUPER0.05.C106.flo.S <- separate(MLMSUPER0.05.C106.flo, SNP.T, c("SNP","Trait.name"),sep="/")


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
flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="C116+3293" & flo.1$Method=="MLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.C116")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.C116/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMSUPER0.05.C116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMSUPER0.05.C116$filename)

MLMSUPER0.05.C116$filename <- gsub('MLMSUPER0.05.C116/Trait.', '', MLMSUPER0.05.C116$filename)
MLMSUPER0.05.C116$filename <- gsub('.csv', '', MLMSUPER0.05.C116$filename)

levels(as.factor(MLMSUPER0.05.C116$filename))

MLMSUPER0.05.C116.1 <- unite(MLMSUPER0.05.C116,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMSUPER0.05.C116.flo <- merge(flo.2,MLMSUPER0.05.C116.1,by="SNP.T")
MLMSUPER0.05.C116.flo.S <- separate(MLMSUPER0.05.C116.flo, SNP.T, c("SNP","Trait.name"),sep="/")

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
flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="C125+2562" & flo.1$Method=="MLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.C125")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.C125/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMSUPER0.05.C125 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMSUPER0.05.C125$filename)

MLMSUPER0.05.C125$filename <- gsub('MLMSUPER0.05.C125/Trait.', '', MLMSUPER0.05.C125$filename)
MLMSUPER0.05.C125$filename <- gsub('.csv', '', MLMSUPER0.05.C125$filename)

levels(as.factor(MLMSUPER0.05.C125$filename))

MLMSUPER0.05.C125.1 <- unite(MLMSUPER0.05.C125,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMSUPER0.05.C125.flo <- merge(flo.2,MLMSUPER0.05.C125.1,by="SNP.T")
MLMSUPER0.05.C125.flo.S <- separate(MLMSUPER0.05.C125.flo, SNP.T, c("SNP","Trait.name"),sep="/")


###for OWA 106
###for OWA 106
###for OWA 106
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYO.106.4834.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDO.106.4834.csv")
flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="O106+4834" & flo.1$Method=="MLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.O106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.O106/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMSUPER0.05.O106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMSUPER0.05.O106$filename)

MLMSUPER0.05.O106$filename <- gsub('MLMSUPER0.05.O106/Trait.', '', MLMSUPER0.05.O106$filename)
MLMSUPER0.05.O106$filename <- gsub('.csv', '', MLMSUPER0.05.O106$filename)

levels(as.factor(MLMSUPER0.05.O106$filename))

MLMSUPER0.05.O106.1 <- unite(MLMSUPER0.05.O106,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMSUPER0.05.O106.flo <- merge(flo.2,MLMSUPER0.05.O106.1,by="SNP.T")
MLMSUPER0.05.O106.flo.S <- separate(MLMSUPER0.05.O106.flo, SNP.T, c("SNP","Trait.name"),sep="/")

##for OWA 116
###for OWA 116
###for OWA 116
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYO.116.4232.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDO.116.4232.csv")
flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="O116+4232" & flo.1$Method=="MLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.O116")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.O116/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMSUPER0.05.O116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMSUPER0.05.O116$filename)

MLMSUPER0.05.O116$filename <- gsub('MLMSUPER0.05.O116/Trait.', '', MLMSUPER0.05.O116$filename)
MLMSUPER0.05.O116$filename <- gsub('.csv', '', MLMSUPER0.05.O116$filename)

levels(as.factor(MLMSUPER0.05.O116$filename))

MLMSUPER0.05.O116.1 <- unite(MLMSUPER0.05.O116,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMSUPER0.05.O116.flo <- merge(flo.2,MLMSUPER0.05.O116.1,by="SNP.T")
MLMSUPER0.05.O116.flo.S <- separate(MLMSUPER0.05.O116.flo, SNP.T, c("SNP","Trait.name"),sep="/")



##for OWA 126
###for OWA 126
###for OWA 126
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYO.126.3659.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDO.126.3659.csv")
flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="O126+3659" & flo.1$Method=="MLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.O126")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.O126/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMSUPER0.05.O126 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMSUPER0.05.O126$filename)

MLMSUPER0.05.O126$filename <- gsub('MLMSUPER0.05.O126/Trait.', '', MLMSUPER0.05.O126$filename)
MLMSUPER0.05.O126$filename <- gsub('.csv', '', MLMSUPER0.05.O126$filename)

levels(as.factor(MLMSUPER0.05.O126$filename))

MLMSUPER0.05.O126.1 <- unite(MLMSUPER0.05.O126,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMSUPER0.05.O126.flo <- merge(flo.2,MLMSUPER0.05.O126.1,by="SNP.T")
MLMSUPER0.05.O126.flo.S <- separate(MLMSUPER0.05.O126.flo, SNP.T, c("SNP","Trait.name"),sep="/")


##for OWA 135
###for OWA 135
###for OWA 135
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYO.135.2855.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDO.135.2855.csv")
flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="O135+2855" & flo.1$Method=="MLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.O135")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.O135/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMSUPER0.05.O135 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMSUPER0.05.O135$filename)

MLMSUPER0.05.O135$filename <- gsub('MLMSUPER0.05.O135/Trait.', '', MLMSUPER0.05.O135$filename)
MLMSUPER0.05.O135$filename <- gsub('.csv', '', MLMSUPER0.05.O135$filename)

levels(as.factor(MLMSUPER0.05.O135$filename))

MLMSUPER0.05.O135.1 <- unite(MLMSUPER0.05.O135,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMSUPER0.05.O135.flo <- merge(flo.2,MLMSUPER0.05.O135.1,by="SNP.T")
MLMSUPER0.05.O135.flo.S <- separate(MLMSUPER0.05.O135.flo, SNP.T, c("SNP","Trait.name"),sep="/")
###combine all of the result
MLMSUPER0.05.PVE.a <- do.call("rbind",list(MLMSUPER0.05.C106.flo.S,MLMSUPER0.05.C116.flo.S,MLMSUPER0.05.C125.flo.S,MLMSUPER0.05.O106.flo.S,
                                        MLMSUPER0.05.O116.flo.S,MLMSUPER0.05.O126.flo.S,MLMSUPER0.05.O135.flo.S,MLMCUPERim.flo.S,
                                        MLMSUPER0.05.F106.flo.S,MLMSUPER0.05.F116.flo.S,MLMSUPER0.05.F122.flo.S))

write.csv(MLMSUPER0.05.PVE.a,file="Result.8.3/MLMSUPER0.05.PS.PVE.a .csv")
###combine 
MLMCMLMSUPER0.05.PVE.a <- do.call("rbind", list(MLMSUPER0.05.PVE.a, CMLMSUPER0.05.PVE.a))

write.csv(MLMCMLMSUPER0.05.PVE.a,file="Result.8.3/MLMCMLMSUPER0.05.PS.PVE.a.csv")

