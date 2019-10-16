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
flo.1 <- read.csv("Result.8.3/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
levels(flo.1$Ind.SNP)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
levels(flo.1$Trait.name)
flo <- flo.1[flo.1$Ind.SNP=="+36088im" & flo.1$Method=="FarmCPU",]
levels(flo$Ind.SNP)
str(flo)
names(flo)
flo <- droplevels(flo)
levels(flo$Ind.SNP)
levels(as.factor(flo$Method))
levels(as.factor(flo$Trait.name))
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="fprin"),]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)

#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.7.15.116im")
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.7.25.116im")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.36088im")
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
fileNames <- Sys.glob(paste("FarmCPUCV0.05.36088im/*.", extension, sep = ""))

#myY <- read.csv("data/myYimputedSNP19.csv")
myGD <- read.csv("data/myGDimputedSNP19.csv")

mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCVim <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCVim$filename)

FarmCPUCVim$filename <- gsub('FarmCPUCV0.05.36088im/Trait.', '', FarmCPUCVim$filename)
FarmCPUCVim$filename <- gsub('.csv', '', FarmCPUCVim$filename)

levels(as.factor(FarmCPUCVim$filename))

FarmCPUCVim.1 <- unite(FarmCPUCVim,"SNP.T",SNP,filename,sep="/")
flo.2<- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCVim.flo <- merge(flo.2,FarmCPUCVim.1,by="SNP.T")
FarmCPUCVim.flo.S <- separate(FarmCPUCVim.flo, SNP.T, c("SNP","Trait.name"),sep="/")



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
flo.1 <- read.csv("Result.8.3/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="C106+4202" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="fprin"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.C106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.C106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV.C106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV.C106$filename)

FarmCPUCV.C106$filename <- gsub('FarmCPUCV0.05.C106/Trait.', '', FarmCPUCV.C106$filename)
FarmCPUCV.C106$filename <- gsub('.csv', '', FarmCPUCV.C106$filename)

levels(as.factor(FarmCPUCV.C106$filename))

FarmCPUCV.C106.1 <- unite(FarmCPUCV.C106,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV.C106.flo <- merge(flo.2,FarmCPUCV.C106.1,by="SNP.T")
FarmCPUCV.C106.flo.S <- separate(FarmCPUCV.C106.flo, SNP.T, c("SNP","Trait.name"),sep="/")


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
names(myY)
myGD <- read.csv("data/myGDC.116.3293.csv")
flo.1 <- read.csv("Result.8.3/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="C116+3293" & flo.1$Method=="FarmCPU",]
str(flo)
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="fprin"),]

flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.C116")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.C116/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  str(sample)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV.C116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV.C116$filename)

FarmCPUCV.C116$filename <- gsub('FarmCPUCV0.05.C116/Trait.', '', FarmCPUCV.C116$filename)
FarmCPUCV.C116$filename <- gsub('.csv', '', FarmCPUCV.C116$filename)

levels(as.factor(FarmCPUCV.C116$filename))

FarmCPUCV.C116.1 <- unite(FarmCPUCV.C116,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV.C116.flo <- merge(flo.2,FarmCPUCV.C116.1,by="SNP.T")
FarmCPUCV.C116.flo.S <- separate(FarmCPUCV.C116.flo, SNP.T, c("SNP","Trait.name"),sep="/")

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
flo.1 <- read.csv("Result.8.3/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="C125+2562" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="fprin"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.C125")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.C125/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV.C125 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV.C125$filename)

FarmCPUCV.C125$filename <- gsub('FarmCPUCV0.05.C125/Trait.', '', FarmCPUCV.C125$filename)
FarmCPUCV.C125$filename <- gsub('.csv', '', FarmCPUCV.C125$filename)

levels(as.factor(FarmCPUCV.C125$filename))

FarmCPUCV.C125.1 <- unite(FarmCPUCV.C125,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV.C125.flo <- merge(flo.2,FarmCPUCV.C125.1,by="SNP.T")
FarmCPUCV.C125.flo.S <- separate(FarmCPUCV.C125.flo, SNP.T, c("SNP","Trait.name"),sep="/")


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
flo.1 <- read.csv("Result.8.3/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"

levels(flo.1$Ind.SNP)
flo <- flo.1[flo.1$Ind.SNP=="O135+2855" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="fprin"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.O135")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.O135/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV.O135 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV.O135$filename)

FarmCPUCV.O135$filename <- gsub('FarmCPUCV0.05.O135/Trait.', '', FarmCPUCV.O135$filename)
FarmCPUCV.O135$filename <- gsub('.csv', '', FarmCPUCV.O135$filename)

levels(as.factor(FarmCPUCV.O135$filename))

FarmCPUCV.O135.1 <- unite(FarmCPUCV.O135,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV.O135.flo <- merge(flo.2,FarmCPUCV.O135.1,by="SNP.T")
FarmCPUCV.O135.flo.S <- separate(FarmCPUCV.O135.flo, SNP.T, c("SNP","Trait.name"),sep="/")
###combine all of the result


FarmCPUCV.PVE.a <- do.call("rbind",list(FarmCPUCV.C106.flo.S,FarmCPUCV.C116.flo.S,FarmCPUCV.C125.flo.S,
                                        FarmCPUCVim.flo.S,FarmCPUCV.O135.flo.S))
write.csv(FarmCPUCV.PVE.a,file="Result.8.3/FarmCPU0.05.PS.PVE.a.1.csv")
FarmCPU0.05.PVE.PS.PSNO.a <- do.call("rbind",list(FarmCPUCV.PVE.a, FarmCPU0.05.PVE.a))
write.csv(FarmCPU0.05.PVE.PS.PSNO.a,file="Result.8.3/FarmCPU0.05.PVE.PS.PSNO.a.1.csv")


###combine 
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
flo.1 <- read.csv("Result.8.3/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="F106+4185" & flo.1$Method=="FarmCPU",]

flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="fprin"),]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.F106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.F106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV.F106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV.F106$filename)

FarmCPUCV.F106$filename <- gsub('FarmCPUCV0.05.F106/Trait.', '', FarmCPUCV.F106$filename)
FarmCPUCV.F106$filename <- gsub('.csv', '', FarmCPUCV.F106$filename)

levels(as.factor(FarmCPUCV.F106$filename))

FarmCPUCV.F106.1 <- unite(FarmCPUCV.F106,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV.F106.flo <- merge(flo.2,FarmCPUCV.F106.1,by="SNP.T")
FarmCPUCV.F106.flo.S <- separate(FarmCPUCV.F106.flo, SNP.T, c("SNP","Trait.name"),sep="/")

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

flo.1 <- read.csv("Result.8.3/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"

flo <- flo.1[flo.1$Ind.SNP=="F116+3077" & flo.1$Method=="FarmCPU",]

flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="fprin"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.F116")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.F116/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV.F116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV.F116$filename)

FarmCPUCV.F116$filename <- gsub('FarmCPUCV0.05.F116/Trait.', '', FarmCPUCV.F116$filename)
FarmCPUCV.F116$filename <- gsub('.csv', '', FarmCPUCV.F116$filename)

levels(as.factor(FarmCPUCV.F116$filename))

FarmCPUCV.F116.1 <- unite(FarmCPUCV.F116,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV.F116.flo <- merge(flo.2,FarmCPUCV.F116.1,by="SNP.T")
FarmCPUCV.F116.flo.S <- separate(FarmCPUCV.F116.flo, SNP.T, c("SNP","Trait.name"),sep="/")




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


flo.1 <- read.csv("Result.8.3/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="F122+2272" & flo.1$Method=="FarmCPU",]


flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="fprin"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.F122")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.F122/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV.F122 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV.F122$filename)

FarmCPUCV.F122$filename <- gsub('FarmCPUCV0.05.F122/Trait.', '', FarmCPUCV.F122$filename)
FarmCPUCV.F122$filename <- gsub('.csv', '', FarmCPUCV.F122$filename)

levels(as.factor(FarmCPUCV.F122$filename))

FarmCPUCV.F122.1 <- unite(FarmCPUCV.F122,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV.F122.flo <- merge(flo.2,FarmCPUCV.F122.1,by="SNP.T")
FarmCPUCV.F122.flo.S <- separate(FarmCPUCV.F122.flo, SNP.T, c("SNP","Trait.name"),sep="/")

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
flo.1 <- read.csv("Result.8.3/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"
flo <- flo.1[flo.1$Ind.SNP=="O106+4834" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="fprin"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.O106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.O106/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV.O106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV.O106$filename)

FarmCPUCV.O106$filename <- gsub('FarmCPUCV0.05.O106/Trait.', '', FarmCPUCV.O106$filename)
FarmCPUCV.O106$filename <- gsub('.csv', '', FarmCPUCV.O106$filename)

levels(as.factor(FarmCPUCV.O106$filename))

FarmCPUCV.O106.1 <- unite(FarmCPUCV.O106,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV.O106.flo <- merge(flo.2,FarmCPUCV.O106.1,by="SNP.T")
FarmCPUCV.O106.flo.S <- separate(FarmCPUCV.O106.flo, SNP.T, c("SNP","Trait.name"),sep="/")

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
flo.1 <- read.csv("Result.8.3/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"

flo <- flo.1[flo.1$Ind.SNP=="O116+4232" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="fprin"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.O116")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.O116/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV.O116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV.O116$filename)

FarmCPUCV.O116$filename <- gsub('FarmCPUCV0.05.O116/Trait.', '', FarmCPUCV.O116$filename)
FarmCPUCV.O116$filename <- gsub('.csv', '', FarmCPUCV.O116$filename)

levels(as.factor(FarmCPUCV.O116$filename))

FarmCPUCV.O116.1 <- unite(FarmCPUCV.O116,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV.O116.flo <- merge(flo.2,FarmCPUCV.O116.1,by="SNP.T")
FarmCPUCV.O116.flo.S <- separate(FarmCPUCV.O116.flo, SNP.T, c("SNP","Trait.name"),sep="/")

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
flo.1 <- read.csv("Result.8.3/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
colnames(flo.1)[colnames(flo.1)=="Name"] <- "SNP"

flo <- flo.1[flo.1$Ind.SNP=="O126+3659" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="fprin"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.O126")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.O126/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV.O126 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV.O126$filename)

FarmCPUCV.O126$filename <- gsub('FarmCPUCV0.05.O126/Trait.', '', FarmCPUCV.O126$filename)
FarmCPUCV.O126$filename <- gsub('.csv', '', FarmCPUCV.O126$filename)

levels(as.factor(FarmCPUCV.O126$filename))

FarmCPUCV.O126.1 <- unite(FarmCPUCV.O126,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV.O126.flo <- merge(flo.2,FarmCPUCV.O126.1,by="SNP.T")
FarmCPUCV.O126.flo.S <- separate(FarmCPUCV.O126.flo, SNP.T, c("SNP","Trait.name"),sep="/")



