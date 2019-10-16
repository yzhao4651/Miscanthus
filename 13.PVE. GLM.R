###import the data with GLM as method 
flo.maf.t.7.14 <- read.csv("data/flo.maf.t.7.14.csv")
str(flo.maf.t.7.14)
flo <- flo.maf.t.7.14[flo.maf.t.7.14$Method=="GLM",]
levels(flo$Trait.name)[levels(flo$Trait.name)=="HD_50"] <- "HD_50."
levels(flo$Trait.name)[levels(flo$Trait.name)=="FD_50"] <- "FD_50."
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GLM.7.15.106im")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait)>0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

###import data 
flo.maf.t.7.14 <- read.csv("data/flo.maf.t.7.14.csv")
str(flo.maf.t.7.14)
flo <- flo.maf.t.7.14[flo.maf.t.7.14$Method=="GLM" & flo.maf.t.7.14$filename=="106",]
levels(flo$Trait.name)[levels(flo$Trait.name)=="HD_50"] <- "HD_50."
levels(flo$Trait.name)[levels(flo$Trait.name)=="FD_50"] <- "FD_50."
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GLM.7.15.106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait)>0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

flo.maf.t.7.14 <- read.csv("data/flo.maf.t.7.14.csv")
str(flo.maf.t.7.14)
flo <- flo.maf.t.7.14[flo.maf.t.7.14$Method=="GLM" & flo.maf.t.7.14$filename=="116",]
levels(flo$Trait.name)[levels(flo$Trait.name)=="HD_50"] <- "HD_50."
levels(flo$Trait.name)[levels(flo$Trait.name)=="FD_50"] <- "FD_50."
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GLM.7.15.116")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait)>0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

flo.maf.t.7.14 <- read.csv("data/flo.maf.t.7.14.csv")
str(flo.maf.t.7.14)
flo <- flo.maf.t.7.14[flo.maf.t.7.14$Method=="GLM" & flo.maf.t.7.14$filename=="122",]
levels(flo$Trait.name)[levels(flo$Trait.name)=="HD_50"] <- "HD_50."
levels(flo$Trait.name)[levels(flo$Trait.name)=="FD_50"] <- "FD_50."
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GLM.7.15.122")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait)>0) {
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
    ftdGDSNPs2 <- ftdGDSNPs[,-1]
    rownames(ftdGDSNPs2) <- ftdGDSNPs[,1]
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
fileNames <- Sys.glob(paste("GLM.7.15.106im/*.", extension, sep = ""))

myY <- read.csv("data/myYimputedSNP19.csv")
myGD <- read.csv("data/myGDimputedSNP19.csv")

mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
resultPVEF116im <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)




setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("GLM.7.15.106/*.", extension, sep = ""))
myY <- read.csv("data/myYf.106.4185.csv")
myGD <- read.csv("data/myGDf.106.4185.csv")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
resultPVEF106 <- plyr::rbind.fill(mzList)

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("GLM.7.15.116/*.", extension, sep = ""))
myY <- read.csv("data/myYf.116.3077.csv")
myGD <- read.csv("data/myGDf.116.3077.csv")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
resultPVEF116 <- plyr::rbind.fill(mzList)

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("GLM.7.15.122/*.", extension, sep = ""))
myY <- read.csv("data/myYf.122.2272.csv")
myGD <- read.csv("data/myGDf.122.2272.csv")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
resultPVEF122 <- plyr::rbind.fill(mzList)

resultF.T <- do.call("rbind",list(resultPVEF116im,resultPVEF106,resultPVEF116,resultPVEF122))

PVE.T.7.15.F.S <- do.call("rbind",list(resultF.T,result.T))
write.csv(PVE.T.7.15.F.S,file="PVE/Flo.GLM.PVE.T.csv")



