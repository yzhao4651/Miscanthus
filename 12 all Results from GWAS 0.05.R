##for flowering data each one of methods were separated 
###check ECMMLCV and ECMMLNOCV
##Result from CMLMM method 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
extension <- "csv"
fileNames <- Sys.glob(paste("GAPIT005/*ECMMLCV0.05/*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value <  0.05)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
ECMMLCVa <- plyr::ldply(mzList, data.frame)
str(ECMMLCVa)
levels(ECMMLCVa$filename)
ECMMLCVa$filename <- gsub('GAPIT005/', '', ECMMLCVa$filename)
ECMMLCVa$filename <- gsub('.GWAS.Results.csv', '', ECMMLCVa$filename)
ECMMLCVa$filename <- gsub('ECMMLCV0.05', '.PS_Y', ECMMLCVa$filename)
ECMMLCVa$filename <- gsub('/GAPIT.', '.GAPIT.', ECMMLCVa$filename)
levels(as.factor(ECMMLCVa$filename))
##install packages to get function unit
#install.packages(c("tidyr", "devtools"))
library(tidyr)
ECMMLCVa.1 <- separate(ECMMLCVa, filename, c("Ind.SNP","PS","Software","Method","Trait.name"), "\\.")
str(ECMMLCVa.1)
levels(as.factor(ECMMLCVa$nobs))
ECMMLCVa.1$Ind.SNP <- gsub("^$", '+36088im', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("C125", 'C125+2562', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("F116", 'F116+3077', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("F122", 'F122+2272', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("O106", 'O106+4834', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("O116", 'O116+4232', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("O126", 'O126+3659', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("O135", 'O135+2855', ECMMLCVa.1$Ind.SNP)
str(ECMMLCVa.1)
ECMMLCVa.1$PVE <- (ECMMLCVa.1$Rsquare.of.Model.with.SNP-ECMMLCVa.1$Rsquare.of.Model.without.SNP)*100
#ECMMLCVa.2 <- ECMMLCVa.1[,c(1:15)]
#str(ECMMLCVa.2)
str(ECMMLCVa.1)
levels(as.factor(ECMMLCVa.1$Ind.SNP))
write.csv(ECMMLCVa.1, file="Result.8.3/ECMML0.05.PS.csv")


##ECMMLNOCV without CV
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
extension <- "csv"
fileNames <- Sys.glob(paste("GAPIT005/*ECMML0.05/*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value <  0.05)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
ECMMLNOCVa <-  plyr::ldply(mzList, data.frame)
str(ECMMLNOCVa)
levels(ECMMLNOCVa$filename)
ECMMLNOCVa$filename <- gsub('GAPIT005/', '', ECMMLNOCVa$filename)
ECMMLNOCVa$filename <- gsub('.GWAS.Results.csv', '', ECMMLNOCVa$filename)
ECMMLNOCVa$filename <- gsub('ECMML0.05', '.PS_N', ECMMLNOCVa$filename)
ECMMLNOCVa$filename <- gsub('/GAPIT.', '.GAPIT.', ECMMLNOCVa$filename)
levels(as.factor(ECMMLNOCVa$filename))
ECMMLNOCVa.1 <- separate(ECMMLNOCVa, filename, c("Ind.SNP","PS","Software","Method","Trait.name"), "\\.")
str(ECMMLNOCVa.1)
levels(as.factor(ECMMLNOCVa$nobs))
ECMMLNOCVa.1$Ind.SNP <- gsub("^$", '+36088im', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("C125", 'C125+2562', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("F116", 'F116+3077', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("F122", 'F122+2272', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("O106", 'O106+4834', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("O116", 'O116+4232', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("O126", 'O126+3659', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("O135", 'O135+2855', ECMMLNOCVa.1$Ind.SNP)
str(ECMMLNOCVa.1)
ECMMLNOCVa.1$PVE <- (ECMMLNOCVa.1$Rsquare.of.Model.with.SNP-ECMMLNOCVa.1$Rsquare.of.Model.without.SNP)*100
write.csv(ECMMLNOCVa.1, file="Result.8.3/ECMML0.05.Adj.P.PSNO.csv")
ECMMLNOCVa.2 <- ECMMLNOCVa.1[,c(1:6,9:15,17)]
str(ECMMLNOCVa.2)
write.csv(ECMMLNOCVa.2, file="Result.8.3/ECMML.PSNO.csv")


##MLMSUPPER and CMLMMSUPPER with CV 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
#setwd("GAPIT005/MLMSUPPER/")
extension <- "csv"
fileNames <- Sys.glob(paste("GAPIT005/*MLMSUPPERCV0.05/*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value <  0.05)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
SUPERCVa <-  plyr::ldply(mzList, data.frame)
str(SUPERCVa)
levels(SUPERCVa$filename)
str(SUPERCVa)
levels(SUPERCVa$filename)
SUPERCVa$filename <- gsub('GAPIT005/', '', SUPERCVa$filename)
SUPERCVa$filename <- gsub('.GWAS.Results.csv', '', SUPERCVa$filename)
SUPERCVa$filename <- gsub('CMLMSUPPERCV0.05', '.CMLM+SUPER.PS_Y', SUPERCVa$filename)
SUPERCVa$filename <- gsub('MLMSUPPERCV0.05', '.MLM+SUPER.PS_Y', SUPERCVa$filename)
SUPERCVa$filename <- gsub('/GAPIT.', '.GAPIT.', SUPERCVa$filename)
levels(as.factor(SUPERCVa$filename))
str(SUPERCVa.1)
SUPERCVa.1 <- separate(SUPERCVa, filename, c("Ind.SNP","Method","PS","Software","Method2","Trait.name"), "\\.")
levels(as.factor(SUPERCVa.1$nobs))
levels(as.factor(SUPERCVa.1$Ind.SNP))
str(SUPERCVa.1)

SUPERCVa.1$Ind.SNP <- gsub("C125", 'C125+2562', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("F116", 'F116+3077', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("F122", 'F122+2272', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("O106", 'O106+4834', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("O116", 'O116+4232', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("O126", 'O126+3659', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("O135", 'O135+2855', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("^$", "+36088im", SUPERCVa.1$Ind.SNP)

levels(as.factor(SUPERCVa.1$Ind.SNP))
str(SUPERCVa.1)
write.csv(SUPERCVa.1, file="Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
#SUPERCVa.1$PVE <- (SUPERCVa.1$Rsquare.of.Model.with.SNP-SUPERCVa.1$Rsquare.of.Model.without.SNP)*100
SUPERCVa.2 <- SUPERCVa.1[,c(1:6,9:16)]
str(SUPERCVa.2)

write.csv(SUPERCVa.2, file="Result.8.3/MLMCMLMSUPER.PS.csv")




##MLMSUPPER and CMLM NO CV 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
#setwd("GAPIT005/MLMSUPPER/")
extension <- "csv"
fileNames <- Sys.glob(paste("GAPIT005/*MLMSUPPER0.05/*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value <  0.05)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
SUPERNOCVa <- plyr::ldply(mzList, data.frame)
#MLMCMLMSUPPERa <- do.call("rbind", mzList)
str(SUPERNOCVa)
levels(SUPERNOCVa$filename)
str(SUPERNOCVa)
levels(SUPERNOCVa$filename)
SUPERNOCVa$filename <- gsub('GAPIT005/', '', SUPERNOCVa$filename)
SUPERNOCVa$filename <- gsub('.GWAS.Results.csv', '', SUPERNOCVa$filename)
SUPERNOCVa$filename <- gsub('CMLMSUPPER0.05', '.CMLM+SUPER.PS_N', SUPERNOCVa$filename)
SUPERNOCVa$filename <- gsub('MLMSUPPER0.05', '.MLM+SUPER.PS_N', SUPERNOCVa$filename)
SUPERNOCVa$filename <- gsub('/GAPIT.', '.GAPIT.', SUPERNOCVa$filename)
levels(as.factor(SUPERNOCVa$filename))
str(SUPERNOCVa)
SUPERNOCVa.1 <- separate(SUPERNOCVa, filename, c("Ind.SNP","Method","PS","Software","Method2","Trait.name"), "\\.")
levels(as.factor(SUPERNOCVa.1$nobs))
levels(as.factor(SUPERNOCVa.1$Ind.SNP))

SUPERNOCVa.1$Ind.SNP <- gsub("C125", 'C125+2562', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("F116", 'F116+3077', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("F122", 'F122+2272', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("O106", 'O106+4834', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("O116", 'O116+4232', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("O126", 'O126+3659', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("O135", 'O135+2855', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("^$", "+36088im", SUPERNOCVa.1$Ind.SNP)

levels(as.factor(SUPERNOCVa.1$Ind.SNP))
str(SUPERNOCVa.1)
#SUPERNOCVa.1$PVE <- (SUPERNOCVa.1$Rsquare.of.Model.with.SNP-SUPERNOCVa.1$Rsquare.of.Model.without.SNP)*100

write.csv(SUPERNOCVa.1, file="Result.8.3/MLMCMLMSUPER0.05.Adj.P.PSNO.csv")

SUPERNOCVa.2 <- SUPERNOCVa.1[,c(1:6,9:16)]
str(SUPERNOCVa.2)
write.csv(SUPERNOCVa.2, file="Result.8.3/MLMCMLMSUPER.PSNO.csv")



###GLM and MLM and MLMM with CV
###GLM with CV
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
extension <- "csv"
fileNames <- Sys.glob(paste("GAPIT005/*GLMMLMMLMMCV0.05/GAPIT.*LM.*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value <  0.05)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
GLMMLMMLMMCVa <- plyr::ldply(mzList, data.frame)
str(GLMMLMMLMMCVa)
levels(GLMMLMMLMMCVa$filename)
levels(as.factor(GLMMLMMLMMCVa$nobs))

GLMMLMMLMMCVa$filename <- gsub('GAPIT005/', '', GLMMLMMLMMCVa$filename)
GLMMLMMLMMCVa$filename <- gsub('.GWAS.Results.csv', '', GLMMLMMLMMCVa$filename)
GLMMLMMLMMCVa$filename <- gsub('GLMMLMMLMMCV0.05', '.PS_Y', GLMMLMMLMMCVa$filename)
GLMMLMMLMMCVa$filename <- gsub('/GAPIT.', '.GAPIT.', GLMMLMMLMMCVa$filename)
levels(as.factor(GLMMLMMLMMCVa$filename))
str(GLMMLMMLMMCVa)
GLMMLMMLMMCVa.1 <- separate(GLMMLMMLMMCVa, filename, c("Ind.SNP","PS","Software","Method","Trait.name"), "\\.")
levels(as.factor(GLMMLMMLMMCVa.1$nobs))
levels(as.factor(GLMMLMMLMMCVa.1$Ind.SNP))

GLMMLMMLMMCVa.1$Ind.SNP <- gsub("C125", 'C125+2562', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("F116", 'F116+3077', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("F122", 'F122+2272', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("O106", 'O106+4834', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("O116", 'O116+4232', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("O126", 'O126+3659', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("O135", 'O135+2855', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("^$", "+36088im", GLMMLMMLMMCVa.1$Ind.SNP)

levels(as.factor(GLMMLMMLMMCVa.1$Ind.SNP))
str(GLMMLMMLMMCVa.1)
GLMMLMMLMMCVa.1$PVE <- (GLMMLMMLMMCVa.1$Rsquare.of.Model.with.SNP-GLMMLMMLMMCVa.1$Rsquare.of.Model.without.SNP)*100
write.csv(GLMMLMMLMMCVa.1, file="Result.8.3/GLMMLMMLMM0.05.Adj.P.PS.csv")
GLMMLMMLMMCVa.2 <- GLMMLMMLMMCVa.1[,c(1:6,9:15,17)]
str(GLMMLMMLMMCVa.2)
write.csv(GLMMLMMLMMCVa.2, file="Result.8.3/GLMMLMMLMM.PS.csv")


###GLM and MLM, MLMMwithout CV
###withiut CV
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
extension <- "csv"
fileNames <- Sys.glob(paste("GAPIT005/*GLMMLMMLMM0.05/GAPIT.*LM.*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value <  0.05)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
GLMMLMMLMMNOCVa <- plyr::ldply(mzList, data.frame)
str(GLMMLMMLMMNOCVa)
levels(GLMMLMMLMMNOCVa$filename)
levels(as.factor(GLMMLMMLMMNOCVa$nobs))

GLMMLMMLMMNOCVa$filename <- gsub('GAPIT005/', '', GLMMLMMLMMNOCVa$filename)
GLMMLMMLMMNOCVa$filename <- gsub('.GWAS.Results.csv', '', GLMMLMMLMMNOCVa$filename)
GLMMLMMLMMNOCVa$filename <- gsub('GLMMLMMLMM0.05', '.PS_N', GLMMLMMLMMNOCVa$filename)
GLMMLMMLMMNOCVa$filename <- gsub('/GAPIT.', '.GAPIT.', GLMMLMMLMMNOCVa$filename)
levels(as.factor(GLMMLMMLMMNOCVa$filename))
str(GLMMLMMLMMNOCVa)
GLMMLMMLMMNOCVa.1 <- separate(GLMMLMMLMMNOCVa, filename, c("Ind.SNP","PS","Software","Method","Trait.name"), "\\.")
levels(as.factor(GLMMLMMLMMNOCVa.1$nobs))
levels(as.factor(GLMMLMMLMMNOCVa.1$Ind.SNP))

GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("C125", 'C125+2562', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("F116", 'F116+3077', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("F122", 'F122+2272', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("O106", 'O106+4834', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("O116", 'O116+4232', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("O126", 'O126+3659', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("O135", 'O135+2855', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("^$", "+36088im", GLMMLMMLMMNOCVa.1$Ind.SNP)

levels(as.factor(GLMMLMMLMMNOCVa.1$Ind.SNP))
str(GLMMLMMLMMNOCVa.1)
GLMMLMMLMMNOCVa.1$PVE <- (GLMMLMMLMMNOCVa.1$Rsquare.of.Model.with.SNP-GLMMLMMLMMNOCVa.1$Rsquare.of.Model.without.SNP)*100
write.csv(GLMMLMMLMMNOCVa.1, file="Result.8.3/GLMMLMMLMM0.05.Adj.P.PSNO.csv")
GLMMLMMLMMNOCVa.2 <- GLMMLMMLMMNOCVa.1[,c(1:6,9:15,17)]
str(GLMMLMMLMMNOCVa.2)
write.csv(GLMMLMMLMMNOCVa.2, file="Result.8.3/GLMMLMMLMM.PSNO.csv")



##FarmCPU with CV
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
extension <- "csv"
fileNames <- Sys.glob(paste("GAPIT005/*FarmCPUCV0.05/FarmCPU.*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value < 0.05 )
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}

#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
FarmCPUCVa <- plyr::ldply(mzList, data.frame)
str(FarmCPUCVa)
levels(FarmCPUCVa$filename)
levels(FarmCPUCVa$nobs)

FarmCPUCVa$filename <- gsub('GAPIT005/', '', FarmCPUCVa$filename)
FarmCPUCVa$filename <- gsub('.GWAS.Results.csv', '', FarmCPUCVa$filename)
FarmCPUCVa$filename <- gsub('FarmCPUCV0.05', '.PS_Y', FarmCPUCVa$filename)
FarmCPUCVa$filename <- gsub('/FarmCPU.', '.FarmCPU.FarmCPU.', FarmCPUCVa$filename)
levels(as.factor(FarmCPUCVa$filename))
str(FarmCPUCVa)
FarmCPUCVa.1 <- separate(FarmCPUCVa, filename, c("Ind.SNP","PS","Software","Method","Trait.name"), "\\.")
levels(as.factor(FarmCPUCVa.1$bs))
levels(as.factor(FarmCPUCVa.1$Ind.SNP))

FarmCPUCVa.1$Ind.SNP <- gsub("C125", 'C125+2562', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("F116", 'F116+3077', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("F122", 'F122+2272', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("O106", 'O106+4834', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("O116", 'O116+4232', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("O126", 'O126+3659', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("O135", 'O135+2855', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("^$", "+36088im", FarmCPUCVa.1$Ind.SNP)

levels(as.factor(FarmCPUCVa.1$Ind.SNP))
str(FarmCPUCVa.1)
#FarmCPUCVa.1$PVE <- (FarmCPUCVa.1$Rsquare.of.Model.with.SNP-FarmCPUCVa.1$Rsquare.of.Model.without.SNP)*100
write.csv(FarmCPUCVa.1,file="Result.8.3/FarmCPU0.05.Adj.P.PS.csv")


FarmCPUCVa.2 <- FarmCPUCVa.1[,c(1:6,8:12)]
str(FarmCPUCVa.2)

write.csv(FarmCPUCVa.2,file="Result.8.3/FarmCPU.PS.csv")


##FarmCPU NO CV 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
extension <- "csv"
fileNames <- Sys.glob(paste("GAPIT005/*FarmCPU0.05/FarmCPU.*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value < 0.05 )
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
FarmCPUNOCVa <- plyr::ldply(mzList, data.frame)
str(FarmCPUNOCVa)
levels(FarmCPUNOCVa$filename)
levels(as.factor(FarmCPUNOCVa$nobs))

FarmCPUNOCVa$filename <- gsub('GAPIT005/', '', FarmCPUNOCVa$filename)
levels(as.factor(FarmCPUNOCVa$filename))
FarmCPUNOCVa$filename <- gsub('.GWAS.Results.csv', '', FarmCPUNOCVa$filename)
levels(as.factor(FarmCPUNOCVa$filename))
FarmCPUNOCVa$filename <- gsub('FarmCPU0.05', '.PS_N', FarmCPUNOCVa$filename)
levels(as.factor(FarmCPUNOCVa$filename))
FarmCPUNOCVa$filename <- gsub('/', '.FarmCPU.', FarmCPUNOCVa$filename)
levels(as.factor(FarmCPUNOCVa$filename))
str(FarmCPUNOCVa)
FarmCPUNOCVa.1 <- separate(FarmCPUNOCVa, filename, c("Ind.SNP","PS","Software","Method","Trait.name"), "\\.")
levels(as.factor(FarmCPUNOCVa.1$nobs))
levels(as.factor(FarmCPUNOCVa.1$Ind.SNP))

FarmCPUNOCVa.1$Ind.SNP <- gsub("C125", 'C125+2562', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("F116", 'F116+3077', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("F122", 'F122+2272', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("O106", 'O106+4834', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("O116", 'O116+4232', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("O126", 'O126+3659', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("O135", 'O135+2855', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("^$", "+36088im", FarmCPUNOCVa.1$Ind.SNP)

levels(as.factor(FarmCPUNOCVa.1$Ind.SNP))
str(FarmCPUNOCVa.1)
write.csv(FarmCPUNOCVa.1,file="Result.8.3/FarmCPU0.05.Adj.P.PSNO.csv")
#FarmCPUNOCVa.1$PVE <- (FarmCPUNOCVa.1$Rsquare.of.Model.with.SNP-FarmCPUNOCVa.1$Rsquare.of.Model.without.SNP)*100
FarmCPUNOCVa.2 <- FarmCPUNOCVa.1[,c(1:5,7:8,10:14)]
str(FarmCPUNOCVa.2)
write.csv(FarmCPUNOCVa.2,file="Result.8.3/FarmCPU.PSNO.csv")



###rmlmm with CV
###rmlmm with CV
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/Resultfloall1.1/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$X.log10.P.!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMCVFim <- do.call("rbind", mzList)
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait13"] <- "TFN"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait14"] <- "FNMain"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait15"] <- "FNsmall"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait19"] <- "HD_1"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait20"] <- "FD_1"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait21"] <- "HD_50"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait22"] <- "FD_50"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait23"] <- "HW_1"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait24"] <- "FW_1"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait25"] <- "HW_50"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait26"] <- "FW_50"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait27"] <- "GHW_1"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait28"] <- "GFW_1"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait29"] <- "GHW_50"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait30"] <- "GFW_50"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait31"] <- "HM_1"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait32"] <- "FM_1"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait33"] <- "HM_50"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait34"] <- "FM_50"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait35"] <- "fprind"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait36"] <- "fprinW"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait37"] <- "fprinGW"
levels(rMLMMCVFim$Trait.name)[levels(rMLMMCVFim$Trait.name)=="Trait38"] <- "fprinM"
rMLMMCVFim <- droplevels(rMLMMCVFim)
str(rMLMMCVFim)

###f106,116,122 rmlmm
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/F*.1/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$X.log10.P.!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMCVfnm <- do.call("rbind", mzList)
str(rMLMMCVfnm)
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait1"] <- "TFN"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait2"] <- "FNMain"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait3"] <- "FNsmall"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait4"] <- "HD_1"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait5"] <- "FD_1"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait6"] <- "HD_50"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait7"] <- "FD_50"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait8"] <- "HW_1"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait9"] <- "FW_1"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait10"] <- "HW_50"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait11"] <- "FW_50"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait12"] <- "GHW_1"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait13"] <- "GFW_1"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait14"] <- "GHW_50"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait15"] <- "GFW_50"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait16"] <- "HM_1"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait17"] <- "FM_1"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait18"] <- "HM_50"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait19"] <- "FM_50"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait20"] <- "fprind"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait21"] <- "fprinW"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait22"] <- "fprinGW"
levels(rMLMMCVfnm$Trait.name)[levels(rMLMMCVfnm$Trait.name)=="Trait23"] <- "fprinM"

str(rMLMMCVfnm)
rMLMMCVfloall <- do.call("rbind",list(rMLMMCVFim,rMLMMCVfnm))
levels(rMLMMCVfloall$filename)
levels(as.factor(rMLMMCVfloall$nobs))

rMLMMCVfloall$filename <- gsub('mrMLMM2', 'mrMLMM/2', rMLMMCVfloall$filename)
rMLMMCVfloall$filename <- gsub('\\.1', '/PS_Y', rMLMMCVfloall$filename)
#rMLMMCVfloall$filename <- gsub('Resultfloall1', 'F116+36088im', rMLMMCVfloall$filename)
#rMLMMCVfloall$filename <- gsub('/GAPIT.', '.GAPIT.', rMLMMCVfloall$filename)
levels(as.factor(rMLMMCVfloall$filename))
str(rMLMMCVfloall)
rMLMMCVfloall.1 <- separate(rMLMMCVfloall, filename, c("Software","2","Ind.SNP","PS","discard"), "/")
levels(as.factor(rMLMMCVfloall.1$nobs))
levels(as.factor(rMLMMCVfloall.1$Ind.SNP))

rMLMMCVfloall.1$Ind.SNP <- gsub("C125", 'C125+2562', rMLMMCVfloall.1$Ind.SNP)
rMLMMCVfloall.1$Ind.SNP <- gsub("C106", 'C106+4202', rMLMMCVfloall.1$Ind.SNP)
rMLMMCVfloall.1$Ind.SNP <- gsub("C116", 'C116+3293', rMLMMCVfloall.1$Ind.SNP)
rMLMMCVfloall.1$Ind.SNP <- gsub("F106", 'F106+4185', rMLMMCVfloall.1$Ind.SNP)
rMLMMCVfloall.1$Ind.SNP <- gsub("F116", 'F116+3077', rMLMMCVfloall.1$Ind.SNP)
rMLMMCVfloall.1$Ind.SNP <- gsub("F122", 'F122+2272', rMLMMCVfloall.1$Ind.SNP)
rMLMMCVfloall.1$Ind.SNP <- gsub("Resultfloall1", "F116+36088im", rMLMMCVfloall.1$Ind.SNP)

levels(as.factor(rMLMMCVfloall.1$Ind.SNP))
str(rMLMMCVfloall.1)
colnames(rMLMMCVfloall.1)[colnames(rMLMMCVfloall.1)=="RS."] <- "SNP"
colnames(rMLMMCVfloall.1)[colnames(rMLMMCVfloall.1)=="Marker.Position..bp."] <- "Position"
rMLMMCVfloall.1$P.value <- 10^-(rMLMMCVfloall.1$X.log10.P.)
colnames(rMLMMCVfloall.1)[colnames(rMLMMCVfloall.1)=="r2...."] <- "PVE"
str(rMLMMCVfloall.1)

#rMLMMCVfloall.1$PVE <- (rMLMMCVfloall.1$Rsquare.of.Model.with.SNP-rMLMMCVfloall.1$Rsquare.of.Model.without.SNP)*100
rMLMMCVfloall.2 <- rMLMMCVfloall.1[,c(4:6,2:3,7:8,10:15,17:18,20)]
str(rMLMMCVfloall.2)
write.csv(rMLMMCVfloall.2,file="Result.8.3/rMLMMflo.PS.csv")


####  flo NO NO CV
###rmlmm without CV
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/Resultfloall1.2/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$X.log10.P.!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMNOCVFim <- do.call("rbind", mzList)
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait13"] <- "TFN"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait14"] <- "FNMain"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait15"] <- "FNsmall"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait19"] <- "HD_1"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait20"] <- "FD_1"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait21"] <- "HD_50"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait22"] <- "FD_50"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait23"] <- "HW_1"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait24"] <- "FW_1"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait25"] <- "HW_50"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait26"] <- "FW_50"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait27"] <- "GHW_1"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait28"] <- "GFW_1"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait29"] <- "GHW_50"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait30"] <- "GFW_50"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait31"] <- "HM_1"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait32"] <- "FM_1"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait33"] <- "HM_50"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait34"] <- "FM_50"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait35"] <- "fprind"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait36"] <- "fprinW"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait37"] <- "fprinGW"
levels(rMLMMNOCVFim$Trait.name)[levels(rMLMMNOCVFim$Trait.name)=="Trait38"] <- "fprinM"
rMLMMNOCVFim <- droplevels(rMLMMNOCVFim)
str(rMLMMNOCVFim)


###f106, 116, 122 rmlmm NO CV
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/F*.2/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$X.log10.P.!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMNOCVFnm <- do.call("rbind", mzList)
str(rMLMMNOCVFnm)

levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait1"] <- "TFN"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait2"] <- "FNMain"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait3"] <- "FNsmall"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait4"] <- "HD_1"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait5"] <- "FD_1"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait6"] <- "HD_50"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait7"] <- "FD_50"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait8"] <- "HW_1"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait9"] <- "FW_1"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait10"] <- "HW_50"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait11"] <- "FW_50"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait12"] <- "GHW_1"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait13"] <- "GFW_1"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait14"] <- "GHW_50"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait15"] <- "GFW_50"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait16"] <- "HM_1"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait17"] <- "FM_1"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait18"] <- "HM_50"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait19"] <- "FM_50"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait20"] <- "fprind"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait21"] <- "fprinW"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait22"] <- "fprinGW"
levels(rMLMMNOCVFnm$Trait.name)[levels(rMLMMNOCVFnm$Trait.name)=="Trait23"] <- "fprinM"

str(rMLMMNOCVFnm)
rMLMMNOCVfloall <- do.call("rbind",list(rMLMMNOCVFim,rMLMMNOCVFnm))

str(rMLMMNOCVfloall)


str(rMLMMNOCVfloall)
levels(rMLMMNOCVfloall$filename)
levels(as.factor(rMLMMNOCVfloall$nobs))

rMLMMNOCVfloall$filename <- gsub('mrMLMM2', 'mrMLMM/2', rMLMMNOCVfloall$filename)
rMLMMNOCVfloall$filename <- gsub('\\.2', '/PS_N', rMLMMNOCVfloall$filename)
#rMLMMNOCVfloall$filename <- gsub('Resultfloall1', 'F116+36088im', rMLMMNOCVfloall$filename)
#rMLMMNOCVfloall$filename <- gsub('/GAPIT.', '.GAPIT.', rMLMMNOCVfloall$filename)
levels(as.factor(rMLMMNOCVfloall$filename))
str(rMLMMNOCVfloall)
rMLMMNOCVfloall.1 <- separate(rMLMMNOCVfloall, filename, c("Software","2","Ind.SNP","PS","discard"), "/")
levels(as.factor(rMLMMNOCVfloall.1$nobs))
levels(as.factor(rMLMMNOCVfloall.1$Ind.SNP))

rMLMMNOCVfloall.1$Ind.SNP <- gsub("C125", 'C125+2562', rMLMMNOCVfloall.1$Ind.SNP)
rMLMMNOCVfloall.1$Ind.SNP <- gsub("C106", 'C106+4202', rMLMMNOCVfloall.1$Ind.SNP)
rMLMMNOCVfloall.1$Ind.SNP <- gsub("C116", 'C116+3293', rMLMMNOCVfloall.1$Ind.SNP)
rMLMMNOCVfloall.1$Ind.SNP <- gsub("F106", 'F106+4185', rMLMMNOCVfloall.1$Ind.SNP)
rMLMMNOCVfloall.1$Ind.SNP <- gsub("F116", 'F116+3077', rMLMMNOCVfloall.1$Ind.SNP)
rMLMMNOCVfloall.1$Ind.SNP <- gsub("F122", 'F122+2272', rMLMMNOCVfloall.1$Ind.SNP)
rMLMMNOCVfloall.1$Ind.SNP <- gsub("Resultfloall1", "F116+36088im", rMLMMNOCVfloall.1$Ind.SNP)

levels(as.factor(rMLMMNOCVfloall.1$Ind.SNP))
str(rMLMMNOCVfloall.1)
colnames(rMLMMNOCVfloall.1)[colnames(rMLMMNOCVfloall.1)=="RS."] <- "SNP"
colnames(rMLMMNOCVfloall.1)[colnames(rMLMMNOCVfloall.1)=="Marker.Position..bp."] <- "Position"
rMLMMNOCVfloall.1$P.value <- 10^-(rMLMMNOCVfloall.1$X.log10.P.)
colnames(rMLMMNOCVfloall.1)[colnames(rMLMMNOCVfloall.1)=="r2...."] <- "PVE"
str(rMLMMNOCVfloall.1)
#rMLMMNOCVfloall.1$PVE <- (rMLMMNOCVfloall.1$Rsquare.of.Model.with.SNP-rMLMMNOCVfloall.1$Rsquare.of.Model.without.SNP)*100
rMLMMNOCVfloall.2 <- rMLMMNOCVfloall.1[,c(4:6,2:3,7:8,10:15,17:18,20)]
str(rMLMMNOCVfloall.2)
write.csv(rMLMMNOCVfloall.2,file="Result.8.3/rMLMMflo.PSNO.csv")




###rmlmm  CV SRD with imputed SNP
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/ResultSRDall1.1/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$X.log10.P.!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMCVSRDim <- do.call("rbind", mzList)
levels(rMLMMCVSRDim$Trait.name)[levels(rMLMMCVSRDim$Trait.name)=="Trait1"] <- "SRD"
str(rMLMMCVSRDim)


###rmlmm CV OWA with imputed SNP
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/ResultOWA3all1.1/1_Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$X.log10.P.!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMCVOWAim <- do.call("rbind", mzList)
levels(rMLMMCVOWAim$Trait.name)[levels(rMLMMCVOWAim$Trait.name)=="Trait1"] <- "OWA"
str(rMLMMCVOWAim)


###rmlmm CV OWA with non missing value
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/OWA*.1/1_Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$X.log10.P.!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMCVOWAnm <- do.call("rbind", mzList)
levels(rMLMMCVOWAnm$Trait.name)[levels(rMLMMCVOWAnm$Trait.name)=="Trait1"] <- "OWA"
str(rMLMMCVOWAnm)


###rmlmm Clum CV with imputed SNPs
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/Resultculmupdate.1/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$X.log10.P.!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMCulmCVim <- do.call("rbind", mzList)
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait2"] <- "CmDW_g"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait3"] <- "Cml_cm"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait4"] <- "CmD_BI_mm"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait5"] <- "CmD_LI_mm"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait6"] <- "CmN." 
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait7"] <- "Bcirc_cm"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait8"] <- "Yld_kg"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait9"] <- "SDW_kg"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait10"] <- "CCirc_cm"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait11"] <- "Lg"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait12"] <- "GS"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait16"] <- "FD"
#levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait17"] <- "SRD"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait18"] <- "ADD"
rMLMMCulmCVim <- droplevels(rMLMMCulmCVim)
str(rMLMMCulmCVim)

###Clum CV with non missing
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/C*.1/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  colnames(sample)[colnames(sample)=="Marker.position..bp."] <- "Marker.Position..bp."
  if (nrow(sample)>0){
    mz.idx = which(sample$X.log10.P.!= Inf)
    mz1 = sample[mz.idx, ]
    mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
  }
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMCVCulmnm <- plyr::ldply(mzList, data.frame)
str(rMLMMCVCulmnm)

levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait1"] <- "CmDW_g"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait2"] <- "Cml_cm"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait3"] <- "CmD_BI_mm"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait4"] <- "CmD_LI_mm"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait5"] <- "CmN." 
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait6"] <- "Bcirc_cm"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait7"] <- "Yld_kg"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait8"] <- "SDW_kg"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait9"] <- "CCirc_cm"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait10"] <- "Lg"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait11"] <- "GS"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait12"] <- "FD"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait13"] <- "SRD"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait14"] <- "ADD"
str(rMLMMCVCulmnm)

###combina all of the Results together 

rMLMMClumCV <- plyr::ldply(list(rMLMMCulmCVim,rMLMMCVCulmnm,rMLMMCVOWAim,rMLMMCVOWAnm,rMLMMCVSRDim), data.frame)

str(rMLMMClumCV)
levels(rMLMMClumCV$filename)
levels(as.factor(rMLMMClumCV$nobs))

rMLMMClumCV$filename <- gsub('mrMLMM2', 'mrMLMM/2', rMLMMClumCV$filename)
rMLMMClumCV$filename <- gsub('\\.1', '/PS_Y', rMLMMClumCV$filename)
#rMLMMClumCV$filename <- gsub('Resultfloall1', 'F116+36088im', rMLMMClumCV$filename)
#rMLMMClumCV$filename <- gsub('/GAPIT.', '.GAPIT.', rMLMMClumCV$filename)
levels(as.factor(rMLMMClumCV$filename))
str(rMLMMClumCV)
rMLMMClumCV.1 <- separate(rMLMMClumCV, filename, c("Software","2","Ind.SNP","PS","discard"), "/")
levels(as.factor(rMLMMClumCV.1$nobs))
levels(as.factor(rMLMMClumCV.1$Ind.SNP))

rMLMMClumCV.1$Ind.SNP <- gsub("C125", 'C125+2562', rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("C106", 'C106+4202', rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("C116", 'C116+3293', rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("f106", 'F106+4185', rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("f116", 'F116+3077', rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("f122", 'F122+2272', rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("Resultculmupdate", "124+36088im", rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("ResultSRDall1", "122+35256im", rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("ResultOWA3all1", "143+35256im", rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("OWA106", "O106+4834", rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("OWA116", "O116+4232", rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("OWA126", "O126+3659", rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("OWA135", "O135+2855", rMLMMClumCV.1$Ind.SNP)

levels(as.factor(rMLMMClumCV.1$Ind.SNP))
str(rMLMMClumCV.1)
colnames(rMLMMClumCV.1)[colnames(rMLMMClumCV.1)=="RS."] <- "SNP"
colnames(rMLMMClumCV.1)[colnames(rMLMMClumCV.1)=="Marker.Position..bp."] <- "Position"
rMLMMClumCV.1$P.value <- 10^-(rMLMMClumCV.1$X.log10.P.)
colnames(rMLMMClumCV.1)[colnames(rMLMMClumCV.1)=="r2...."] <- "PVE"
str(rMLMMClumCV.1)
#rMLMMClumCV.1$PVE <- (rMLMMClumCV.1$Rsquare.of.Model.with.SNP-rMLMMClumCV.1$Rsquare.of.Model.without.SNP)*100
rMLMMClumCV.2 <- rMLMMClumCV.1[,c(4:6,2:3,7:8,10:15,17:18,21)]
str(rMLMMClumCV.2)
write.csv(rMLMMClumCV.2,file="Result.8.3/rMLMMCulm.PS.csv")



###rmlmm NO CV SRD
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/ResultSRDall1.2/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$X.log10.P.!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMSRDNOCVim <- do.call("rbind", mzList)

levels(rMLMMSRDNOCVim$Trait.name)[levels(rMLMMSRDNOCVim$Trait.name)=="Trait1"] <- "SRD"
str(rMLMMSRDNOCVim)


###rmlmm NO CV OWA with imputed SNP
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/ResultOWA3all1.2/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$X.log10.P.!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMOWANOCVim <- do.call("rbind", mzList)

levels(rMLMMOWANOCVim$Trait.name)[levels(rMLMMOWANOCVim$Trait.name)=="Trait1"] <- "OWA"
str(rMLMMOWANOCVim)


###rmlmm NO CV OWA with non missing value
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/OWA*.2/1_Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$X.log10.P.!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMOWANOCVnm <- do.call("rbind", mzList)
levels(rMLMMOWANOCVnm$Trait.name)[levels(rMLMMOWANOCVnm$Trait.name)=="Trait1"] <- "OWA"
str(rMLMMOWANOCVnm)


###rmlmm Clum NO CV with imputed SNPs
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/Resultculmupdate.2/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  colnames(sample)[colnames(sample)=="Marker.position..bp."] <- "Marker.Position..bp."
  mz.idx = which(sample$X.log10.P.!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMCulmNOCVim <- do.call("rbind", mzList)
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait2"] <- "CmDW_g"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait3"] <- "Cml_cm"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait4"] <- "CmD_BI_mm"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait5"] <- "CmD_LI_mm"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait6"] <- "CmN." 
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait7"] <- "Bcirc_cm"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait8"] <- "Yld_kg"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait9"] <- "SDW_kg"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait10"] <- "CCirc_cm"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait11"] <- "Lg"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait12"] <- "GS"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait16"] <- "FD"
#levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait17"] <- "SRD"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait18"] <- "ADD"
rMLMMCulmNOCVim <- droplevels(rMLMMCulmNOCVim)
str(rMLMMCulmNOCVim)


###Clum NO CV non missing
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/C*.2/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  colnames(sample)[colnames(sample)=="Marker.position..bp."] <- "Marker.Position..bp."
  mz.idx = which(sample$X.log10.P.!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMNOCVCulmnm <- plyr::ldply(mzList, data.frame)
str(rMLMMNOCVCulmnm)
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait1"] <- "CmDW_g"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait2"] <- "Cml_cm"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait3"] <- "CmD_BI_mm"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait4"] <- "CmD_LI_mm"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait5"] <- "CmN." 
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait6"] <- "Bcirc_cm"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait7"] <- "Yld_kg"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait8"] <- "SDW_kg"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait9"] <- "CCirc_cm"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait10"] <- "Lg"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait11"] <- "GS"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait12"] <- "FD"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait13"] <- "SRD"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait14"] <- "ADD"
str(rMLMMNOCVCulmnm)

###Combine all of the results NO CV
rMLMMClumNOCV <- plyr::ldply(list(rMLMMCulmNOCVim,rMLMMNOCVCulmnm,rMLMMOWANOCVim,rMLMMOWANOCVnm,rMLMMSRDNOCVim), data.frame)
str(rMLMMClumNOCV)
levels(rMLMMClumNOCV$filename)
levels(as.factor(rMLMMClumNOCV$nobs))

rMLMMClumNOCV$filename <- gsub('mrMLMM2', 'mrMLMM/2', rMLMMClumNOCV$filename)
rMLMMClumNOCV$filename <- gsub('\\.2', '/PS_N', rMLMMClumNOCV$filename)
#rMLMMClumNOCV$filename <- gsub('Resultfloall1', 'F116+36088im', rMLMMClumNOCV$filename)
#rMLMMClumNOCV$filename <- gsub('/GAPIT.', '.GAPIT.', rMLMMClumNOCV$filename)
levels(as.factor(rMLMMClumNOCV$filename))
str(rMLMMClumNOCV)
rMLMMClumNOCV.1 <- separate(rMLMMClumNOCV, filename, c("Software","2","Ind.SNP","PS","discard"), "/")
levels(as.factor(rMLMMClumNOCV.1$nobs))
levels(as.factor(rMLMMClumNOCV.1$Ind.SNP))

rMLMMClumNOCV.1$Ind.SNP <- gsub("C125", 'C125+2562', rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("C106", 'C106+4202', rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("C116", 'C116+3293', rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("f106", 'F106+4185', rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("f116", 'F116+3077', rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("f122", 'F122+2272', rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("Resultculmupdate", "124+36088im", rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("ResultSRDall1", "122+35256im", rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("ResultOWA3all1", "143+35256im", rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("OWA106", "O106+4834", rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("OWA116", "O116+4232", rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("OWA126", "O126+3659", rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("OWA135", "O135+2855", rMLMMClumNOCV.1$Ind.SNP)

levels(as.factor(rMLMMClumNOCV.1$Ind.SNP))
str(rMLMMClumNOCV.1)
colnames(rMLMMClumNOCV.1)[colnames(rMLMMClumNOCV.1)=="RS."] <- "SNP"
colnames(rMLMMClumNOCV.1)[colnames(rMLMMClumNOCV.1)=="Marker.Position..bp."] <- "Position"
rMLMMClumNOCV.1$P.value <- 10^-(rMLMMClumNOCV.1$X.log10.P.)
colnames(rMLMMClumNOCV.1)[colnames(rMLMMClumNOCV.1)=="r2...."] <- "PVE"
str(rMLMMClumNOCV.1)
#rMLMMClumNOCV.1$PVE <- (rMLMMClumNOCV.1$Rsquare.of.Model.with.SNP-rMLMMClumNOCV.1$Rsquare.of.Model.without.SNP)*100
rMLMMClumNOCV.2 <- rMLMMClumNOCV.1[,c(4:6,2:3,7:8,10:15,17:18,21)]
str(rMLMMClumNOCV.2)
write.csv(rMLMMClumNOCV.2,file="Result.8.3/rMLMMCulm.PSNO.csv")



###result from rrBLUP
###result from rrBLUP
####import result from RRBLUP all of imputed SNPs
rrBLUPimall <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
str(rrBLUPimall)
names(rrBLUPimall)
rrBLUPimall.1 <- rrBLUPimall[, -c(4,20,42,58)]

str(rrBLUPimall.1)
rrBLUPimall.2 <- rrBLUPimall.1[rowSums(rrBLUPimall.1[4:39])!=0,]
names(rrBLUPimall.2 )

rrBLUPimall.3 <- rrBLUPimall.2[, c(1:3,40:75)]
str(rrBLUPimall.3)
str(rrBLUPimall.3)
names(rrBLUPimall.3)

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(rrBLUPimall.3, 4,39)
str(plot.p.adjusted)
data <- plot.p.adjusted
rrblupList = list()
for(i in 4:39){
  rrBLUP <- data[,c(1:3,i)][data[[i+36]] < 0.05, ]
  if(nrow(rrBLUP)>=1){
    rrblupList[[i]] <- data.frame(rrBLUP, Trait.name = rep(colnames(data[i]), nrow(rrBLUP)))
    colnames(rrblupList[[i]])[colnames(rrblupList[[i]])== colnames(data[i]) ] <- "P.value"
  }
}
rrblupSNPim <- plyr::ldply(rrblupList, data.frame)
rrblupSNPim$filename  <- paste("rrBLUP143+36088im/")
####import result from RRBLUP for flower 106. 116.122

extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPF*/rrBLUP_GWAS_results.flo*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample <- sample[rowSums(sample[4:26])!=0,]
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  plot.p.adjusted <- adj_P_function(sample, 27, 49)
  data <- plot.p.adjusted[,c(1:3,27:72)]
  rrblupList = list()
  for(j in 4:26){
    rrBLUP <- data[,c(1:3,j)][data[[j+23]] < 0.05, ]
    if(nrow(rrBLUP)>=1){
      rrblupList[[j]] <- data.frame(rrBLUP, Trait.name = rep(colnames(data[j]), nrow(rrBLUP)))
      colnames(rrblupList[[j]])[colnames(rrblupList[[j]])== colnames(data[j]) ] <- "P.value"
    }
  }
  rrblupSNPfnon <- plyr::ldply(rrblupList, data.frame)
  #rrblupSNPfnon <- plyr::rbind.fill(rrblupList)
  mzList[[i]] = data.frame(rrblupSNPfnon, filename = rep(fileNames[i], nrow(rrblupSNPfnon)))  
}
rrblupSNPfnm <- plyr::ldply(mzList, data.frame)

####import result from OWA
####import result from OWA
extension <- "csv"
#rrBLUPimall <- read.csv("rrBLUPOWA/rrBLUP_GWAS_results.*.csv",row.names = 1)
fileNames <- Sys.glob(paste("rrBLUPOWA/rrBLUP_GWAS_results.*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample <- sample[rowSums(sample[4])!=0,]
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  plot.p.adjusted <- adj_P_function(sample, 5, 5)
  data <- plot.p.adjusted[,c(1:3,5:6)]
  rrblupList = list()
  for(j in 4:4){
    rrBLUP <- data[,c(1:3,j)][data[[j+1]] < 0.05, ]
    if(nrow(rrBLUP)>=1){
      rrblupList[[j]] <- data.frame(rrBLUP, Trait.name = rep(colnames(data[j]), nrow(rrBLUP)))
      colnames(rrblupList[[j]])[colnames(rrblupList[[j]])== colnames(data[j]) ] <- "P.value"
    }
  }
  rrblupSNPfnon <- plyr::ldply(rrblupList, data.frame)
  #rrblupSNPfnon <- plyr::rbind.fill(rrblupList)
  mzList[[i]] = data.frame(rrblupSNPfnon, filename = rep(fileNames[i], nrow(rrblupSNPfnon)))  
}
rrblupSNPOWA <- plyr::ldply(mzList, data.frame)


####import result from SRD
####import result from SRD
extension <- "csv"
#rrBLUPimall <- read.csv("rrBLUPOWA/rrBLUP_GWAS_results.*.csv",row.names = 1)
fileNames <- Sys.glob(paste("rrBLUPSRD/rrBLUP_GWAS_results.*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample <- sample[rowSums(sample[4])!=0,]
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  plot.p.adjusted <- adj_P_function(sample, 5, 5)
  data <- plot.p.adjusted[,c(1:3,5:6)]
  rrblupList = list()
  for(j in 4:4){
    rrBLUP <- data[,c(1:3,j)][data[[j+1]] < 0.05, ]
    if(nrow(rrBLUP)>=1){
      rrblupList[[j]] <- data.frame(rrBLUP, Trait.name = rep(colnames(data[j]), nrow(rrBLUP)))
      colnames(rrblupList[[j]])[colnames(rrblupList[[j]])== colnames(data[j]) ] <- "P.value"
    }
  }
  rrblupSNPfnon <- plyr::ldply(rrblupList, data.frame)
  #rrblupSNPfnon <- plyr::rbind.fill(rrblupList)
  mzList[[i]] = data.frame(rrblupSNPfnon, filename = rep(fileNames[i], nrow(rrblupSNPfnon)))  
}
rrblupSNPSRD <- plyr::ldply(mzList, data.frame)

####import result from Culm
####import result from Culm
extension <- "csv"
#rrBLUPimall <- read.csv("rrBLUPOWA/rrBLUP_GWAS_results.*.csv",row.names = 1)
fileNames <- Sys.glob(paste("rrBLUPC*/rrBLUP_GWAS_results.*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample <- sample[rowSums(sample[,c(4:17)])!=0,]
  names(sample)
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  plot.p.adjusted <- adj_P_function(sample, 18, 31)
  data <- plot.p.adjusted[,c(1:3,18:45)]
  names(data)
  rrblupList = list()
  for(j in 4:17){
    rrBLUP <- data[,c(1:3,j)][data[[j+14]] < 0.05, ]
    if(nrow(rrBLUP)>=1){
      rrblupList[[j]] <- data.frame(rrBLUP, Trait.name = rep(colnames(data[j]), nrow(rrBLUP)))
      colnames(rrblupList[[j]])[colnames(rrblupList[[j]])== colnames(data[j]) ] <- "P.value"
    }
  }
  rrblupSNPfnon <- plyr::ldply(rrblupList, data.frame)
  #rrblupSNPfnon <- plyr::rbind.fill(rrblupList)
  mzList[[i]] = data.frame(rrblupSNPfnon, filename = rep(fileNames[i], nrow(rrblupSNPfnon)))  
}
rrblupSNPClumnm <- plyr::ldply(mzList, data.frame)

####combine all of results from rrBLUP together 
rrBLUPalltrait <- plyr::ldply(list(rrblupSNPClumnm,rrblupSNPfnm,rrblupSNPim,rrblupSNPOWA,rrblupSNPSRD),data.frame)

levels(rrBLUPalltrait$Trait.name)
levels(as.factor(rrBLUPalltrait$filename))

rrBLUPalltrait.1 <- separate(rrBLUPalltrait, Trait.name, c("Method","Trait.name"), "\\.")
str(rrBLUPalltrait.1)

rrBLUPalltrait.1$filename <- gsub('rrBLUP', 'rrBLUP/', rrBLUPalltrait.1$filename)
levels(as.factor(rrBLUPalltrait.1$filename))

rrBLUPalltrait.2 <- separate(rrBLUPalltrait.1, filename, c("Software","Ind.SNP","discard","discard2"), "/")
str(rrBLUPalltrait.2)
rrBLUPalltrait.2$Ind.SNP <- gsub('C125', 'C125+2562', rrBLUPalltrait.2$Ind.SNP)
rrBLUPalltrait.2$Ind.SNP <- gsub('C130', 'C130+2001', rrBLUPalltrait.2$Ind.SNP)
rrBLUPalltrait.2$Ind.SNP <- gsub('C116', 'C116+3293', rrBLUPalltrait.2$Ind.SNP)
rrBLUPalltrait.2$Ind.SNP <- gsub('C106', 'C106+4202', rrBLUPalltrait.2$Ind.SNP)
rrBLUPalltrait.2$Ind.SNP <- gsub('F106', 'F106+4185', rrBLUPalltrait.2$Ind.SNP)
str(rrBLUPalltrait.2)
colnames(rrBLUPalltrait.2)[colnames(rrBLUPalltrait.2)=="Name"] <- "SNP"
rrBLUPalltrait.3 <- rrBLUPalltrait.2[1:8]
rrBLUPalltrait.3$PS <- paste("PCA=3")
str(rrBLUPalltrait.3)
write.csv(rrBLUPalltrait.3,file="Result.8.3/rrBLUPalltrait.PCA.csv")

###selet the flo traits
###selet the flo traits
##selet the flo traits
##import



