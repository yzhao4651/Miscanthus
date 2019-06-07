### select the SNP with adjusted_p.values <= 0.05
##set up the path way for result from ECMMLCV
### CMLM(k+Q) model
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
setwd("Result GAPIT1/ECMMLCV/")
extension <- "csv"
fileNames <- Sys.glob(paste("*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$FDR_Adjusted_P.values <= 0.05 & sample$maf >= 0.01)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
resultECMMLCV1 = do.call("rbind", mzList)
str(resultECMMLCV1)
levels(resultECMMLCV1$filename)
write.csv(resultECMMLCV1, file="Result/resultECMMLCV1.csv")

### MLM(k+Q) model
##set up the path way for result from ImputedSNPCV 
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
setwd("Result GAPIT1/ImputedSNPCV/")
extension <- "csv"
fileNames <- Sys.glob(paste("*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$FDR_Adjusted_P.values <= 0.05 & sample$maf >= 0.01)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
resultImputedSNPCV <- do.call("rbind", mzList)
write.csv(resultImputedSNPCV, file="Result/resultImputedSNPCV.csv")

##set up the path way for result from CMLMSUPPER
##CMLMSUPPER(K+Q)
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
setwd("Result GAPIT1/CMLMSUPPER/")
extension <- "csv"
fileNames <- Sys.glob(paste("*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$FDR_Adjusted_P.values <= 0.05 & sample$maf >= 0.01)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
resultCMLMSUPPER <- do.call("rbind", mzList)
str(resultCMLMSUPPER)
write.csv(resultCMLMSUPPER, file="Result/resultCMLMSUPPER.csv")


##set up the path way for result from MLMSUPPER
##MLMSUPPER
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
setwd("Result GAPIT1/MLMSUPPER/")
extension <- "csv"
fileNames <- Sys.glob(paste("*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  mz.idx = which(sample$FDR_Adjusted_P.values <= 0.05 & sample$maf >= 0.01)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
resultMLMSUPPER <- do.call("rbind", mzList)
str(resultMLMSUPPER)
write.csv(resultMLMSUPPER, file="Result/resultMLMSUPPER.csv")




