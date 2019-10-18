
###import all the dataset with MAF 0.05
EPS <- read.csv("Result.8.3/ECMML0.05.PS.csv",row.names=1)

str(EPS)
#chang maf to MAF
colnames(EPS)[colnames(EPS)=="maf"] <- "MAF"
str(EPS)
#EPSNO <- read.csv("Result.8.3/ECMML0.05.Adj.P.PSNO.csv",row.names=1)
#colnames(EPSNO)[colnames(EPSNO)=="maf"] <- "MAF"
#str(EPSNO)
GPS <- read.csv("Result.8.3/GLMMLMMLMM0.05.Adj.P.PS.csv",row.names=1)
colnames(GPS)[colnames(GPS)=="maf"] <- "MAF"
str(GPS)
rMFPS <- read.csv("Result.8.3/rMLMMflo.PS.csv",row.names=1)
str(rMFPS)
rMFPSNO <- read.csv("Result.8.3/rMLMMflo.PSNO.csv",row.names=1)
str(rMFPSNO)
rMCPS <- read.csv("Result.8.3/rMLMMCulm.PS.csv",row.names=1)
str(rMFPS)
rMCPSNO <- read.csv("Result.8.3/rMLMMCulm.PSNO.csv",row.names=1)
str(rMFPSNO)
rrBLUP <- read.csv("Result.8.3/rrBLUPalltrait.PCA.csv",row.names=1)
str(rrBLUP)
FarmPS <- read.csv("Result.8.3/FarmCPU0.05.PVE.PS.PSNO.a.1.csv",row.names=1)
colnames(FarmPS)[colnames(FarmPS)=="maf"] <- "MAF"
str(FarmPS)

MLMCMLM <- read.csv("Result.8.3/MLMCMLMSUPER0.05.PS.PVE.a.csv",row.names=1)
colnames(MLMCMLM)[colnames(MLMCMLM)=="maf"] <- "MAF"
str(MLMCMLM)

###combineing all of them together 
alltrait.method <- plyr::ldply(list(EPS,GPS,rMFPS,rMFPSNO,rMCPS,rMCPSNO,rrBLUP,FarmPS,MLMCMLM))
str(alltrait.method)
###
###change the P.value to with zero decime 
alltrait.method$P.value.1 <- as.numeric(formatC(alltrait.method$P.value, format = "e", digits = 0))
alltrait.method$MAF.1<- as.numeric(formatC(alltrait.method$MAF, format = "e", digits = 3))
alltrait.method$PVE.1 <- as.numeric(formatC(alltrait.method$PVE, format = "e", digits = 3))
##write out the total traits 
write.csv(alltrait.method, file="data/alltrait.updated2.csv")

###select all the flowering trait

###import the all traits 
all <- read.csv("Result.8.3/Alltrait.1.csv")
levels(all$Trait.name)

floall <- all[all$Trait.name =="HD_1" | all$Trait.name =="HD_50" | 
                               all$Trait.name =="FD_1" | all$Trait.name =="FD_50" |
                               all$Trait.name =="FM_1" | all$Trait.name =="FM_50" |
                               all$Trait.name =="HM_1" | all$Trait.name =="HM_50" |
                               all$Trait.name =="FW_1" | all$Trait.name =="FW_50" |
                               all$Trait.name =="HW_1" | all$Trait.name =="HW_50" |
                               all$Trait.name =="GHW_1" | all$Trait.name =="GHW_50" |
                               all$Trait.name =="GFW_1" | all$Trait.name =="GFW_50" |
                               all$Trait.name =="fprind" | all$Trait.name =="fprinGW" |
                               all$Trait.name =="fprinM" | all$Trait.name =="fprinW",]

###write out the flowering trait 

write.csv(floall,file="Result.8.3/flo.8.3.1.updated.csv")

###select the OWA

OWAall <- all[all$Trait.name =="OWA",] 
str(OWAall)
##write out the owa data set 

write.csv(OWAall,file="OWA/owa.csv")

### select the largest MAF and PVE

OWAall <- OWAall[order(OWAall$SNP, -(as.numeric(OWAall$MAF.1))),] #sort by id and reverse of abs(value)
OWAall.MAF.2 <- OWAall[ !duplicated(OWAall$SNP), ] 
OWAall <- OWAall[order(OWAall$SNP, -(as.numeric(OWAall$PVE.1))), ] #sort by id and reverse of abs(value)
OWAall.PVE.2 <- OWAall[ !duplicated(OWAall$SNP), ] 
##merge this two 
OWA.MAF.PVE.Max <- merge(x=OWAall.MAF.2,y=OWAall.PVE.2, by="SNP" )
str(OWA.MAF.PVE.Max)
OWA.MAF.PVE.Max.1 <- OWA.MAF.PVE.Max[,c(1,7,47)]

OWA.MAF.PVE.Max.1$ MAF.x <- as.numeric(formatC(OWA.MAF.PVE.Max.1$ MAF.x, format = "e", digits = 3))
OWA.MAF.PVE.Max.1$ PVE.y <- as.numeric(formatC(OWA.MAF.PVE.Max.1$ PVE.y, format = "e", digits = 3))

write.csv(OWA.MAF.PVE.Max.1,file="OWA/OWA.MAF.PVE.Max.1.csv")





  
  
  
  
  
  
  
  
  
  