####GAPIT and FarmCPU
####import the dataset need for software 
###imported all the data need for GAPIT and FarmCPU
myY <- read.csv("data/myYimputedSNP19.csv")
myGD <- read.csv("data/myGDimputedSNP19.csv")
myGM <- read.csv("data/myGMimputedSNP19.csv")
myQ<- read.csv("data/myQimputedSNP19.csv")

###install the packaged need for GAPIT
source("http://www.bioconductor.org/biocLite.R")
biocLite("multtest")
install.packages("gplots")
install.packages("LDheatmap")
install.packages("genetics")
install.packages("ape")
install.packages("EMMREML")
install.packages("scatterplot3d") #The downloaded link at: http://cran.r-project.org/package=scatterplot3d
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")
####Basic Scenario with Model selection in order to find the optimal number of PCs
setwd("~/Documents/GAPITandFarmCPU/GAPITnm")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ
  PCA.total=3, 
  Model.selection = TRUE
)
####Basic Scenario of Compressed MLM without CV 
setwd("~/Documents/GAPITandFarmCPU/GAPITnmnocv")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM
)
####Basic Scenario of Supper method with CV
setwd("~/Documents/GAPITandFarmCPU/GAPITCVCMLM")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER
    sangwich.bottom="CMLM" #options are GLM,MLM,CMLM, FaST and SUPER
)
####Basic Scenario of Supper method with CV
setwd("~/Documents/GAPITandFarmCPU/GAPITnmCVsuper")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER
    sangwich.bottom="SUPER" #options are GLM,MLM,CMLM, FaST and SUPER
  )
###########Using ECMLM
  setwd("~/Documents/GAPITandFarmCPU/myGAPITCVecmlm")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    kinship.cluster=c("average", "complete", "ward"),
    kinship.group=c("Mean", "Max"),
    group.from=200,
    group.to=1000000,
    group.by=10
  )
####Basic Scenario of CMLM method without CV
  setwd("~/Documents/GAPITandFarmCPU/GAPITCMLM")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER
    sangwich.bottom="CMLM" #options are GLM,MLM,CMLM, FaST and SUPER
  )
  
  ###########Using ECMLM
  setwd("~/Documents/GAPITandFarmCPU/myGAPITecmlm")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    kinship.cluster=c("average", "complete", "ward"),
    kinship.group=c("Mean", "Max"),
    group.from=200,
    group.to=1000000,
    group.by=10
  )
  setwd("~/Documents/GAPITandFarmCPU/GAPITnmsuper")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    sangwich.top="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER
    sangwich.bottom="SUPER" #options are GLM,MLM,CMLM, FaST and SUPER
  )
  
  #FarmCPU
  #FarmCPU
  #FarmCPU
  install.packages("bigmemory")
  install.packages("biganalytics")
  library("bigmemory")
  library("biganalytics")
  library("bigmemory")
  library("biganalytics")
  library("compiler") #this library is already installed in R
  source("http://zzlab.net/GAPIT/gapit_functions.txt")
  source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
  

setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmhday") 
myFarmCPU <- FarmCPU(
    Y=myY[,1:2],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmfday")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,3)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmhhday")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,4)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmhfday")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,5)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmFW")
 
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,6)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmSDW_kg")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,7)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmLength")
 
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,8)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmOutdi")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,9)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmIndi")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,10)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmnodes")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,11)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCmDW_g")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,12)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmBcirc_cm")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,13)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmFNMain")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,14)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmFNsmall")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,15)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmFN")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,16)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmLg")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,17)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmFD")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,18)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCCirc_cm")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,19)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmSRD")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,20)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmADD")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,21)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmGS")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,22)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
)

  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmflopc1")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,23)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  
  
####below all of these add the CV (Population ) in the model
####below all of these add the CV (Population ) in the model 
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVhday")
 
  myFarmCPU <- FarmCPU(
    Y=myY[,1:2],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVfday")
  
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,3)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVhhday")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,4)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVhfday")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,5)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVFW")

  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,6)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVSDW_kg")
 
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,7)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    CV=myQ[2:9],
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVLength")
 
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,8)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVOutdi")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,9)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVIndi")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,10)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVnodes")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,11)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVCmDW_g")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,12)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVBcirc_cm")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,13)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVFNMain")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,14)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVFNsmall")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,15)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVFN")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,16)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVLg")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,17)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVFD")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,18)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVCCirc_cm")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,19)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVSRD")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,20)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVADD")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,21)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVGS")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,22)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITandFarmCPU/FarmCPUnmCVFprin")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,23)],
    GD=myGD,
    GM=myGM,
    CV=myQ[2:9],
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  
  
  
  
  
  
  