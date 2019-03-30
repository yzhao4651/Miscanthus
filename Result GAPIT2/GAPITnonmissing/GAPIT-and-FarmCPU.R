####GAPIT and FarmCPU
####import the dataset need for software 
myY <- read.csv("~/Documents/whole traits/myYimputedSNP19.csv",na.strings = c("",".","NA"))
str(myY)
myGD <- read.csv("~/Documents/whole traits/myGDimputedSNP19.csv",na.strings = c("",".","NA"))
str(myGD)
myGDGAPIT <- data.frame(myGD,row.names=1)
str(myGDGAPIT)
myGM <- read.csv("~/Documents/whole traits/myGMimputedSNP19.csv",na.strings = c("",".","NA"))
str(myGM)
myQ<- read.csv("~/Documents/whole traits/myQimputedSNP19.csv"),na.strings = c("",".","NA")
str(myQ)

####GAPIT
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
  
####Basic Scenario of Compressed MLM with CV
  setwd("~/Documents/GAPITandFarmCPU/GAPITnonmissingcv")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGDGAPIT,
    GM=myGM,
    CV=myQ
  )
  
####Basic Scenario of Compressed MLM without CV 
  setwd("~/Documents/GAPITandFarmCPU/GAPITnonmissingnocv")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM
  )
  
  
####Basic Scenario of CMLM method without CV
  setwd("~/Documents/GAPITandFarmCPU/GAPITnonmissingsuper")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    sangwich.top="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER
  )
  
####Basic Scenario of SUPER method without CV
  setwd("~/Documents/GAPITandFarmCPU/GAPITCMLM")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER
    sangwich.bottom="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER
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
  str(myY)
  
  
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmFW")
  myFarmCPU <- FarmCPU(
    Y=myY[,1:2],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmLength")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,3)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmOutdi")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,4)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd(~/Documents/GAPITandFarmCPU//FarmCPUnmIndi")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,5)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd(~/Documents/GAPITandFarmCPU/FarmCPUnmnodes")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,6)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd(~/Documents/GAPITandFarmCPU/FarmCPUnmBC")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,7)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmFN")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,8)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmhday")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,9)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmfday")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,10)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmhhday")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,11)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmhfday")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,12)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmflopc1")
  myFarmCPU <- FarmCPU(
    Y=myY[, c(1,13)],
    GD=myGD,
    GM=myGM,
    maxLoop=10,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )
  
  
  str(myQ)
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmCVFW")
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
  
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmCVLength")
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
  
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmCVOutdi")
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
  
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmCVIndi")
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
  setwd("~/Documents/GAPITandFarmCPU/FarmCPUnmCVnodes")
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
  
  setwd("C:/Users/Admin/Desktop/Farmcpu19/FarmCPUnmCVBC")
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
  
  setwd("C:/Users/Admin/Desktop/Farmcpu19/FarmCPUnmCVFN")
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
  
  setwd("C:/Users/Admin/Desktop/Farmcpu19/FarmCPUnmCVhday")
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
  
  setwd("C:/Users/Admin/Desktop/Farmcpu19/FarmCPUnmCVfday")
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
  
  setwd("C:/Users/Admin/Desktop/Farmcpu19/FarmCPUnmCVhhday")
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
  
  setwd("C:/Users/Admin/Desktop/Farmcpu19/FarmCPUnmCVhfday")
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
  setwd("C:/Users/Admin/Desktop/Farmcpu19/FarmCPUnmCVflopc1")
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
  
  
  
  
  
  
  
  
  