####GAPIT 
###imported all the data need for GAPIT
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

####Basic Scenario with Model selection in order to find the optimal number of PCs for different traits. 
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYC.106.4202.csv")
myGD <- read.csv("data/myGDC.106.4202.csv")
myGM <- read.csv("data/myGMC.106.4202.csv")
myQ<- read.csv("data/myQC.106.4202.csv")
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/ECMMLNOCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C106ECMMLNOCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
  group.from=100,
  group.to=5000, 
  group.by=10
)
####Enhanced Compression with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/ECMMLCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C106ECMMLCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
    group.from=100,
    group.to=5000, 
    group.by=10
  )
###########"GLM","MLM","MLMM","FarmCPU" without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMMLMMLMMFarmCPU")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C106GLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMMLMMLMMFarmCPUCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C106GLMMLMMLMMFarmCPUCV")
   myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM")
  )
####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
####becasue i saw the example does not use the PCA, but CV in the manual 
##########Using MLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/MLMSUPPER")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C106MLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
###########Using CMLM SUPPER with CV
  ###
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
  myY <- read.csv("data/myYC.116.3293.csv")
  myGD <- read.csv("data/myGDC.116.3293.csv")
  myGM <- read.csv("data/myGMC.116.3293.csv")
  myQ<- read.csv("data/myQC.116.3293.csv")
  
  ####Basic Scenario with Model selection in order to find the optimal number of PCs for different traits. 
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/ECMMLNOCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C116ECMMLNOCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
    group.from=100,
    group.to=5000, 
    group.by=10
  )
  
  ####Enhanced Compression with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/ECMMLCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C116ECMMLCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
    group.from=100,
    group.to=5000, 
    group.by=10
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMMLMMLMMFarmCPU")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C116GLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMMLMMLMMFarmCPUCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C116GLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM")
  )
  
  ####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
  ####becasue i saw the example does not use the PCA, but CV in the manual 
  ##########Using MLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/MLMSUPPER")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C116MLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )

  ###
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
  myY <- read.csv("data/myYC.125.2562.csv")
  myGD <- read.csv("data/myGDC.125.2562.csv")
  myGM <- read.csv("data/myGMC.125.2562.csv")
  myQ<- read.csv("data/myQC.125.2562.csv")
  
  ####Basic Scenario with Model selection in order to find the optimal number of PCs for different traits. 
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/ECMMLNOCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C125ECMMLNOCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
    group.from=100,
    group.to=5000, 
    group.by=10
  )
  
  ####Enhanced Compression with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/ECMMLCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C125ECMMLCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
    group.from=100,
    group.to=5000, 
    group.by=10
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMMLMMLMMFarmCPU")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C125GLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMMLMMLMMFarmCPUCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C125GLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM")
  )
  
  ####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
  ####becasue i saw the example does not use the PCA, but CV in the manual 
  ##########Using MLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/MLMSUPPER")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C125MLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )

  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
  myY <- read.csv("data/myYf.106.4185.csv")
  myGD <- read.csv("data/myGDf.106.4185.csv")
  myGM <- read.csv("data/myGMf.106.4185.csv")
  myQ<- read.csv("data/myQf.106.4185.csv")
  ####Basic Scenario with Model selection in order to find the optimal number of PCs for different traits. 
  
  ####Enhanced Compression without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITflo106/ECMMLNOCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f106ECMMLNOCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
    group.from=100,
    group.to=5000, 
    group.by=10
  )
  
  ####Enhanced Compression with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITflo106/ECMMLCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f106ECMMLCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
    group.from=100,
    group.to=5000, 
    group.by=10
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" without CV
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPITflo106/GLMMLMMLMMFarmCPU")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f106GLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPITflo106/GLMMLMMLMMFarmCPUCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f106GLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM")
  )
  
  ####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
  ####becasue i saw the example does not use the PCA, but CV in the manual 
  ##########Using MLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITflo106/MLMSUPPER")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f106MLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
  ###
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
  myY <- read.csv("data/myYf.116.3077.csv")
  myGD <- read.csv("data/myGDf.116.3077.csv")
  myGM <- read.csv("data/myGMf.116.3077.csv")
  myQ<- read.csv("data/myQf.116.3077.csv")
  
  
  ####Enhanced Compression without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITflo106/ECMMLNOCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f116ECMMLNOCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
    group.from=100,
    group.to=5000, 
    group.by=10
  )
  
  ####Enhanced Compression with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITflo106/ECMMLCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f116ECMMLCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
    group.from=100,
    group.to=5000, 
    group.by=10
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" without CV
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPITflo106/GLMMLMMLMMFarmCPU")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f116GLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPITflo106/GLMMLMMLMMFarmCPUCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f116GLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM")
  )
  
  ####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
  ####becasue i saw the example does not use the PCA, but CV in the manual 
  ##########Using MLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITflo106/MLMSUPPER")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f116MLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
  myY <- read.csv("data/myYf.122.2272.csv")
  myGD <- read.csv("data/myGDf.122.2272.csv")
  myGM <- read.csv("data/myGMf.122.2272.csv")
  myQ<- read.csv("data/myQf.122.2272.csv")
  
  
  ####Enhanced Compression without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITflo106/ECMMLNOCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f122ECMMLNOCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
    group.from=100,
    group.to=5000, 
    group.by=10
  )
  
  ####Enhanced Compression with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITflo106/ECMMLCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f122ECMMLCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
    group.from=100,
    group.to=5000, 
    group.by=10
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" without CV
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPITflo106/GLMMLMMLMMFarmCPU")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f122GLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPITflo106/GLMMLMMLMMFarmCPUCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f122GLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM")
  )
  
  ####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
  ####becasue i saw the example does not use the PCA, but CV in the manual 
  ##########Using MLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITflo106/MLMSUPPER")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/f122MLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  