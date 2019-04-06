####GAPIT 
###imported all the data need for GAPIT
myY <- read.csv("data/myYCulm.106.4220.csv")
myGD <- read.csv("data/myGDCulm.106.4220.csv")
myGM <- read.csv("data/myGMCulm.106.4220.csv")
myQ<- read.csv("data/myQCulm.106.4220.csv")

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
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/modelselection")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  PCA.total=3, 
  Model.selection = TRUE
)
###through the model selection, the optimal number of PCs is=0, expect the trait of FNMain, PCs=1, 
###for convienence, I set up the PC=0 for all of the traits

####Basic Scenario using PCs=0 without incoporating CV:Population Structure
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/SNPNOCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM
)

####Basic Scenario using PCs=0 with incoporate CV:Population Structure
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/SNPCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ
)

####Enhanced Compression without CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/ECMMLNOCV")
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
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
####becasue i saw the example does not use the PCA, but CV in the manual 
##########Using MLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/MLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
###########Using GLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="GLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
###########Using CMLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/CMLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
###########Using SUPPER SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/SUPPERSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )

  ###
  myY <- read.csv("data/myYCulm.115.3154.csv")
  myGD <- read.csv("data/myGDCulm.115.3154.csv")
  myGM <- read.csv("data/myGMCulm.115.3154.csv")
  myQ<- read.csv("data/myQCulm.115.3154.csv")
  
  ####Basic Scenario with Model selection in order to find the optimal number of PCs for different traits. 
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm115/modelselection")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    PCA.total=3, 
    Model.selection = TRUE
  )
  ###through the model selection, the optimal number of PCs is=0, expect the trait of FNMain, PCs=1, 
  ###for convienence, I set up the PC=0 for all of the traits
  
  ####Basic Scenario using PCs=0 without incoporating CV:Population Structure
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm115/SNPNOCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM
  )
  
  ####Basic Scenario using PCs=0 with incoporate CV:Population Structure
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm115/SNPCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ
  )
  
  ####Enhanced Compression without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm115/ECMMLNOCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm115/ECMMLCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm115/GLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm115/GLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
  ####becasue i saw the example does not use the PCA, but CV in the manual 
  ##########Using MLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm115/MLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
  ###########Using GLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm115/GLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="GLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
  ###########Using CMLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm115/CMLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
  ###########Using SUPPER SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm115/SUPPERSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
  
  ###
  myY <- read.csv("data/myYCulm.113.3489.csv")
  myGD <- read.csv("data/myGDCulm.113.3489.csv")
  myGM <- read.csv("data/myGMCulm.113.3489.csv")
  myQ<- read.csv("data/myQCulm.113.3489.csv")
  
  ####Basic Scenario with Model selection in order to find the optimal number of PCs for different traits. 
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm113/modelselection")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    PCA.total=3, 
    Model.selection = TRUE
  )
  ###through the model selection, the optimal number of PCs is=0, expect the trait of FNMain, PCs=1, 
  ###for convienence, I set up the PC=0 for all of the traits
  
  ####Basic Scenario using PCs=0 without incoporating CV:Population Structure
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm113/SNPNOCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM
  )
  
  ####Basic Scenario using PCs=0 with incoporate CV:Population Structure
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm113/SNPCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ
  )
  
  ####Enhanced Compression without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm113/ECMMLNOCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm113/ECMMLCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm113/GLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm113/GLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
  ####becasue i saw the example does not use the PCA, but CV in the manual 
  ##########Using MLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm113/MLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
  ###########Using GLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm113/GLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="GLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
  ###########Using CMLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm113/CMLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
  ###########Using SUPPER SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm113/SUPPERSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
  myY <- read.csv("data/myYCulm.123.2670.csv")
  myGD <- read.csv("data/myGDCulm.123.2670.csv")
  myGM <- read.csv("data/myGMCulm.123.2670.csv")
  myQ<- read.csv("data/myQCulm.123.2670.csv")
  
  ####Basic Scenario with Model selection in order to find the optimal number of PCs for different traits. 
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm123/modelselection")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    PCA.total=3, 
    Model.selection = TRUE
  )
  ###through the model selection, the optimal number of PCs is=0, expect the trait of FNMain, PCs=1, 
  ###for convienence, I set up the PC=0 for all of the traits
  
  ####Basic Scenario using PCs=0 without incoporating CV:Population Structure
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm123/SNPNOCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM
  )
  
  ####Basic Scenario using PCs=0 with incoporate CV:Population Structure
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm123/SNPCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ
  )
  
  ####Enhanced Compression without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm123/ECMMLNOCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm123/ECMMLCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm123/GLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm123/GLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
  ####becasue i saw the example does not use the PCA, but CV in the manual 
  ##########Using MLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm123/MLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
  ###########Using GLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm123/GLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="GLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
  ###########Using CMLM SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm123/CMLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  
  ###########Using SUPPER SUPPER with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm123/SUPPERSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )