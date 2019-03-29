####GAPIT 
###imported all the data need for GAPIT
myY <- read.csv("data/myY.129.4033.csv")
myGD <- read.csv("data/myGD.129.4033.csv")
myGM <- read.csv("data/myGM.129.4033.csv")
myQ<- read.csv("data/myQ.129.4033.csv")

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
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/129modelselection")
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
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/129SNPNOCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM
)

####Basic Scenario using PCs=0 with incoporate CV:Population Structure
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/129SNPCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ
)

####Enhanced Compression without CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/129ECMMLNOCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/129ECMMLCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/129GLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/129GLMMLMMLMMFarmCPUCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/129MLMSUPPER")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/129GLMSUPPER")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/129CMLMSUPPER")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/129SUPPERSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )