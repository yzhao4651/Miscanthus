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

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYO.106.4834.csv")
str(myY)
myGD <- read.csv("data/myGDO.106.4834.csv")
str(myGD)
myGM <- read.csv("data/myGMO.106.4834.csv")
str(myGM)
myQ<- read.csv("data/myQO.106.4834.csv")
str(myQ)

####Enhanced Compression without CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAECMMLNOCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/106OWAECMMLNOCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAECMMLCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/106OWAECMMLCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAGLMMLMMLMMFarmCPU")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/106OWAGLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAGLMMLMMLMMFarmCPUCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/106OWAGLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
###CMLM, SUPPER
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/146CMLMSUPPER")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/106OWACMLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  ####GAPIT 
  ###imported all the data need for GAPIT
  myY <- read.csv("data/myYO.116.4232.csv")
  str(myY)
  myGD <- read.csv("data/myGDO.116.4232.csv")
  str(myGD)
  myGM <- read.csv("data/myGMO.116.4232.csv")
  str(myGM)
  myQ<- read.csv("data/myQO.116.4232.csv") 
  str(myQ)
  ####Enhanced Compression without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAECMMLNOCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/116OWAECMMLNOCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAECMMLCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/116OWAECMMLCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAGLMMLMMLMMFarmCPU")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/116OWAGLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAGLMMLMMLMMFarmCPUCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/116OWAGLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###CMLM, SUPPER
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/146CMLMSUPPER")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/116OWACMLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  

  ####GAPIT 
  ###imported all the data need for GAPIT
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
  myY <- read.csv("data/myYO.126.3659.csv")
  str(myY)
  myGD <- read.csv("data/myGDO.126.3659.csv")
  str(myGD)
  myGM <- read.csv("data/myGMO.135.3628.csv")
  str(myGM)
  myQ<- read.csv("data/myQO.126.3659.csv")
  str(myQ)
  ####Enhanced Compression without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAECMMLNOCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/126OWAECMMLNOCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAECMMLCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/126OWAECMMLCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAGLMMLMMLMMFarmCPU")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/126OWAGLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAGLMMLMMLMMFarmCPUCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/126OWAGLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###CMLM, SUPPER
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/146CMLMSUPPER")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/126OWACMLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )

  ####GAPIT 
  ###imported all the data need for GAPIT
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
  myY <- read.csv("data/myYO.135.2855.csv")
  str(myY)
  myGD <- read.csv("data/myGDO.135.2855.csv")
  str(myGD)
  myGM <- read.csv("data/myGMO.135.2855.csv")
  str(myGM)
  myQ<- read.csv("data/myQO.135.2855.csv") 
  str(myQ)
  ####Enhanced Compression without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAECMMLNOCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/135OWAECMMLNOCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAECMMLCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/135OWAECMMLCV")
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
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAGLMMLMMLMMFarmCPU")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/135OWAGLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/106OWAGLMMLMMLMMFarmCPUCV")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/135OWAGLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
  ###CMLM, SUPPER
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/146CMLMSUPPER")
  setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/135OWACMLMSUPPER")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ,
    sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
    sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
    LD=0.1
  )
  