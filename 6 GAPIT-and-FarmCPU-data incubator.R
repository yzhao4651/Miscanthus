####GAPIT 
###imported all the data need for GAPIT
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

###########"GLM","MLM","MLMM","FarmCPU" without CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/GLMMLMMLMMFarmCPU")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
###########"GLM","MLM","MLMM","FarmCPU" with CV
  setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT/GLMMLMMLMMFarmCPUCV")
  myGAPIT <- GAPIT(
    Y=myY,
    GD=myGD,
    GM=myGM,
    CV=myQ, 
    model=c("GLM","MLM","MLMM","FarmCPU")
  )
  
