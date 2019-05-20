### Multiple images for flowering related traits
###qeustion: 
###1: The path I set up i do not know if it is hard for you. 
###I do not know if there is a better way to work for us?
###2: I did trait one by one, do you have any idea that can do it one time for many traits?
###HD_1
###this one is for HD_1 with single locus analysis 
###import trait HD_1
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.HD_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.HD_1"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.HD_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.HD_1"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.HD_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.HD_1"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.HD_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.HD_1"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,60)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")

###using this function to adjust the p value 
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)

###and get all of the SNPS suitable for the p-value
p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.HD_1 < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.HD_1"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("HD_1")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.HD_1 < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.HD_1"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("HD_1")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.HD_1 < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.HD_1"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("HD_1")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.HD_1 < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.HD_1"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("HD_1")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.HD_1 < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.HD_1"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("HD_1")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/HD_1/HD_1.1.csv")
###plot QQ plot inone image
setwd("Allimages/HD_1")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=1.13e-5,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1.13e-5),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1.13e-5),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
       signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for HD_1 with multiple locus analysis 
###this one is for HD_1 with multiple locus analysis 
###import trait HD_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.HD_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.HD_1"
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.HD_1.GWAS.Results.csv")
MLMM <- read.csv("Result GAPIT1/MLMM/GAPIT.MLMM.HD_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.HD_1"
intermediate <- read.csv("Resultfloall1/4_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.HD_1"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.HD_1"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.HD_1"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.HD_1))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1
## get the final data sets with suitable p-value 
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.HD_1"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("HD_1")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/4_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait19"] <- "HD_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/HD_1/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/4_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait19"] <- "HD_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/HD_1/mrmlm.final.csv")
}
setwd("Allimages/HD_1")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###FD_1
###this one is for FD_1 with single locus analysis 
###import trait FD_1
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FD_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FD_1"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FD_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FD_1"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.FD_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FD_1"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FD_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FD_1"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,61)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
str(plot)
min(plot[4:8])
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.FD_1 < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.FD_1"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("FD_1")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.FD_1 < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.FD_1"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("FD_1")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.FD_1 < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.FD_1"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("FD_1")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.FD_1 < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.FD_1"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("FD_1")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.FD_1 < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.FD_1"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("FD_1")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/FD_1/FD_1.1.csv")


###plot QQ plot inone image
setwd("Allimages/FD_1")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=1.76E-06,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1.76E-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1.76E-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
       signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FD_1 with multiple locus analysis 
###this one is for FD_1 with multiple locus analysis 
###import trait FD_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FD_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FD_1"
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FD_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FD_1"
intermediate <- read.csv("Resultfloall1/5_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.FD_1"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.FD_1"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.FD_1"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.FD_1))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1
str(plot)
## get the final data set with suitabbl p-value:
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,7] < 0.05 | plot.p.adjusted[,8] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.FD_1"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("FD_1")
  P.less.0.05 <- P.less.0.05[,c(10,9,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/5_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait20"] <- "FD_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/FD_1/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/5_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait20"] <- "FD_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/FD_1/mrmlm.final.csv")
}
setwd("Allimages/FD_1")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###HD_50.
###this one is for HD_50. with single locus analysis 
###import traitHD_50.
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
#setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.HD_50..GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.HD_50."
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.HD_50..GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.HD_50."
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.HD_50..GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.HD_50."
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.HD_50..GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.HD_50."
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,62)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
str(plot)
min(plot[4:8])
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)


###plot QQ plot inone image
setwd("Allimages/HD_50")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=4.640832e-06,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(4.640832e-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(4.640832e-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
       signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for HD_50. with multiple locus analysis 
###this one is for HD_50. with multiple locus analysis 
###import trait HD_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.HD_50..GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.HD_50."
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.HD_50..GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.HD_50."
intermediate <- read.csv("Resultfloall1/6_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.HD_50."
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.HD_50."
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.HD_50."
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.HD_50.))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.HD_50."] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("HD_50")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/6_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait21"] <- "HD_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/HD_50/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/6_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait21"] <- "HD_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/HD_50/mrmlm.final.csv")
}

setwd("Allimages/HD_50")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###FD_50.
###this one is for FD_50. with single locus analysis 
###import trait FD_50.
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FD_50..GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FD_50."
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FD_50..GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FD_50."
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.FD_50..GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FD_50."
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FD_50..GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FD_50."
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,63)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
str(plot)
min(plot[4:8])
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.FD_50. < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.FD_50."] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("FD_50")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.FD_50. < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.FD_50."] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("FD_50")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.FD_50. < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.FD_50."] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("FD_50")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.FD_50. < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.FD_50."] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("FD_50")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.FD_50. < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.FD_50."] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("FD_50")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/FD_50./FD_50.1.csv")


###plot QQ plot inone image
setwd("Allimages/FD_50")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=8.85e-7,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(8.85e-7,9.83e-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(8.85e-7,9.83e-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FD_50. with multiple locus analysis 
###this one is for FD_50. with multiple locus analysis 
###import trait FD_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FD_50..GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FD_50."
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FD_50..GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FD_50."
intermediate <- read.csv("Resultfloall1/7_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.FD_50."
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.FD_50."
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.FD_50."
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.FD_50.))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.FD_50."] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("FD_50")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/7_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait22"] <- "FD_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/FD_50/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/7_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait22"] <- "FD_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/FD_50/mrmlm.final.csv")
}

setwd("Allimages/FD_50")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###HW_1
###this one is for HW_1 with single locus analysis 
###import trait HW_1
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.HW_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.HW_1"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.HW_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.HW_1"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.HW_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.HW_1"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.HW_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.HW_1"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,64)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.HW_1 < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.HW_1"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("HW_1")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.HW_1 < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.HW_1"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("HW_1")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.HW_1 < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.HW_1"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("HW_1")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.HW_1 < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.HW_1"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("HW_1")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.HW_1 < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.HW_1"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("HW_1")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/HW_1/HW_1.1.csv")

###plot QQ plot inone image
setwd("Allimages/HW_1")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=8E-06,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1.1e-6, 8.4E-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1.1e-6, 8.4E-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for HW_1 with multiple locus analysis 
###this one is for HW_1 with multiple locus analysis 
###import trait HW_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.HW_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.HW_1"
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.HW_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.HW_1"
intermediate <- read.csv("Resultfloall1/8_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.HW_1"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.HW_1"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.HW_1"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.HW_1))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.HW_1"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("HW_1")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/8_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait23"] <- "HW_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/HW_1/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/8_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait23"] <- "HW_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/HW_1/mrmlm.final.csv")
}

setwd("Allimages/HW_1")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###FW_1
###this one is for FW_1 with single locus analysis 
###import trait FW_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FW_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FW_1"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FW_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FW_1"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.FW_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FW_1"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FW_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FW_1"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,65)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.FW_1 < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.FW_1"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("FW_1")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.FW_1 < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.FW_1"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("FW_1")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.FW_1 < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.FW_1"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("FW_1")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.FW_1 < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.FW_1"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("FW_1")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.FW_1 < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.FW_1"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("FW_1")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/FW_1/FW_1.1.csv")

###plot QQ plot inone image
setwd("Allimages/FW_1")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=1.33E-05,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(2.96E-06,1.33E-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(2.96E-06,1.33E-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FW_1 with multiple locus analysis 
###this one is for FW_1 with multiple locus analysis 
###import trait FW_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FW_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FW_1"
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FW_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FW_1"
intermediate <- read.csv("Resultfloall1/9_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.FW_1"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.FW_1"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.FW_1"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.FW_1))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)
str(plot.p.adjusted)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.FW_1"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("FW_1")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/9_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait24"] <- "FW_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/FW_1/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/9_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait24"] <- "FW_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/FW_1/mrmlm.final.csv")
}

setwd("Allimages/FW_1")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###HW_50.
###this one is for HW_50. with single locus analysis 
###import trait HW_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.HW_50.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.HW_50."
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.HW_50.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.HW_50."
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.HW_50.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.HW_50."
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.HW_50.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.HW_50."
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,66)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)
str(plot.p.adjusted)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,9] < 0.05 | 
                        plot.p.adjusted[,10] < 0.05 | plot.p.adjusted[,11] < 0.05 |
                        plot.p.adjusted[,12] < 0.05 | plot.p.adjusted[,13] < 0.05)
str(P.less.0.05)
write.csv(P.less.0.05, file="Allimages/HW_50/HW_50.p.less.0.05.csv")


###plot QQ plot inone image
setwd("Allimages/HW_50")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=1e-6,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for HW_50. with multiple locus analysis 
###this one is for HW_50. with multiple locus analysis 
###import trait HW_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.HW_50.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.HW_50."
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.HW_50.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.HW_50."
intermediate <- read.csv("Resultfloall1/10_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.HW_50."
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.HW_50."
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.HW_50."
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.HW_50.))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.HW_50."] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("HW_50")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/10_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait25"] <- "HW_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/HW_50/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/10_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait25"] <- "HW_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/HW_50/mrmlm.final.csv")
}

setwd("Allimages/HW_50")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###FW_50.
###this one is for FW_50. with single locus analysis 
###import trait FW_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FW_50.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FW_50."
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FW_50.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FW_50."
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.FW_50.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FW_50."
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FW_50.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FW_50."
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,67)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.FW_50. < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.FW_50."] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("FW_50")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.FW_50. < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.FW_50."] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("FW_50")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.FW_50. < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.FW_50."] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("FW_50")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.FW_50. < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.FW_50."] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("FW_50")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.FW_50. < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.FW_50."] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("FW_50")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/FW_50/FW_50.1.csv")

setwd("Allimages/FW_50")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=3.7E-06,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(3.7E-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(3.7E-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FW_50. with multiple locus analysis 
###this one is for FW_50. with multiple locus analysis 
###import trait FW_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FW_50.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FW_50."
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FW_50.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FW_50."
intermediate <- read.csv("Resultfloall1/11_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.FW_50."
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.FW_50."
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.FW_50."
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.FW_50.))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.FW_50."] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("FW_50")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/11_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait26"] <- "FW_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/FW_50/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/11_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait26"] <- "FW_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/FW_50/mrmlm.final.csv")
}

setwd("Allimages/FW_50")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###GHW_1
###this one is for GHW_1 with single locus analysis 
###import trait GHW_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.GHW_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.GHW_1"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.GHW_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.GHW_1"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.GHW_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.GHW_1"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.GHW_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.GHW_1"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,68)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.GHW_1 < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.GHW_1"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("GHW_1")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.GHW_1 < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.GHW_1"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("GHW_1")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.GHW_1 < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.GHW_1"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("GHW_1")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.GHW_1 < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.GHW_1"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("GHW_1")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.GHW_1 < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.GHW_1"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("GHW_1")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/GHW_1/GHW_1.1.csv")

setwd("Allimages/GHW_1")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=4.57E-06,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(4.57E-06,1.43E-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(4.57E-06,1.43E-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for GHW_1 with multiple locus analysis 
###this one is for GHW_1 with multiple locus analysis 
###import trait GHW_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.GHW_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.GHW_1"
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.GHW_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.GHW_1"
intermediate <- read.csv("Resultfloall1/12_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.GHW_1"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.GHW_1"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.GHW_1"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.GHW))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.GHW_1"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("GHW_1")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/12_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait27"] <- "GHW_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/GHW_1/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/12_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait27"] <- "GHW_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/GHW_1/mrmlm.final.csv")
}

setwd("Allimages/GHW_1")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###GFW_1
###this one is for GFW_1 with single locus analysis 
###import trait GFW_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.GFW_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.GFW_1"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.GFW_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.GFW_1"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.GFW_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.GFW_1"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.GFW_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.GFW_1"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,69)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)
str(plot.p.adjusted)


p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.GFW_1 < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.GFW_1"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("GFW_1")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.GFW_1 < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.GFW_1"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("GFW_1")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.GFW_1 < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.GFW_1"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("GFW_1")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.GFW_1 < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.GFW_1"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("GFW_1")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.GFW_1 < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.GFW_1"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("GFW_1")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/GFW_1/GFW_1.1.csv")

###plot QQ plot inone image
setwd("Allimages/GFW_1")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")

CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=8.73E-06,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1.2e-6,8.73E-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1.2e-6,8.73E-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for GFW_1 with multiple locus analysis 
###this one is for GFW_1 with multiple locus analysis 
###import trait GFW_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.GFW_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.GFW_1"
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.GFW_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.GFW_1"
intermediate <- read.csv("Resultfloall1/13_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.GFW_1"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.GFW_1"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.GFW_1"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.GFW_1))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.GFW_1"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("GFW_1")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/13_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait28"] <- "GFW_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/GFW_1/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/13_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait28"] <- "GFW_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/GFW_1/mrmlm.final.csv")
}

setwd("Allimages/GFW_1")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###GHW_50.
###this one is for GHW_50. with single locus analysis 
###import trait GHW_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.GHW_50.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.GHW_50."
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.GHW_50.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.GHW_50."
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.GHW_50.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.GHW_50."
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.GHW_50.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.GHW_50."
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,70)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)
str(plot.p.adjusted)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,9] < 0.05 | 
                        plot.p.adjusted[,10] < 0.05 | plot.p.adjusted[,11] < 0.05 |
                        plot.p.adjusted[,12] < 0.05 | plot.p.adjusted[,13] < 0.05)
str(P.less.0.05)
write.csv(P.less.0.05, file="Allimages/GHW_50/GHW_50.p.less.0.05.csv")
min(plot[4:8])

###plot QQ plot inone image
setwd("Allimages/GHW_50")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")

CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=4e-06,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(4e-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(4e-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for GHW_50. with multiple locus analysis 
###this one is for GHW_50. with multiple locus analysis 
###import trait GHW_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.GHW_50.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.GHW_50."
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.GHW_50.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.GHW_50."
intermediate <- read.csv("Resultfloall1/14_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.GHW_50."
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.GHW_50."
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.GHW_50."
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.GHW_50.))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,7] < 0.05 | 
                        plot.p.adjusted[,8] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.GHW_50."] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("GHW_50")
  P.less.0.05 <- P.less.0.05[,c(10,9,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/14_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait29"] <- "GHW_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/GHW_50/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/14_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait29"] <- "GHW_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/GHW_50/mrmlm.final.csv")
}

setwd("Allimages/GHW_50")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###GFW_50.
###this one is for GFW_50. with single locus analysis 
###import trait GFW_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.GFW_50.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.GFW_50."
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.GFW_50.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.GFW_50."
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.GFW_50.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.GFW_50."
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.GFW_50.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.GFW_50."
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,71)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.GFW_50. < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.GFW_50."] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("GFW_50")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.GFW_50. < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.GFW_50."] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("GFW_50")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.GFW_50. < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.GFW_50."] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("GFW_50")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.GFW_50. < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.GFW_50."] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("GFW_50")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.GFW_50. < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.GFW_50."] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("GFW_50")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/GFW_50./GFW_50.1.csv")



###plot QQ plot inone image
setwd("Allimages/GFW_50")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")

###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=2.42E-06,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(8.82E-07,2.42E-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(8.82E-07,2.42E-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for GFW_50. with multiple locus analysis 
###this one is for GFW_50. with multiple locus analysis 
###import trait GFW_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.GFW_50.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.GFW_50."
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.GFW_50.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.GFW_50."
intermediate <- read.csv("Resultfloall1/15_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.GFW_50."
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.GFW_50."
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.GFW_50."
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.GFW_50.))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)
str(plot.p.adjusted)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.GFW_50."] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("GFW_50")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/14_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait29"] <- "GFW_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/GFW_50/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/14_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait29"] <- "GFW_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/GFW_50/mrmlm.final.csv")
}

setwd("Allimages/GFW_50")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###HM_1
###this one is for HM_1 with single locus analysis 
###import trait HM_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.HM_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.HM_1"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.HM_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.HM_1"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.HM_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.HM_1"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.HM_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.HM_1"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,72)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.HM_1 < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.HM_1"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("HM_1")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.HM_1 < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.HM_1"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("HM_1")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.HM_1 < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.HM_1"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("HM_1")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.HM_1 < 0.03),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.HM_1"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("HM_1")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.HM_1 < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.HM_1"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("HM_1")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/HM_1/HM_1.1.csv")



###plot QQ plot inone image
setwd("Allimages/HM_1")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=4.26E-06,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(4.26E-06,1.24E-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(4.26E-06,1.24E-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for HM_1 with multiple locus analysis 
###this one is for HM_1 with multiple locus analysis 
###import trait HM_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.HM_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.HM_1"
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.HM_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.HM_1"
intermediate <- read.csv("Resultfloall1/16_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.HM_1"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.HM_1"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.HM_1"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.HM_1))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1
str(plot)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,7] < 0.05 | 
                        plot.p.adjusted[,8] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.HM_1"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("HM_1")
  P.less.0.05 <- P.less.0.05[,c(10,9,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/15_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait30"] <- "HM_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/HM_1/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/15_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait30"] <- "HM_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/HM_1/mrmlm.final.csv")
}
str(plot)
setwd("Allimages/HM_1")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)



###FM_1
###this one is for FM_1 with single locus analysis 
###import trait FM_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FM_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FM_1"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FM_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FM_1"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.FM_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FM_1"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FM_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FM_1"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,73)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)
str(plot.p.adjusted)

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.FM_1 < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.FM_1"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("FM_1")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.FM_1 < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.FM_1"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("FM_1")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.FM_1 < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.FM_1"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("FM_1")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.FM_1 < 0.03),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.FM_1"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("FM_1")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.FM_1 < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.FM_1"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("FM_1")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/FM_1/FM_1.1.csv")


###plot QQ plot inone image
setwd("Allimages/FM_1")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=4.46E-06,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(4.46E-06,1.37E-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(4.46E-06,1.37E-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FM_1 with multiple locus analysis 
###this one is for FM_1 with multiple locus analysis 
###import trait FM_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FM_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FM_1"
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FM_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FM_1"
intermediate <- read.csv("Resultfloall1/17_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.FM_1"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.FM_1"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.FM_1"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.FM_1))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.FM_1"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("FM_1")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/16_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait31"] <- "FM_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/FM_1/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/16_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait31"] <- "FM_1"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/FM_1/mrmlm.final.csv")
}

setwd("Allimages/FM_1")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###HM_50.
###this one is for HM_50. with single locus analysis 
###import traitHM_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.HM_50.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.HM_50."
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.HM_50.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.HM_50."
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.HM_50.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.HM_50."
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.HM_50.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.HM_50."
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,74)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)
str(plot.p.adjusted)

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.HM_50. < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.HM_50."] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("HM_50")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.HM_50. < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.HM_50."] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("HM_50")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.HM_50. < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.HM_50."] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("HM_50")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.HM_50. < 0.03),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.HM_50."] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("HM_50")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.HM_50. < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.HM_50."] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("HM_50")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/HM_50/HM_50.1.csv")

###plot QQ plot inone image
setwd("Allimages/HM_50")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=3.5e-07,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(3.5e-07),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(3.5e-07),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for HM_50. with multiple locus analysis 
###this one is for HM_50. with multiple locus analysis 
###import trait HM_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.HM_50.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.HM_50."
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.HM_50.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.HM_50."
intermediate <- read.csv("Resultfloall1/18_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.HM_50."
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.HM_50."
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.HM_50."
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.HM_50.))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.HM_50."] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("HM_50")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/17_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait32"] <- "HM_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/HM_50/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/17_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait32"] <- "HM_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/HM_50/mrmlm.final.csv")
}

setwd("Allimages/HM_50")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###FM_50.
###this one is for FM_50. with single locus analysis 
###import trait FM_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FM_50.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FM_50."
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FM_50.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FM_50."
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.FM_50.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FM_50."
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FM_50.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FM_50."
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,75)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)
str(plot.p.adjusted)

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.FM_50. < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.FM_50."] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("FM_50")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.FM_50. < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.FM_50."] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("FM_50")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.FM_50. < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.FM_50."] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("FM_50")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.FM_50. < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.FM_50."] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("FM_50")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.FM_50. < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.FM_50."] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("FM_50")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/FM_50/FM_50.1.csv")

###plot QQ plot inone image
setwd("Allimages/FM_50")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=2.91E-05,
       signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(7.57E-06,2.91E-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(7.57E-06,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FM_50. with multiple locus analysis 
###this one is for FM_50. with multiple locus analysis 
###import trait FM_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FM_50.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FM_50."
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FM_50.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FM_50."
intermediate <- read.csv("Resultfloall1/19_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.FM_50."
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.FM_50."
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.FM_50."
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.FM_50.))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.FM_50."] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("FM_50")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/18_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait33"] <- "FM_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/FM_50/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/18_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait33"] <- "FM_50"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/FM_50/mrmlm.final.csv")
}

setwd("Allimages/FM_50")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###fprind
###this one is for fprind with single locus analysis 
###import trait fprind
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.fprind.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.fprind"

MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.fprind.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.fprind"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.fprind.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.fprind"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.fprind.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.fprind"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,76)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.ajusted <- adj_P_function(plot,4,8)
P.less.0.05 <- subset(plot.p.ajusted, plot.p.ajusted[9] < 0.05 | 
                        plot.p.ajusted[10] < 0.05 | plot.p.ajusted[11] < 0.05 |
                        plot.p.ajusted[12] < 0.05 | plot.p.ajusted[13] < 0.05)
str(P.less.0.05)
###write the result
write.csv(P.less.0.05, file="Allimages/fprind/fprind.csv")

###plot QQ plot inone image
setwd("Allimages/fprind")
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for fprind with multiple locus analysis 
###this one is for fprind with multiple locus analysis 
###import trait fprind
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprind.GWAS.Results.csv")
colnames(FarmCPU)[colnames(FarmCPU)=="Name"] <- "SNP"
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.fprind"
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.fprind.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.fprind"
intermediate <- read.csv("Resultfloall1/20_intermediate result.csv")
levels(intermediate$Trait.name)
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.fprind"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.fprind"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.fprind"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.fprind))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1
str(plot)
## get the final data set with suitabbl p-value:
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) == 0){
  mrmlm.final <- read.csv("Resultfloall1/20_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait35"] <- "fprind"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/fprind/mrmlm.final.csv")
} else {
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.fprind"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("fprind")
  P.less.0.05 <- P.less.0.05[,c(11, 10, 1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/20_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait35"] <- "fprind"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/fprind/mrmlm.final2.csv")
}
setwd("Allimages/fprindCV")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###fprinW
###this one is for fprinW with single locus analysis 
###import trait fprinW
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.fprinW.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.fprinW"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.fprinW.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.fprinW"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.fprinW.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.fprinW"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.fprinW.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.fprinW"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,77)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot,4,8)

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.fprinW < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.fprinW"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("fprinW")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.fprinW < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.fprinW"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("fprinW")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.fprinW < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.fprinW"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("fprinW")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.fprinW < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.fprinW"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("fprinW")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.fprinW < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.fprinW"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("fprinW")
str(p.rrBLUP)
library(lessR)
total <- Merge(p.SUPER, p.rrBLUP)
#install.packages("tidyverse")
library(tidyverse)
total <- tibble::rowid_to_column(p.SUPER, "ID")
total <- total[,c(1,6:7,2:5)]
str(total)
###write the result
write.csv(total, file="Allimages/fprinW/fprinW.1.csv")


###write the result
write.csv(P.less.0.05, file="Allimages/fprinW/fprinW.csv")
###plot QQ plot inone image
setwd("Allimages/fprinW")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=4.1e-06,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,4.1e-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,4.1e-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for fprinW with multiple locus analysis 
###this one is for fprinW with multiple locus analysis 
###import trait fprinW
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinW.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.fprinW"
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.fprinW.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.fprinW"
intermediate <- read.csv("Resultfloall1/21_intermediate result.csv")
levels(intermediate$Trait.name) 
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.fprinW"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.fprinW"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.fprinW"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.fprinW))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1
str(plot)
## get the final data set with suitabbl p-value:
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) == 0){
  mrmlm.final <- read.csv("Resultfloall1/21_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == levels(intermediate$Trait.name)] <- "fprinW"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/fprinW/mrmlm.final.csv")
} else {
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.fprinW"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("fprinW")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/21_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait36"] <- "fprinW"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/fprinW/mrmlm.final2.csv")
}
setwd("Allimages/fprinW")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
#library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
       signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###fprinGW
###this one is for fprinGW with single locus analysis 
###import trait fprinGW
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.fprinGW.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.fprinGW"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.fprinGW.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.fprinGW"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.fprinGW.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.fprinGW"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.fprinGW.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.fprinGW"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,78)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot,4,8)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[9] < 0.05 | 
                        plot.p.adjusted[10] < 0.05 | plot.p.adjusted[11] < 0.05 |
                        plot.p.adjusted[12] < 0.05 | plot.p.adjusted[13] < 0.05)
str(P.less.0.05)
###write the result
write.csv(P.less.0.05, file="Allimages/fprinGW/fprinGW.csv")
###plot QQ plot inone image
setwd("Allimages/fprinGW")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=min(plot[4:8]-0.1),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(min(plot[4:8])),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(min(plot[4:8])),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
       signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for fprinGW with multiple locus analysis 
###this one is for fprinGW with multiple locus analysis 
###import trait fprinGW
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinGW.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.fprinGW"
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.fprinGW.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.fprinGW"
intermediate <- read.csv("Resultfloall1/22_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.fprinGW"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.fprinGW"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.fprinGW"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.fprinGW))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

## get the final data set with suitabbl p-value:
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) == 0){
  mrmlm.final <- read.csv("Resultfloall1/22_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == levels(intermediate$Trait.name)] <- "fprinGW"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/fprinGW/mrmlm.final.csv")
} else {
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.fprinGW"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("fprinGW")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/22_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait37"] <- "fprinGW"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/fprinGW/mrmlm.final2.csv")
}

setwd("Allimages/fprinGW")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
       signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###fprinM
###this one is for fprinM with single locus analysis 
###import trait fprinM
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.fprinM.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.fprinM"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.fprinM.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.fprinM"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.fprinM.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.fprinM"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.fprinM.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.fprinM"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,79)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
min(plot[4:8])
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot,4,8)
str(plot.p.adjusted)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[9] < 0.05 | 
                        plot.p.adjusted[10] < 0.05 | plot.p.adjusted[11] < 0.05 |
                        plot.p.adjusted[12] < 0.05 | plot.p.adjusted[13] < 0.05)
str(P.less.0.05)
###write the result
write.csv(P.less.0.05, file="Allimages/fprinM/fprinM.csv")

###plot QQ plot inone image
setwd("Allimages/fprinM")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
###plot QQ plot inone image
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=2.5e-6,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(2.5e-6),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(2.5e-6),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
       signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for fprinM with multiple locus analysis 
###this one is for fprinM with multiple locus analysis 
###import trait fprinM
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinM.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.fprinM"
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.fprinM.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.fprinM"
intermediate <- read.csv("Resultfloall1/23_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.fprinM"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.fprinM"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.fprinM"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.fprinM))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

## get the final data set with suitabbl p-value:
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) == 0){
  mrmlm.final <- read.csv("Resultfloall1/23_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == levels(intermediate$Trait.name)] <- "fprinM"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/fprinM/mrmlm.final.csv")
} else {
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.fprinM"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("fprinM")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/23_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait38"] <- "fprinM"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/fprinM/mrmlm.final2.csv")
}

setwd("Allimages/fprinM")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
       signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###TFN.
###this one is for TFN. with single locus analysis 
###import trait TFN.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.TFN..GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.TFN."
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.TFN..GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.TFN."
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.TFN..GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.TFN."
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.TFN..GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.TFN."
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,54)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
##get the 
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.ajusted <- adj_P_function(plot,4,8)
P.less.0.05 <- subset(plot.p.ajusted, plot.p.ajusted[9] < 0.05 | 
                        plot.p.ajusted[10] < 0.05 | plot.p.ajusted[11] < 0.05 |
                        plot.p.ajusted[12] < 0.05 | plot.p.ajusted[13] < 0.05)
str(P.less.0.05)
###write the result
write.csv(P.less.0.05, file="Allimages/TFN/TFN.csv")

###plot QQ plot inone image
setwd("Allimages/TFN")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for TFN. with multiple locus analysis 
###this one is for TFN. with multiple locus analysis 
###import trait TFN.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.TFN..GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.TFN."
### if this one does not work get the data set, the second path will works
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.TFN..GWAS.Results.csv")
MLMM <- read.csv("Result GAPIT1/MLMM/GAPIT.MLMM.TFN..GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.TFN."
intermediate <- read.csv("Resultfloall1/1_intermediate result.csv")
str(intermediate)
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.TFN."
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.TFN."
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.TFN."
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.TFN.))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

str(plot)
## get the final data set with suitabbl p-value:
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.TFN."] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("TFN")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/1_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait13"] <- "TFN"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/TFN/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/1_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait13"] <- "TFN"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/TFN/mrmlm.final.csv")
}
setwd("Allimages/TFN")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###FNMain
###this one is for FNMain with single locus analysis 
###import trait FNMain
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FNMain.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FNMain"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FNMain.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FNMain"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.FNMain.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FNMain"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FNMain.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FNMain"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,55)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.ajusted <- adj_P_function(plot,4,8)
P.less.0.05 <- subset(plot.p.ajusted, plot.p.ajusted[9] < 0.05 | 
                        plot.p.ajusted[10] < 0.05 | plot.p.ajusted[11] < 0.05 |
                        plot.p.ajusted[12] < 0.05 | plot.p.ajusted[13] < 0.05)
str(P.less.0.05)
###write the result
write.csv(P.less.0.05, file="Allimages/FNMain/FNMain.csv")
min(plot[4:8])
###plot QQ plot inone image
setwd("Allimages/FNMain")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=1e-5,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-5),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-5),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FNMain with multiple locus analysis 
###this one is for FNMain with multiple locus analysis 
###import trait FNMain
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FNMain.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FNMain"
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FNMain.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FNMain"
intermediate <- read.csv("Resultfloall1/2_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.FNMain"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.FNMain"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.FNMain"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.FNMain))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1
## get the final data set with suitabbl p-value:
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.FNMain"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("FNMain")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/2_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait14"] <- "FNMain"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/FNMain/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/2_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait14"] <- "FNMain"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/FNMain/mrmlm.final.csv")
}
setwd("Allimages/FNMain")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###FNsmall
###this one is for FNsmall with single locus analysis 
###import trait FNsmall
##this one for desktop
#setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
##get the current path way
#getwd()
##set up the 
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FNsmall.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FNsmall"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FNsmall.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FNsmall"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.FNsmall.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FNsmall"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FNsmall.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FNsmall"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,56)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
str(plot)
min(plot[4:8])
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 8)
str(plot.p.adjusted)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,9] < 0.05 | 
                        plot.p.adjusted[,10] < 0.05 | plot.p.adjusted[,11] < 0.05 |
                        plot.p.adjusted[,12] < 0.05 | plot.p.adjusted[,13] < 0.05)
str(P.less.0.05)
write.csv(P.less.0.05, file="Allimages/FNsmall/FNsmall.p.less.0.05.csv")

###plot QQ plot inone image
setwd("Allimages/FNsmall")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=1.4e-05,
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1.4e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1.4e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
       signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FNsmall with multiple locus analysis 
###this one is for FNsmall with multiple locus analysis 
###import trait FNsmall
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FNsmall.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FNsmall"
MLMM <- read.csv("Result GAPIT1/MLMM/GAPIT.MLMM.FNsmall.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FNsmall"
intermediate <- read.csv("Resultfloall1/3_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.FNsmall"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.FNsmall"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.FNsmall"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.FNsmall))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1
## get the final data set with suitabbl p-value:
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

str(plot.p.adjusted)
min(plot.p.adjusted[8:9])
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.FNsmall"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("FNsmall")
  P.less.0.05 <- P.less.0.05[,c(11,10,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("Resultfloall1/3_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait15"] <- "FNsmall"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  ##combine the result from above two
  ## Vertical
  #install.packages("lessR")
  library(lessR)
  total <- Merge(mrmlm.final, P.less.0.05)
  #install.packages("tidyverse")
  library(tidyverse)
  total <- tibble::rowid_to_column(total, "ID")
  ###write the result
  write.csv(total, file="Allimages/FNsmall/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("Resultfloall1/3_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait15"] <- "FNsmall"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/FNsmall/mrmlm.final.csv")
}
setwd("Allimages/FNsmall")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.05/nrow(plot)),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.05/nrow(plot)),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)