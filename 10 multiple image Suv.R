
###Suvr
###this one is for Surv with single locus analysis 
###import trait Surv
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.Surv.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.Surv"
str(GLM)
MLM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.Surv.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.Surv"
CMLM <- read.csv("Result GAPIT1/ECMMLCV/GAPIT.CMLM.Surv.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.Surv"
MLMSUPER <- read.csv("Result GAPIT1/MLMSUPPER/GAPIT.SUPER.Surv.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.Surv"
gwasResultsqq <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,42)]
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
P.less.0.04 <- subset(plot.p.ajusted, plot.p.ajusted[9] < 0.04 | 
                        plot.p.ajusted[10] < 0.04 | plot.p.ajusted[11] < 0.04 |
                        plot.p.ajusted[12] < 0.04 | plot.p.ajusted[13] < 0.04)
str(P.less.0.04)
###write the result
write.csv(P.less.0.04, file="Allimages/Surv/Surv.csv")

###plot QQ plot inone image
setwd("Allimages/Surv")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=1.28837E-05,
       signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(4.04775E-06,1.28837E-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(4.04775E-06,1.28837E-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for Surv with multiple locus analysis 
###this one is for Surv with multiple locus analysis 
###import trait Surv
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.Surv.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.Surv"
str(FarmCPU)
MLMM <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.Surv.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.Surv"
intermediate <- read.csv("mrMLMM2/ResultSuvall1/1_intermediate result.csv")
str(intermediate)
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.Surv"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.Surv"
###pKWmEB
pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
pKWmEB <- pKWmEB[2:3]
colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.Surv"
subset(pKWmEB,is.na(pKWmEB$pKWmEB.Surv))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,pKWmEB), type = "left", by="SNP")
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
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,9] < 0.05 | 
                        plot.p.adjusted[,10] < 0.05)
str(P.less.0.05)
if (nrow(P.less.0.05) > 0){
  colnames(P.less.0.05)[colnames(P.less.0.05)=="FarmCPU.Surv"] <- "P.value"
  P.less.0.05$Method <- paste("FarmCPU")
  P.less.0.05$Trait.name <- paste("Surv")
  P.less.0.05 <- P.less.0.05[,c(12,11,1:4)]
  ###The result from mrMLM
  mrmlm.final <- read.csv("mrMLMM2/ResultSuvall1/1_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait1"] <- "Surv"
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
  write.csv(total, file="Allimages/Surv/mrmlm.final2.csv")
} else {
  mrmlm.final <- read.csv("mrMLMM2/ResultSuvall1/1_Final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait1"] <- "Surv"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.NAME.)
  str(mrmlm.final)
  mrmlm.final <- mrmlm.final[,c(2:6, 15)]
  write.csv(mrmlm.final, file="Allimages/Surv/mrmlm.final.csv")
}
setwd("Allimages/Surv")
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


