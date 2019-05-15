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
FarmCPU <- read.csv("FarmCPUresultall/FarmCPU.fprind.GWAS.Results.csv")
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
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait13"] <- "fprind"
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
plot.p.ajusted <- adj_P_function(plot,4,8)
P.less.0.05 <- subset(plot.p.ajusted, plot.p.ajusted[9] < 0.05 | 
                        plot.p.ajusted[10] < 0.05 | plot.p.ajusted[11] < 0.05 |
                        plot.p.ajusted[12] < 0.05 | plot.p.ajusted[13] < 0.05)
str(P.less.0.05)
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
plot.p.ajusted <- adj_P_function(plot,4,8)
P.less.0.05 <- subset(plot.p.ajusted, plot.p.ajusted[9] < 0.05 | 
                        plot.p.ajusted[10] < 0.05 | plot.p.ajusted[11] < 0.05 |
                        plot.p.ajusted[12] < 0.05 | plot.p.ajusted[13] < 0.05)
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
plot.p.ajusted <- adj_P_function(plot,4,8)
str(plot.p.ajusted)
P.less.0.05 <- subset(plot.p.ajusted, plot.p.ajusted[9] < 0.05 | 
                        plot.p.ajusted[10] < 0.05 | plot.p.ajusted[11] < 0.05 |
                        plot.p.ajusted[12] < 0.05 | plot.p.ajusted[13] < 0.05)
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

###combine all of the result from FarmCPU
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("FarmCPUresultall/FarmCPU.fprind.GWAS.Results.csv")
colnames(FarmCPU)[colnames(FarmCPU)=="Name"] <- "SNP"
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.fprind"
FarmCPU.d <- FarmCPU[1:4]
str(FarmCPU.d)

FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinW.GWAS.Results.csv")
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.fprinW"
FarmCPU.W <- FarmCPU[,c(1,4)]
str(FarmCPU.W)

FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinGW.GWAS.Results.csv")
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.fprinGW"
FarmCPU.GW <- FarmCPU[,c(1,4)]
str(FarmCPU.GW)

FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinM.GWAS.Results.csv")
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.fprinM"
FarmCPU.M <- FarmCPU[,c(1,4)]
str(FarmCPU.M)
plot <-  plyr::join_all(list(FarmCPU.d,FarmCPU.W,FarmCPU.GW,FarmCPU.M), type = "left", by="SNP")
str(plot)

setwd("Allimages/FarmCPU.all")
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

