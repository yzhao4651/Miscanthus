###CmDW_g
###this one is for CmDW_g with single locus analysis 
###import trait CmDW_g
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.CmDW_g.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.CmDW_g"
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.CmDW_g.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.CmDW_g"
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.CmDW_g.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.CmDW_g"
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.CmDW_g.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.CmDW_g"
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,18)]
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
p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.CmDW_g < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.CmDW_g"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("CmDW_g")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.CmDW_g < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.CmDW_g"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("CmDW_g")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.CmDW_g < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.CmDW_g"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("CmDW_g")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.CmDW_g < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.CmDW_g"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("CmDW_g")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.CmDW_g < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.CmDW_g"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("CmDW_g")
str(p.rrBLUP)
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result
if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/CmDW_g_1.csv",row.names = F)
}

###plot QQ plot inone image
setwd("AllimagesC125/CmDW_g")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.0000006,
       signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(0.0000006),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(0.0000006),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
       signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for CmDW_g with multiple locus analysis 
###this one is for CmDW_g with multiple locus analysis 
###import trait CmDW_g
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.CmDW_g.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.CmDW_g"
str(FarmCPU)
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.CmDW_g.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.CmDW_g"
intermediate <- read.csv("mrMLMM2/C125/1_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.CmDW_g"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.CmDW_g"

###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.CmDW_g"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.CmDW_g"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.CmDW_g))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA), type = "left", by="SNP")
plot[is.na(plot)] <- 1
## get the final data set with suitabbl p-value:
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,9] < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.CmDW_g"] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("CmDW_g")
  FarmCPU <-  FarmCPU[,c(6,5,1:4)]
  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted[,10] < 0.05),]
  colnames(MLMM)[colnames(MLMM)=="MLMM.CmDW_g"] <- "P.value"
  MLMM$Method <- paste("MLMM")   
  MLMM$Trait.name <- paste("CmDW_g")
  MLMM <- MLMM[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/1_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait1"] <- "CmDW_g"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
  total$SNP.N <- paste("125")
  total <- unique(total[,c(7,1:6)])
  write.csv(total, file="Culmallresults125/CmDW_g_2.csv",row.names = F)
    
  setwd("AllimagesC125/CmDW_g")
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
###Cml_cm
###this one is for Cml_cm with single locus analysis 
###import trait Cml_cm
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.Cml_cm.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.Cml_cm"
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.Cml_cm.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.Cml_cm"
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.Cml_cm.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.Cml_cm"
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.Cml_cm.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.Cml_cm"
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,19)]
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
p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.Cml_cm < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.Cml_cm"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("Cml_cm")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.Cml_cm < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.Cml_cm"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("Cml_cm")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.Cml_cm < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.Cml_cm"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("Cml_cm")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.Cml_cm < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.Cml_cm"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("Cml_cm")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.Cml_cm < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.Cml_cm"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("Cml_cm")
str(p.rrBLUP)
#install.packages("tidyverse")
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result
if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/Cml_cm_1.csv",row.names = F)
}

###plot QQ plot inone image
setwd("AllimagesC125/Cml_cm")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=6.6E-06,
       signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(6.6E-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(6.6E-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
       signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)



###this one is for Cml_cm with multiple locus analysis 
###this one is for Cml_cm with multiple locus analysis 
###import trait Cml_cm
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.Cml_cm.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.Cml_cm"
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.Cml_cm.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.Cml_cm"
intermediate <- read.csv("mrMLMM2/C125/2_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.Cml_cm"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.Cml_cm"

###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.Cml_cm"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.Cml_cm"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.Cml_cm))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA), type = "left", by="SNP")
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

  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,9] < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.Cml_cm"] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("Cml_cm")
  FarmCPU <-  FarmCPU[,c(6,5,1:4)]
  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted[,10] < 0.05),]   
  colnames(MLMM)[colnames(MLMM)=="MLMM.FNsmall"] <- "P.value" 
  MLMM$Method <- paste("MLMM")   MLMM$Trait.name <- paste("FNsmall")
  MLMM <- MLMM[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/2_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait2"] <- "Cml_cm"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]

    total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
    total$SNP.N <- paste("125")
    total <- unique(total[,c(7,1:6)])
    write.csv(total, file="Culmallresults125/Cml_cm_2.csv",row.names = F)
    setwd("AllimagesC125/Cml_cm")
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
    
    


###CmD_BI_mm
###this one is for CmD_BI_mm with single locus analysis 
###import trait CmD_BI_mm
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.CmD_BI_mm.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.CmD_BI_mm"
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.CmD_BI_mm.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.CmD_BI_mm"
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.CmD_BI_mm.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.CmD_BI_mm"
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.CmD_BI_mm.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.CmD_BI_mm"
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,20)]
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

##using this one
p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.CmD_BI_mm < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.CmD_BI_mm"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("CmD_BI_mm")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.CmD_BI_mm < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.CmD_BI_mm"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("CmD_BI_mm")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.CmD_BI_mm < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.CmD_BI_mm"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("CmD_BI_mm")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.CmD_BI_mm < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.CmD_BI_mm"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("CmD_BI_mm")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.CmD_BI_mm < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.CmD_BI_mm"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("CmD_BI_mm")
str(p.rrBLUP)
#install.packages("tidyverse")
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result
if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/CmD_BI_mm_1.csv",row.names = F)
}

setwd("AllimagesC125/CmD_BI_mm")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=1.21906E-05,
       signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(7.29759E-06,2.54697E-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(7.29759E-06,2.54697E-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###this one is for CmD_BI_mm with multiple locus analysis 
###this one is for CmD_BI_mm with multiple locus analysis 
###import trait CmD_BI_mm
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.CmD_BI_mm.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.CmD_BI_mm"
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.CmD_BI_mm.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.CmD_BI_mm"
intermediate <- read.csv("mrMLMM2/C125/3_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.CmD_BI_mm"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.CmD_BI_mm"
###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.CmD_BI_mm"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.CmD_BI_mm"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.CmD_BI_mm))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA), type = "left", by="SNP")
plot[is.na(plot)] <- 1
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,9] < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.CmD_BI_mm"] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("CmD_BI_mm")
  FarmCPU <-  FarmCPU[,c(6,5,1:4)]
  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted[,10] < 0.05),]   
  colnames(MLMM)[colnames(MLMM)=="MLMM.CmD_BI_mm"] <- "P.value"  
  MLMM$Method <- paste("MLMM")   MLMM$Trait.name <- paste("CmD_BI_mm")   
  MLMM <- MLMM[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/3_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait3"] <- "CmD_BI_mm"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]

    total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
    total$SNP.N <- paste("125")
    total <- unique(total[,c(7,1:6)])
    write.csv(total, file="Culmallresults125/CmD_BI_mm_2.csv",row.names = F)

    setwd("AllimagesC125/CmD_BI_mm")
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
    

###CmD_LI_mm
###this one is for CmD_LI_mm with single locus analysis 
###import trait CmD_LI_mm
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.CmD_LI_mm.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.CmD_LI_mm"
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.CmD_LI_mm.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.CmD_LI_mm"
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.CmD_LI_mm.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.CmD_LI_mm"
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.CmD_LI_mm.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.CmD_LI_mm"
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,21)]
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

##using this one
p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.CmD_LI_mm < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.CmD_LI_mm"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("CmD_LI_mm")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.CmD_LI_mm < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.CmD_LI_mm"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("CmD_LI_mm")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.CmD_LI_mm < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.CmD_LI_mm"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("CmD_LI_mm")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.CmD_LI_mm < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.CmD_LI_mm"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("CmD_LI_mm")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.CmD_LI_mm < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.CmD_LI_mm"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("CmD_LI_mm")
str(p.rrBLUP)
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result
if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/CmD_LI_mm_1.csv",row.names = F)
}

###plot QQ plot inone image
setwd("AllimagesC125/CmD_LI_mm")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=3e-6,
       signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(3e-6),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(3e-6),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
       signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###this one is for CmD_LI_mm with multiple locus analysis 
###this one is for CmD_LI_mm with multiple locus analysis 
###import trait CmD_LI_mm
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.CmD_LI_mm.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.CmD_LI_mm"
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.CmD_LI_mm.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.CmD_LI_mm"
intermediate <- read.csv("mrMLMM2/C125/4_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.CmD_LI_mm"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.CmD_LI_mm"
###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.CmD_LI_mm"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.CmD_LI_mm"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.CmD_LI_mm))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA), type = "left", by="SNP")
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
  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,9] < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.CmD_LI_mm"] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("CmD_LI_mm")
  FarmCPU <-  FarmCPU[,c(6,5,1:4)]
  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted[,10] < 0.05),]
  colnames(MLMM)[colnames(MLMM)=="MLMM.FNsmall"] <- "P.value"
  MLMM$Method <- paste("MLMM")
  MLMM$Trait.name <- paste("FNsmall")
  MLMM <- MLMM[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/4_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait4"] <- "CmD_LI_mm"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
    total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
    total$SNP.N <- paste("125")
    total <- unique(total[,c(7,1:6)])
    write.csv(total, file="Culmallresults125/CmD_LI_mm_2.csv",row.names = F)
    setwd("AllimagesC125/CmD_LI_mm")
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
###CmN.
###this one is for CmN. with single locus analysis 
###import trait CmN.
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.CmN..GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.CmN."
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.CmN..GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.CmN."
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.CmN..GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.CmN."
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.CmN..GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.CmN."
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,22)]
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


##using this one
p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.CmN. < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.CmN."] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("CmN")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.CmN. < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.CmN."] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("CmN")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.CmN. < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.CmN."] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("CmN")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.CmN. < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.CmN."] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("CmN")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.CmN. < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.CmN."] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("CmN")
str(p.rrBLUP)
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result
if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/CmN_1.csv",row.names = F)
}

###plot QQ plot inone image 
setwd("AllimagesC125/CmN")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=5.14401E-09,
       signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(2.08216E-07,5.14401E-09),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(2.08216E-07,5.14401E-09),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###this one is for CmN. with multiple locus analysis 
###this one is for CmN. with multiple locus analysis 
###import trait CmN.
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.CmN..GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.CmN."
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.CmN..GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.CmN."
intermediate <- read.csv("mrMLMM2/C125/5_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.CmN."
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.CmN."

###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.CmN."
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.CmN."
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.CmN.))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA), type = "left", by="SNP")
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

  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.MLMM.CmN < 0.05),]
  colnames(MLMM)[colnames(MLMM)=="MLMM.CmN."] <- "P.value"
  MLMM$Method <- paste("MLMM")
  MLMM$Trait.name <- paste("CmN")
  MLMM <- MLMM[,c(6,5,1:4)]
  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.FarmCPU.CmN < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.CmN."] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("CmN")
  FarmCPU <- FarmCPU[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/5_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait5"] <- "CmN"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]

    total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
    total$SNP.N <- paste("125")
    total <- unique(total[,c(7,1:6)])
    write.csv(total, file="Culmallresults125/CmN_2.csv",row.names = F)

    ###plot QQ plot inone image 
    setwd("AllimagesC125/Bcirc_cm")
    source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
    ###plot QQ plot inone image
    library("CMplot")
    CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=4.4e-6,
           signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
             TRUE,file="jpg",memo="",dpi=300,
           file.output = TRUE, verbose=TRUE)
    ###plot manhattan plot in one image
    CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
           multracks=TRUE, threshold=c(4.4e-6),threshold.lty=c(1,2), 
           threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
           chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
           file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
    ###this one plot manhttan plot in circur 
    CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
           threshold=c(4.4e-6),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
           signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
           bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
    

###Bcirc_cm
###this one is for Bcirc_cm with single locus analysis 
###import trait Bcirc_cm
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.Bcirc_cm.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.Bcirc_cm"
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.Bcirc_cm.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.Bcirc_cm"
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.Bcirc_cm.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.Bcirc_cm"
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.Bcirc_cm.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.Bcirc_cm"
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,23)]
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

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.Bcirc_cm < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.Bcirc_cm"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("Bcirc_cm")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.Bcirc_cm < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.Bcirc_cm"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("Bcirc_cm")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.Bcirc_cm < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.Bcirc_cm"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("Bcirc_cm")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.Bcirc_cm < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.Bcirc_cm"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("Bcirc_cm")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.Bcirc_cm < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.Bcirc_cm"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("Bcirc_cm")
str(p.rrBLUP)
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result
if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/Bcirc_cm_1.csv",row.names = F)
}

setwd("AllimagesC125/Bcirc_cm")
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

###this one is for Bcirc_cm with multiple locus analysis 
###this one is for Bcirc_cm with multiple locus analysis 
###import trait Bcirc_cm
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.Bcirc_cm.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.Bcirc_cm"
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.Bcirc_cm.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.Bcirc_cm"
intermediate <- read.csv("mrMLMM2/C125/6_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.Bcirc_cm"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.Bcirc_cm"
###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.Bcirc_cm"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.Bcirc_cm"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.Bcirc_cm))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA), type = "left", by="SNP")
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

  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,9] < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.Bcirc_cm"] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("Bcirc_cm")
  FarmCPU <-  FarmCPU[,c(6,5,1:4)]
  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted[,10] < 0.05),]
  colnames(MLMM)[colnames(MLMM)=="MLMM.Bcirc_cm"] <- "P.value"
  MLMM$Method <- paste("MLMM")   MLMM$Trait.name <- paste("Bcirc_cm")
  MLMM <- MLMM[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/6_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait6"] <- "Bcirc_cm"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
   names(mrmlm.final) %in% names(FarmCPU)
    total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
    total$SNP.N <- paste("125")
    total <- unique(total[,c(7,1:6)])
    write.csv(total, file="Culmallresults125/Bcirc_cm_2.csv",row.names = F)
    ###plot QQ plot inone image 
    setwd("AllimagesC125/Yld_kg")
    source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
    CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=1.3e-7,
           signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
             TRUE,file="jpg",memo="",dpi=300,
           file.output = TRUE, verbose=TRUE)
    ###plot manhattan plot in one image
    CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
           multracks=TRUE, threshold=c(1.3e-7),threshold.lty=c(1,2), 
           threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
           chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
           file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
    ###this one plot manhttan plot in circur 
    CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
           threshold=c(1.3e-7),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
           signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
           bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

####Yld_kg
###Yld_kg
###this one is for Yld_kg with single locus analysis 
###import trait Yld_kg
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.Yld_kg.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.Yld_kg"
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.Yld_kg.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.Yld_kg"
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.Yld_kg.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.Yld_kg"
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.Yld_kg.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.Yld_kg"
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,24)]
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


p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.Yld_kg < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.Yld_kg"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("Yld_kg")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.Yld_kg < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.Yld_kg"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("Yld_kg")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.Yld_kg < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.Yld_kg"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("Yld_kg")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.Yld_kg < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.Yld_kg"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("Yld_kg")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.Yld_kg < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.Yld_kg"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("Yld_kg")
str(p.rrBLUP)
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result
if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/Yld_kg_1.csv",row.names = F)
}



###this one is for Yld_kg with multiple locus analysis 
###this one is for Yld_kg with multiple locus analysis 
###import trait Yld_kg
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.Yld_kg.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.Yld_kg"
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.Yld_kg.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.Yld_kg"
intermediate <- read.csv("mrMLMM2/C125/7_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.Yld_kg"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.Yld_kg"

###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.Yld_kg"

###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.Yld_kg"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.Yld_kg))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA), type = "left", by="SNP")
plot[is.na(plot)] <- 1
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)



  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,9] < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.Yld_kg"] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("Yld_kg")
  FarmCPU <-  FarmCPU[,c(6,5,1:4)]
  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted[,10] < 0.05),]
  colnames(MLMM)[colnames(MLMM)=="MLMM.Yld_kg"] <- "P.value"
  MLMM$Method <- paste("MLMM")   
  MLMM$Trait.name <- paste("Yld_kg")
  MLMM <- MLMM[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/7_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait7"] <- "Yld_kg"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]

    total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
    total$SNP.N <- paste("125")
    total <- unique(total[,c(7,1:6)])
    write.csv(total, file="Culmallresults125/Yld_kg_2.csv",row.names = F)
    setwd("AllimagesC125/Yld_kg")
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
    
    

###SDW_kg
###this one is for SDW_kg with single locus analysis 
###import traitSDW_kg
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.SDW_kg.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.SDW_kg"
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.SDW_kg.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.SDW_kg"
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.SDW_kg.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.SDW_kg"
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.SDW_kg.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.SDW_kg"
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,25)]
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

##using this one
p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.SDW_kg < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.SDW_kg"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("SDW_kg")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.SDW_kg < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.SDW_kg"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("SDW_kg")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.SDW_kg < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.SDW_kg"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("SDW_kg")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.SDW_kg < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.SDW_kg"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("SDW_kg")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.SDW_kg < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.SDW_kg"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("SDW_kg")
str(p.rrBLUP)
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result
if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/SDW_kg_1.csv",row.names = F)
}

###plot QQ plot inone image 
setwd("AllimagesC125/SDW_kg")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.000175527,
       signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(8.70001E-05,0.000175527),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(8.70001E-05,0.000175527),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###this one is for SDW_kg with multiple locus analysis 
###this one is for SDW_kg with multiple locus analysis 
###import trait SDW_kg
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.SDW_kg.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.SDW_kg"
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.SDW_kg.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.SDW_kg"
intermediate <- read.csv("mrMLMM2/C125/8_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.SDW_kg"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.SDW_kg"
###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.SDW_kg"
###pKWmEB
pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
pKWmEB <- pKWmEB[2:3]
colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.SDW_kg"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.SDW_kg))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA,pKWmEB), type = "left", by="SNP")
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

  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,10] < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.SDW_kg"] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("SDW_kg")
  FarmCPU <-  FarmCPU[,c(6,5,1:4)]
  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted[,11] < 0.05),] 
  colnames(MLMM)[colnames(MLMM)=="MLMM.SDW_kg"] <- "P.value" 
  MLMM$Method <- paste("MLMM")  
  MLMM$Trait.name <- paste("SDW_kg") 
  MLMM <- MLMM[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/8_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait8"] <- "SDW_kg"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
    total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
    total$SNP.N <- paste("125")
    total <- unique(total[,c(7,1:6)])
    write.csv(total, file="Culmallresults125/CCirc_cm_2.csv",row.names = F)

    setwd("AllimagesC125/SDW_kg")
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
    
###CCirc_cm
###this one is for CCirc_cm with single locus analysis 
###import trait CCirc_cm
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.CCirc_cm.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.CCirc_cm"
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.CCirc_cm.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.CCirc_cm"
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.CCirc_cm.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.CCirc_cm"
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.CCirc_cm.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.CCirc_cm"
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,26)]
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
write.csv(P.less.0.05, file="AllimagesC125/CCirc_cm/CCirc_cm.csv")

##using this one
p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.CCirc_cm < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.CCirc_cm"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("CCirc_cm")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.CCirc_cm < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.CCirc_cm"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("CCirc_cm")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.CCirc_cm < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.CCirc_cm"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("CCirc_cm")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.CCirc_cm < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.CCirc_cm"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("CCirc_cm")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.CCirc_cm < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.CCirc_cm"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("CCirc_cm")
str(p.rrBLUP)
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result

if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/CCirc_cm_1.csv",row.names = F)
}

###plot QQ plot inone image 
setwd("AllimagesC125/CCirc_cm")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=2.4E-06,
       signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(2.4E-06),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(2.4E-06),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red"),
       signal.line=1,signal.col=c("red"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)




###this one is for CCirc_cm with multiple locus analysis 
###this one is for CCirc_cm with multiple locus analysis 
###import trait CCirc_cm
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.CCirc_cm.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.CCirc_cm"
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.CCirc_cm.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.CCirc_cm"
intermediate <- read.csv("mrMLMM2/C125/9_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.CCirc_cm"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.CCirc_cm"

###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.CCirc_cm"
###pKWmEB
pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
pKWmEB <- pKWmEB[2:3]
colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.CCirc_cm"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.CCirc_cm))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA,pKWmEB), type = "left", by="SNP")
plot[is.na(plot)] <- 1
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,10] < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.CCirc_cm"] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("CCirc_cm")
  FarmCPU <-  FarmCPU[,c(6,5,1:4)]
  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted[,11] < 0.05),] 
  colnames(MLMM)[colnames(MLMM)=="MLMM.CCirc_cm"] <- "P.value" 
  MLMM$Method <- paste("MLMM")   
  MLMM$Trait.name <- paste("CCirc_cm")  
  MLMM <- MLMM[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/9_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait9"] <- "CCirc_cm"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]

    total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
    total$SNP.N <- paste("125")
    total <- unique(total[,c(7,1:6)])
    write.csv(total, file="Culmallresults125/CCirc_cm_2.csv",row.names = F)
    
    setwd("AllimagesC125/CCirc_cm")
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
    

###Lg
###this one is for Lg with single locus analysis 
###import trait Lg
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.Lg.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.Lg"
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.Lg.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.Lg"
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.Lg.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.Lg"
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.Lg.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.Lg"
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,27)]
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
p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.Lg < 0.005),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.Lg"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("Lg")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.Lg < 0.005),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.Lg"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("Lg")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.Lg < 0.005),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.Lg"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("Lg")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.Lg < 0.0000005),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.Lg"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("Lg")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.Lg < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.Lg"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("Lg")
str(p.rrBLUP)
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result

if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/Lg_1.csv",row.names = F)
}

###plot QQ plot inone image 
setwd("AllimagesC125/Lg")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=10E-14,
       signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(10E-14,9.6E-11),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(10E-14,9.6E-11),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)



###this one is for Lg with multiple locus analysis 
###this one is for Lg with multiple locus analysis 
###import trait Lg
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.Lg.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.Lg"
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.Lg.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.Lg"
intermediate <- read.csv("mrMLMM2/C125/10_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.Lg"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.Lg"

###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.Lg"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.Lg"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.Lg))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA), type = "left", by="SNP")
plot[is.na(plot)] <- 1
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)

  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,9] < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.Lg"] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("Lg")
  FarmCPU <- FarmCPU[,c(6,5,1:4)]
  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted[,10] < 0.05),]
  colnames(MLMM)[colnames(MLMM)=="MLMM.Lg"] <- "P.value"
  MLMM$Method <- paste("MLMM")   
  MLMM$Trait.name <- paste("Lg")  
  MLMM <- MLMM[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/10_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait10"] <- "Lg"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]

    total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
    total$SNP.N <- paste("125")
    total <- unique(total[,c(7,1:6)])
    write.csv(total, file="Culmallresults125/Lg_2.csv",row.names = F)
    setwd("AllimagesC125/Lg")
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
    
###GS
###this one is for GS with single locus analysis 
###import trait GS
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.GS.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.GS"
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.GS.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.GS"
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.GS.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.GS"
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.GS.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.GS"
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,28)]
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
p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.Lg < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.GS"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("GS")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.Lg < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.GS"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("GS")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.GS < 0.01),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.GS"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("GS")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.GS < 0.000000000000000001),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.GS"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("GS")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.GS < 0.01),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.GS"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("GS")
str(p.rrBLUP)
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result
if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/GS_1.csv",row.names = F)
}

###plot QQ plot inone image
setwd("AllimagesC125/GS")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=9.68434E-24,
       signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(9.68434E-24,9.16386E-22),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(9.68434E-24,9.16386E-22),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###this one is for GS with multiple locus analysis 
###this one is for GS with multiple locus analysis 
###import trait GS
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.GS.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.GS"
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.GS.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.GS"
intermediate <- read.csv("mrMLMM2/C125/11_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.GS"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.GS"

###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.GS"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.GS"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.GS))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA), type = "left", by="SNP")
plot[is.na(plot)] <- 1
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)

  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted[,10] < 0.05),]
  colnames(MLMM)[colnames(MLMM)=="MLMM.GS"] <- "P.value"
  MLMM$Method <- paste("MLMM")
  MLMM$Trait.name <- paste("GS")
  MLMM <- MLMM[,c(6,5,1:4)]
  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,9] < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.GS"] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("GS")
  FarmCPU <- FarmCPU[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/11_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait11"] <- "GS"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  names(mrmlm.final) %in% names(FarmCPU)
    total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
    total$SNP.N <- paste("125")
    total <- unique(total[,c(7,1:6)])
    write.csv(total, file="Culmallresults125/GS_2.csv",row.names = F)

    setwd("AllimagesC125/GS")
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
    

###FD
###this one is for FD with single locus analysis 
###import trait FD
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FD.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FD"
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FD.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FD"
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.FD.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FD"
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.FD.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FD"
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,29)]
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
p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.GS < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.FD"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("FD")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.FD < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.FD"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("FD")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.FD < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.FD"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("FD")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.FD < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.FD"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("FD")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.FD < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.FD"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("FD")
str(p.rrBLUP)
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result
if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/FD_1.csv",row.names = F)
}

###plot QQ plot inone image
setwd("AllimagesC125/FD")
#library("CMplot")
source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=9.35197E-05,
       signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(9.35197E-05,0.000493303),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(9.35197E-05,0.000493303),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","blue"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)



###this one is for FD with multiple locus analysis 
###this one is for FD with multiple locus analysis 
###import trait FD
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.FD.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FD"
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.FD.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FD"
intermediate <- read.csv("mrMLMM2/C125/12_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.FD"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.FD"
###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.FD"
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA), type = "left", by="SNP")
plot[is.na(plot)] <- 1
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)

P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05)

  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,9] < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.FD"] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("FD")
  FarmCPU <-  FarmCPU[,c(6,5,1:4)]
  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted[,10] < 0.05),]
  colnames(MLMM)[colnames(MLMM)=="MLMM.FD"] <- "P.value"
  MLMM$Method <- paste("MLMM")
  MLMM$Trait.name <- paste("FD")
  MLMM <- MLMM[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/12_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait12"] <- "FD"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]

    total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
    total$SNP.N <- paste("125")
    total <- unique(total[,c(7,1:6)])
    write.csv(total, file="Culmallresults125/FD_2.csv",row.names = F)

    setwd("AllimagesC125/FD")
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
    
    
###SRD
###this one is for SRD with single locus analysis 
###import trait SRD
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.SRD.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.SRD"
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.SRD.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.SRD"
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.SRD.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.SRD"
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.SRD.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.SRD"
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,30)]
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
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,9] < 0.05 | 
                        plot.p.adjusted[,10] < 0.05 | plot.p.adjusted[,11] < 0.05 |
                        plot.p.adjusted[,12] < 0.05 | plot.p.adjusted[,13] < 0.05)
str(P.less.0.05)

p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.FD < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.SRD"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("SRD")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.SRD < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.SRD"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("SRD")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.SRD < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.SRD"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("SRD")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.SRD < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.SRD"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("SRD")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.SRD < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.SRD"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("SRD")
str(p.rrBLUP)
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result
if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/SRD_1.csv",row.names = F)
}

setwd("AllimagesC125/SRD")
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


###this one is for SRD with multiple locus analysis 
###this one is for SRD with multiple locus analysis 
###import trait SRD
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.SRD.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.SRD"
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.SRD.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.SRD"
intermediate <- read.csv("mrMLMM2/C125/13_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.SRD"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.SRD"

###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.SRD"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.SRD"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.SRD))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA), type = "left", by="SNP")
plot[is.na(plot)] <- 1
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)


  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,9] < 0.05),]
  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.MLMM.SRD < 0.05),]
  colnames(MLMM)[colnames(MLMM)=="MLMM.SRD"] <- "P.value"
  MLMM$Method <- paste("MLMM")
  MLMM$Trait.name <- paste("SRD")
  MLMM <- MLMM[,c(6,5,1:4)]
  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.FarmCPU.SRD < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.SRD"] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("SRD")
  FarmCPU <- FarmCPU[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/13_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait13"] <- "SRD"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
    total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
    total$SNP.N <- paste("125")
    total <- unique(total[,c(7,1:6)])
    write.csv(total, file="Culmallresults125/SRD_2.csv",row.names = F)
 
    ###plot QQ plot inone image
    setwd("AllimagesC125/SRD")
    #library("CMplot")
    source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")
    ###plot QQ plot inone image
    library("CMplot")
    CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.0000033,
           signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
             TRUE,file="jpg",memo="",dpi=300,
           file.output = TRUE, verbose=TRUE)
    ###plot manhattan plot in one image
    CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
           multracks=TRUE, threshold=c(0.0000033),threshold.lty=c(1,2), 
           threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
           chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
           file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
    ###this one plot manhttan plot in circur 
    CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
           threshold=c(0.0000033),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
           signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
           bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###ADD
###this one is for ADD with single locus analysis 
###import trait ADD
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.GLM.ADD.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.ADD"
str(GLM)
MLM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPUCV/GAPIT.MLM.ADD.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.ADD"
CMLM <- read.csv("Result GAPIT1/C125ECMMLCV/GAPIT.CMLM.ADD.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.ADD"
MLMSUPER <- read.csv("Result GAPIT1/C125MLMSUPPER/GAPIT.SUPER.ADD.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.ADD"
gwasResultsqq <- read.csv("rrBLUPC125/rrBLUP_GWAS_results.C125.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,31)]
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
p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.ADD < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.ADD"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("ADD")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.ADD < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.ADD"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("ADD")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.ADD < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.ADD"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("ADD")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.ADD < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.ADD"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("ADD")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.ADD < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.ADD"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("ADD")
str(p.rrBLUP)
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("125")
total <- total[,c(7,6,5,1:4)]
###write the result

if (nrow(total)>0) {
  write.csv(total, file="Culmallresults125/ADD_1.csv",row.names = F)
}

###this one is for ADD with multiple locus analysis 
###this one is for ADD with multiple locus analysis 
###import trait ADD
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.FarmCPU.ADD.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.ADD"
MLMM <- read.csv("Result GAPIT1/C125GLMMLMMLMMFarmCPU/GAPIT.MLMM.ADD.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.ADD"
intermediate <- read.csv("mrMLMM2/C125/14_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)
colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.ADD"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.ADD"
###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.ADD"
###pKWmEB
#pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
#pKWmEB <- pKWmEB[2:3]
#colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.ADD"
#subset(pKWmEB,is.na(pKWmEB$pKWmEB.ADD))
###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM,FASTmrEMMA), type = "left", by="SNP")
plot[is.na(plot)] <- 1
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 5)


  FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,9] < 0.05),]
  colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.ADD"] <- "P.value"
  FarmCPU$Method <- paste("FarmCPU")
  FarmCPU$Trait.name <- paste("ADD")
  FarmCPU <- FarmCPU[,c(6,5,1:4)]

  MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted[,10] < 0.05),]   
  colnames(MLMM)[colnames(MLMM)=="MLMM.ADD"] <- "P.value"   
  MLMM$Method <- paste("MLMM")   
  MLMM$Trait.name <- paste("ADD")   
  MLMM <- MLMM[,c(6,5,1:4)]
  mrmlm.final <- read.csv("mrMLMM2/C125/14_final result.csv")
  colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
  colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
  levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait14"] <- "ADD"
  mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
  mrmlm.final <- mrmlm.final[,c(2:6,15)]
  #total <-  plyr::rbind.fill(list(mrmlm.final,FarmCPU, MLMM))
  names(FarmCPU) %in% names(mrmlm.final)
    total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
    total$SNP.N <- paste("125")
    total <- unique(total[,c(7,1:6)])
    write.csv(total, file="Culmallresults125/ADD_2.csv",row.names = F)
    setwd("AllimagesC125/ADD")
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
    

