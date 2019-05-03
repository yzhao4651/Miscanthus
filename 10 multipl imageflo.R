
### Multiple image together for flowering related traits
### The path i set up probably hard for you.  
### Every time, I set up the path at beginning,i need full path for import the result. 
### how to edit it?
### I did trait one by one, do you have any idea that i can do it on one time. 
###TFN.
###this one is for TFN. with single locus analysis 
###import trait TFN.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/TFN.")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.TFN..GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.TFN."
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.TFN..GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.TFN."
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.TFN..GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.TFN."
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.TFN..GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.TFN."
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,54)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for TFN. with multiple locus analysis 
###this one is for TFN. with multiple locus analysis 
###import trait TFN.
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.TFN..GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.TFN."
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.TFN..GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.TFN."
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/1_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###FNMain
###this one is for FNMain with single locus analysis 
###import trait FNMain
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/FNMain")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FNMain.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FNMain"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FNMain.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FNMain"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.FNMain.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FNMain"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FNMain.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FNMain"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,55)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FNMain with multiple locus analysis 
###this one is for FNMain with multiple locus analysis 
###import trait FNMain
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FNMain.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FNMain"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FNMain.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FNMain"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/2_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###FNsmall
###this one is for FNsmall with single locus analysis 
###import trait FNsmall
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/FNsmall")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FNsmall.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FNsmall"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FNsmall.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FNsmall"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.FNsmall.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FNsmall"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FNsmall.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FNsmall"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,56)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FNsmall with multiple locus analysis 
###this one is for FNsmall with multiple locus analysis 
###import trait FNsmall
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FNsmall.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FNsmall"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FNsmall.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FNsmall"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/3_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###HD_1
###this one is for HD_1 with single locus analysis 
###import trait HD_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/HD_1")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.HD_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.HD_1"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.HD_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.HD_1"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.HD_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.HD_1"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.HD_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.HD_1"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,60)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for HD_1 with multiple locus analysis 
###this one is for HD_1 with multiple locus analysis 
###import trait HD_1
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.HD_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.HD_1"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.HD_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.HD_1"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/4_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###FD_1
###this one is for FD_1 with single locus analysis 
###import trait FD_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/FD_1")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FD_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FD_1"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FD_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FD_1"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.FD_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FD_1"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FD_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FD_1"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,61)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FD_1 with multiple locus analysis 
###this one is for FD_1 with multiple locus analysis 
###import trait FD_1
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FD_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FD_1"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FD_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FD_1"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/5_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###HD_50.
###this one is for HD_50. with single locus analysis 
###import traitHD_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/HD_50.")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.HD_50..GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.HD_50."
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.HD_50..GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.HD_50."
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.HD_50..GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.HD_50."
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.HD_50..GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.HD_50."
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,62)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for HD_50. with multiple locus analysis 
###this one is for HD_50. with multiple locus analysis 
###import trait HD_50.
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.HD_50..GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.HD_50."
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.HD_50..GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.HD_50."
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/6_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###FD_50.
###this one is for FD_50. with single locus analysis 
###import trait FD_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/FD_50.")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FD_50..GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FD_50."
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FD_50..GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FD_50."
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.FD_50..GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FD_50."
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FD_50..GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FD_50."
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,63)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FD_50. with multiple locus analysis 
###this one is for FD_50. with multiple locus analysis 
###import trait FD_50.
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FD_50..GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FD_50."
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FD_50..GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FD_50."
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/7_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)




###HW_1
###this one is for HW_1 with single locus analysis 
###import trait HW_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/HW_1")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.HW_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.HW_1"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.HW_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.HW_1"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.HW_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.HW_1"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.HW_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.HW_1"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,64)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for HW_1 with multiple locus analysis 
###this one is for HW_1 with multiple locus analysis 
###import trait HW_1
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.HW_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.HW_1"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.HW_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.HW_1"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/8_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###FW_1
###this one is for FW_1 with single locus analysis 
###import trait FW_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/FW_1")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FW_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FW_1"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FW_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FW_1"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.FW_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FW_1"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FW_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FW_1"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,65)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FGW_1 with multiple locus analysis 
###this one is for FGW_1 with multiple locus analysis 
###import trait FGW_1
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FW_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FW_1"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FW_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FW_1"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/9_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###HW_50.
###this one is for HW_50. with single locus analysis 
###import trait HW_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/HW_50.")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.HW_50.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.HW_50."
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.HW_50.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.HW_50."
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.HW_50.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.HW_50."
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.HW_50.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.HW_50."
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,66)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for HW_50. with multiple locus analysis 
###this one is for HW_50. with multiple locus analysis 
###import trait HW_50.
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.HW_50.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.HW_50."
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.HW_50.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.HW_50."
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/10_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###FW_50.
###this one is for FW_50. with single locus analysis 
###import trait FW_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/FW_50.")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FW_50.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FW_50."
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FW_50.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FW_50."
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.FW_50.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FW_50."
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FW_50.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FW_50."
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,67)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FW_50. with multiple locus analysis 
###this one is for FW_50. with multiple locus analysis 
###import trait FW_50.
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FW_50.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FW_50."
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FW_50.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FW_50."
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/11_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)



###GHW_1
###this one is for GHW_1 with single locus analysis 
###import trait GHW_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/GHW_1")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.GHW_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.GHW_1"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.GHW_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.GHW_1"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.GHW_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.GHW_1"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.GHW_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.GHW_1"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,68)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for GHW_1 with multiple locus analysis 
###this one is for GHW_1 with multiple locus analysis 
###import trait GHW_1
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.GHW_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.GHW_1"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.GHW_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.GHW_1"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/12_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###GFW_1
###this one is for GFW_1 with single locus analysis 
###import trait GFW_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/GFW_1")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.GFW_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.GFW_1"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.GFW_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.GFW_1"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.GFW_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.GFW_1"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.GFW_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.GFW_1"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,69)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for GFW_1 with multiple locus analysis 
###this one is for GFW_1 with multiple locus analysis 
###import trait GFW_1
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.GFW_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.GFW_1"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.GFW_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.GFW_1"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/13_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###GHW_50.
###this one is for GHW_50. with single locus analysis 
###import trait GHW_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/GHW_50.")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.GHW_50.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.GHW_50."
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.GHW_50.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.GHW_50."
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.GHW_50.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.GHW_50."
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.GHW_50.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.GHW_50."
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,70)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for GHW_50. with multiple locus analysis 
###this one is for GHW_50. with multiple locus analysis 
###import trait GHW_50.
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.GHW_50.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.GHW_50."
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.GHW_50.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.GHW_50."
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/14_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###GFW_50.
###this one is for GFW_50. with single locus analysis 
###import trait GFW_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/GFW_50.")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.GFW_50.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.GFW_50."
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.GFW_50.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.GFW_50."
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.GFW_50.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.GFW_50."
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.GFW_50.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.GFW_50."
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,71)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for GFW_50. with multiple locus analysis 
###this one is for GFW_50. with multiple locus analysis 
###import trait GFW_50.
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.GFW_50.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.GFW_50."
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.GFW_50.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.GFW_50."
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/15_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###HM_1
###this one is for HM_1 with single locus analysis 
###import trait HM_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/HM_1")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.HM_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.HM_1"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.HM_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.HM_1"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.HM_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.HM_1"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.HM_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.HM_1"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,72)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for HM_1 with multiple locus analysis 
###this one is for HM_1 with multiple locus analysis 
###import trait HM_1
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.HM_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.HM_1"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.HM_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.HM_1"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/16_intermediate result.csv")
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
plot <-  plyr::join_all(list(FarmCPU,MLMM,FASTmrMLM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###FM_1
###this one is for FM_1 with single locus analysis 
###import trait FM_1
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/FM_1")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FM_1.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FM_1"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FM_1.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FM_1"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.FM_1.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FM_1"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FM_1.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FM_1"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,73)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FM_1 with multiple locus analysis 
###this one is for FM_1 with multiple locus analysis 
###import trait FM_1
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FM_1.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FM_1"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FM_1.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FM_1"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/17_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###HM_50.
###this one is for HM_50. with single locus analysis 
###import traitHM_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/HM_50.")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.HM_50.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.HM_50."
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.HM_50.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.HM_50."
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.HM_50.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.HM_50."
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.HM_50.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.HM_50."
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,74)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for HM_50. with multiple locus analysis 
###this one is for HM_50. with multiple locus analysis 
###import trait HM_50.
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.HM_50.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.HM_50."
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.HM_50.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.HM_50."
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/18_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###FM_50.
###this one is for FM_50. with single locus analysis 
###import trait FM_50.
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/FM_50.")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.FM_50.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.FM_50."
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.FM_50.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.FM_50."
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.FM_50.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.FM_50."
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.FM_50.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.FM_50."
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,75)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.7,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for FM_50. with multiple locus analysis 
###this one is for FM_50. with multiple locus analysis 
###import trait FM_50.
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.FM_50.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.FM_50."
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.FM_50.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.FM_50."
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/19_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###fprind
###this one is for fprind with single locus analysis 
###import trait fprind
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/fprind")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.fprind.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.fprind"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.fprind.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.fprind"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.fprind.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.fprind"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.fprind.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.fprind"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,76)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for fprind with multiple locus analysis 
###this one is for fprind with multiple locus analysis 
###import trait fprind
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprind.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.fprind"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.fprind.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.fprind"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/20_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


###fprinW
###this one is for fprinW with single locus analysis 
###import trait fprinW
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/fprinW")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.fprinW.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.fprinW"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.fprinW.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.fprinW"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.fprinW.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.fprinW"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.fprinW.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.fprinW"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,77)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for fprinW with multiple locus analysis 
###this one is for fprinW with multiple locus analysis 
###import trait fprinW
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinW.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.fprinW"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.fprinW.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.fprinW"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/21_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###fprinGW
###this one is for fprinGW with single locus analysis 
###import trait fprinGW
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/fprinGW")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.fprinGW.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.fprinGW"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.fprinGW.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.fprinGW"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.fprinGW.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.fprinGW"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.fprinGW.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.fprinGW"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,78)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for fprinGW with multiple locus analysis 
###this one is for fprinGW with multiple locus analysis 
###import trait fprinGW
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinGW.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.fprinGW"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.fprinGW.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.fprinGW"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/22_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###fprinM
###this one is for fprinM with single locus analysis 
###import trait fprinM
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Allimages/fprinM")
GLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.fprinM.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.fprinM"
str(GLM)
MLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.fprinM.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.fprinM"
CMLM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/ECMMLCV/GAPIT.CMLM.fprinM.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.fprinM"
MLMSUPER <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/MLMSUPPER/GAPIT.SUPER.fprinM.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.fprinM"
gwasResultsqq <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,79)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###plot QQ plot inone image
library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

###this one is for fprinM with multiple locus analysis 
###this one is for fprinM with multiple locus analysis 
###import trait fprinM
FarmCPU <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprinM.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
str(FarmCPU)
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.fprinM"
MLMM <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLMM.fprinM.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.fprinM"
intermediate <- read.csv("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Resultfloall1/23_intermediate result.csv")
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

library("CMplot")
CMplot(plot,plot.type="q",col=c("blue", "orange", "cyan","magenta","green3"),threshold=0.05/nrow(plot),
       signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
         TRUE,file="jpg",memo="",dpi=300,
       file.output = TRUE, verbose=TRUE)
###plot manhattan plot in one image
CMplot(plot, plot.type="m", col=c("magenta", "orange", "cyan","green3","blue"),
       multracks=TRUE, threshold=c(1e-6,5.8e-05),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","blue"),signal.cex=c(1,1),
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)
###this one plot manhttan plot in circur 
CMplot(plot,plot.type="c",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:18,"X"),sep=""),
       threshold=c(1e-6,5.8e-05),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),
       signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


