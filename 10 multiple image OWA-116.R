
###import trait OWA
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
GLM <- read.csv("Result GAPIT1/116OWAGLMMLMMLMMFarmCPUCV/GAPIT.GLM.OWA.GWAS.Results.csv")
GLM <- GLM[1:4]
colnames(GLM)[colnames(GLM)=="P.value"] <- "GLM.OWA"
MLM <- read.csv("Result GAPIT1/116OWAGLMMLMMLMMFarmCPUCV/GAPIT.MLM.OWA.GWAS.Results.csv")
MLM <- MLM[,c(1,4)]
colnames(MLM)[colnames(MLM)=="P.value"] <- "MLM.OWA"
CMLM <- read.csv("Result GAPIT1/116OWAECMMLCV/GAPIT.CMLM.OWA.GWAS.Results.csv")
CMLM <- CMLM[,c(1,4)]
colnames(CMLM)[colnames(CMLM)=="P.value"] <- "CMLM.OWA"
MLMSUPER <- read.csv("Result GAPIT1/116OWAMLMSUPPER/GAPIT.SUPER.OWA.GWAS.Results.csv")
MLMSUPER<- MLMSUPER[,c(1,4)]
colnames(MLMSUPER)[colnames(MLMSUPER)=="P.value"] <- "SUPER.OWA"
gwasResultsqq <- read.csv("rrBLUPS116/rrBLUP_GWAS_results.O116.csv",row.names = 1)
gwasResultsqq <- gwasResultsqq[, c(1,5)]
colnames(gwasResultsqq)[colnames(gwasResultsqq)=="Name"] <- "SNP"
names(gwasResultsqq)
###merge all of them 
plot <-  plyr::join_all(list(GLM,MLM,CMLM,MLMSUPER,gwasResultsqq), by="SNP")
###using this function to adjust the p value 
source("Function/adj_P_function.R")
plot.p.adjusted <- adj_P_function(plot, 4, 8)
###and get all of the SNPS suitable for the p-value
p.GLM <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted$Adj.P.GLM.OWA < 0.05),]
colnames(p.GLM)[colnames(p.GLM)=="GLM.OWA"] <- "P.value"
p.GLM$Method <- paste("GLM")
p.GLM$Trait.name <- paste("OWA")
str(p.GLM)
p.MLM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted$Adj.P.MLM.OWA < 0.05),]
colnames(p.MLM)[colnames(p.MLM)=="MLM.OWA"] <- "P.value"
p.MLM$Method <- paste("MLM")
p.MLM$Trait.name <- paste("OWA")
str(p.MLM)
p.CMLM <- plot.p.adjusted[,c(1:3,6)][which(plot.p.adjusted$Adj.P.CMLM.OWA < 0.05),]
colnames(p.CMLM)[colnames(p.CMLM)=="CMLM.OWA"] <- "P.value"
p.CMLM$Method <- paste("CMLM")
p.CMLM$Trait.name <- paste("OWA")
str(p.CMLM)
p.SUPER <- plot.p.adjusted[,c(1:3,7)][which(plot.p.adjusted$Adj.P.SUPER.OWA < 0.05),]
colnames(p.SUPER)[colnames(p.SUPER)=="SUPER.OWA"] <- "P.value"
p.SUPER$Method <- paste("SUPER")
p.SUPER$Trait.name <- paste("OWA")
str(p.SUPER)
p.rrBLUP <- plot.p.adjusted[,c(1:3,8)][which(plot.p.adjusted$Adj.P.rrBLUP.OWA < 0.05),]
colnames(p.rrBLUP)[colnames(p.rrBLUP)=="rrBLUP.OWA"] <- "P.value"
p.rrBLUP$Method <- paste("rrBLUP")
p.rrBLUP$Trait.name <- paste("OWA")
str(p.rrBLUP)
total <- do.call("rbind", list(p.GLM, p.CMLM, p.SUPER, p.rrBLUP))
total$SNP.N <- paste("116")
total <- total[,c(7,6,5,1:4)]
str(total)
###write the result
if (nrow(total)>0) {
write.csv(total, file="Owaresults106/OWA_1.csv",row.names = F)
}
###plot QQ plot inone image
setwd("Allimage116OWA")
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

###this one is for OWA with multiple locus analysis 
###this one is for OWA with multiple locus analysis 
###import trait OWA
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
FarmCPU <- read.csv("Result GAPIT1/116OWAGLMMLMMLMMFarmCPU/GAPIT.FarmCPU.OWA.GWAS.Results.csv")
FarmCPU <- FarmCPU[1:4]
colnames(FarmCPU)[colnames(FarmCPU)=="P.value"] <- "FarmCPU.OWA"
MLMM <- read.csv("Result GAPIT1/116OWAGLMMLMMLMMFarmCPU/GAPIT.MLMM.OWA.GWAS.Results.csv")
#MLMM <- read.csv("Result GAPIT1/MLMM/GAPIT.MLMM.OWA.GWAS.Results.csv")
MLMM <- MLMM[,c(1,4)]
colnames(MLMM)[colnames(MLMM)=="P.value"] <- "MLMM.OWA"
intermediate <- read.csv("mrMLMM2/OWA106/1_intermediate result.csv")
intermediate$X.log10.P. <- 10^-(intermediate$X.log10.P.)

colnames(intermediate)[colnames(intermediate)=="RS."] <- "SNP"
intermediate <- intermediate[,c(3:4,8)]
levels(intermediate$Method)
###FASTmrMLM
FASTmrMLM <- subset(intermediate,intermediate$Method=="FASTmrMLM")
FASTmrMLM <- FASTmrMLM[2:3]
colnames(FASTmrMLM)[colnames(FASTmrMLM)=="X.log10.P."] <- "FASTmrMLM.OWA"
###mrMLM
mrMLM <- subset(intermediate,intermediate$Method=="mrMLM")
mrMLM <- mrMLM[2:3]
colnames(mrMLM)[colnames(mrMLM)=="X.log10.P."] <- "mrMLM.OWA"
###FASTmrEMMA
FASTmrEMMA <- subset(intermediate,intermediate$Method=="FASTmrEMMA")
FASTmrEMMA <- FASTmrEMMA[2:3]
colnames(FASTmrEMMA)[colnames(FASTmrEMMA)=="X.log10.P."] <- "FASTmrEMMA.OWA"

###pKWmEB
pKWmEB <- subset(intermediate,intermediate$Method=="pKWmEB")
pKWmEB <- pKWmEB[2:3]
colnames(pKWmEB)[colnames(pKWmEB)=="X.log10.P."] <- "pKWmEB.OWA"

###merge all of them 
plot <-  plyr::join_all(list(FarmCPU,MLMM,mrMLM), type = "left", by="SNP")
plot[is.na(plot)] <- 1

source("Function/adj_P_functionb.R")
plot.p.adjusted <- adj_P_function(plot, 4, 5)
str(plot.p.adjusted)

FarmCPU <- plot.p.adjusted[,c(1:4)][which(plot.p.adjusted[,7] < 0.05),]
colnames(FarmCPU)[colnames(FarmCPU)=="FarmCPU.OWA"] <- "P.value"
FarmCPU$Method <- paste("FarmCPU")
FarmCPU$Trait.name <- paste("OWA")
FarmCPU <- FarmCPU[,c(6,5,1:4)]
MLMM <- plot.p.adjusted[,c(1:3,5)][which(plot.p.adjusted[,8] < 0.05),]
colnames(MLMM)[colnames(MLMM)=="MLMM.OWA"] <- "P.value"
MLMM$Method <- paste("MLMM")
MLMM$Trait.name <- paste("OWA")
MLMM <- MLMM[,c(6,5,1:4)]
mrmlm.final <- read.csv("mrMLMM2/OWA106/1_Final result.csv")
colnames(mrmlm.final)[colnames(mrmlm.final) == "RS."] <- "SNP"
colnames(mrmlm.final)[colnames(mrmlm.final) == "Marker.Position..bp."] <- "Position"
levels(mrmlm.final$Trait.name)[levels(mrmlm.final$Trait.name) == "Trait1"] <- "OWA"
mrmlm.final$P.value <- 10^-(mrmlm.final$X.log10.P.)
mrmlm.final <- mrmlm.final[,c(2:6, 15)]
total <- do.call("rbind", list(mrmlm.final,FarmCPU, MLMM))
total$SNP.N <- paste("116")
total <- unique(total[,c(7,1:6)])
###write the result
write.csv(total, file="Owaresults106/OWA_2.csv",row.names = F)
setwd("Allimage116OWA")
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


