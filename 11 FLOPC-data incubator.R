
###combine all of the result from FarmCPU
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
FarmCPU <- read.csv("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.FarmCPU.fprind.GWAS.Results.csv")
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
write.csv(plot, file="Allimages/FarmCPU.all.csv")
str(plot)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="bonferroni", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(plot, 4, 7)
str(plot.p.adjusted)
P.less.0.05 <- subset(plot.p.adjusted, plot.p.adjusted[,8] < 0.05 | 
                        plot.p.adjusted[,9] < 0.05 | plot.p.adjusted[,10] < 0.05|
                        plot.p.adjusted[,11] < 0.05)
str(P.less.0.05)
write.csv(plot, file="Allimages/FarmCPU.p.less.0.05.csv")

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

