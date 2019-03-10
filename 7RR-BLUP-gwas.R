####rrBLUP method
####import the genotype data
myY <- read.csv("data/myYimputedSNP19.csv")
myGD <- read.csv("data/myGDimputedSNP19.csv",row.names=1)
myGM <- read.csv("data/myGMimputedSNP19.csv")
myQ<- read.csv("data/myQimputedSNP19.csv",row.names=1)

###tansform the genotype data set into rrblup data set
myGDt <- t(myGD)-1
#myGDtdataframe <- data.frame(myGDt)
#head(myGDtdataframe)

### change the row name into the first column 
myGDtdataframe <- data.frame(myGDt)
Name <- rownames(myGDtdataframe)
rownames(myGDtdataframe) <- NULL
myGDtdataframe <- cbind(Name,myGDtdataframe)

###combine the genotype with the chromosome position (dataset: myGM)
genorrblup <- merge(myGM,myGDtdataframe,by="Name")
###change the ID of the name of phenotype to same to the ID of the genotype ( like change "_" to "." )
head(myY)
as.character(myY$Taxa)
myY$Taxa <- gsub("-",".", myY$Taxa)
head(myY)
str(myY)
###run rrblup
###run rrblup
install.packages("rrBLUP")


###1:question, I can not open the pdf, even opened it, only one page show image. 
###I tried get off the pdf and dev.off() and run again, but the plot did not show the Plots area. 
 
###2: how to output the result in to csv file? 
###3: I set up pdf(paste("Plots from rrBLUP", 1, ".pdf", sep="")), 
###but i changed the number with 2 next time like pdf(paste("Plots from rrBLUP", 2, ".pdf", sep="")), 
####and run it agian, but the 1 was changed. Would you please give me some advice? thanks 

### this one is without Kin matrix, 
library(rrBLUP)
pdf(paste("Plots from rrBLUP", 1, ".pdf", sep=""))
par(mfrow=c(25,2))
GWAS(myY[1:2], genorrblup, fixed=NULL, K=NULL, n.PC=0,
     min.MAF=0.05, n.core=1, P3D=TRUE, plot=TRUE)
dev.off()


###tryig to calculate the kinship 
###change the first column into row name
myGDk <- data.frame(myGD)-1
dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
kmatrix <- A.mat(myGDk,min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,
      n.core=1,shrink=FALSE,return.imputed=FALSE)
####do GWAS analysis again with Kmatrix
library(rrBLUP)
GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=0,
     min.MAF=0.05, n.core=1, P3D=TRUE, plot=TRUE)


