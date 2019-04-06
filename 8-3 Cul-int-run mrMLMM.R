library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/mrMLMculm2")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("FASTmrMLM"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/FASTmrMLMculm2")
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("FASTmrEMMA"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/FASTmrEMMAculm2")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("pLARmEB"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/pLARmEBculm2")
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/ISIS EM-BLASSOculm2")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("pKWmEB"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/pKWmEBculm2")

####including all of the methods: 

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/Resultculmall2")

####using the genotype produced at the third time
####using the genotype produced at the third time
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/mrMLMculm3")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("FASTmrMLM"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/FASTmrMLMculm3")
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("FASTmrEMMA"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/FASTmrEMMAculm3")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("pLARmEB"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/pLARmEBculm3")
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/ISIS EM-BLASSOculm3")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("pKWmEB"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/pKWmEBculm3")

####including all of the methods: 

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMculm3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/Resultculmall3")