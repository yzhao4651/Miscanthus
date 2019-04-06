###this one for flo
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/mrMLMflo2")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("FASTmrMLM"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/FASTmrMLMflo2")
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("FASTmrEMMA"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/FASTmrEMMAflo2")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("pLARmEB"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/pLARmEBflo2")
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/ISIS EM-BLASSOflo2")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("pKWmEB"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/pKWmEBflo2")

####including all of the methods: 
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo2.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/Resultfloall2")

###using the the third time produced phenotype
###using the the third time produced phenotype
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/mrMLMflo3")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("FASTmrMLM"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/FASTmrMLMflo3")
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("FASTmrEMMA"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/FASTmrEMMAflo3")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("pLARmEB"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/pLARmEBflo3")
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/ISIS EM-BLASSOflo3")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("pKWmEB"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/pKWmEBflo3")

####including all of the methods: 

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMflo3.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMflo3.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:23,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "mrMLMM2/Resultfloall3")
