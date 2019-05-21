####import the Suve data of 2015
Suv15 <- read.csv("data/15Suv.csv")
str(Suv15)
colnames(Suv15)[which(colnames(Suv15)=="Plot")] <- "plot"
str(Suv15)
Suv15$Surv.25[(Suv15$Surv.25) == 1] <- 0
Suv15$Surv.25[(Suv15$Surv.25) == 2] <- 3
all <- read.csv("data/trait1718.3.16.19withoutliers.csv")
all17 <- subset(all, all$Year==2017)
all17 <- all17[order(all17$plot1),]
str(all17)
all1517 <- plyr::join_all(list(all17,Suv15), by="plot")
all1517$Surv.25

all1517$OWA <- all1517$Surv.25 - all1517$Surv
str(all1517$OWA)
levels(as.factor(all1517$OWA))
all1517$OWA[(all1517$OWA)==-1] <- NA
all1517$OWA[(all1517$OWA)== -2] <- 1
all1517$OWA[(all1517$OWA)== 2] <- 0
all1517$OWA[(all1517$OWA)== 1] <- 1


all18 <- subset(all, all$Year==2018)

all18 <- all18[order(all18$plot1),]
str(all18)
all1518 <- plyr::join_all(list(all18,Suv15), by="plot")
all1518$Surv.25

all1518$OWA <- all1518$Surv.25 - all1518$Surv
str(all1518$OWA)
levels(as.factor(all1518$OWA))
all1518$OWA[(all1518$OWA)==-1] <- NA
all1518$OWA[(all1518$OWA)== -2] <- 1
all1518$OWA[(all1518$OWA)== 2] <- 0
all1518$OWA[(all1518$OWA)== 1] <- 1
###meger two data set from vertical 
install.packages("lessR")
library(lessR)
trait1718.3.16.19withoutliers1 <- Merge(all1518,all1517)
trait1718.3.16.19withoutliers1$Surv <- trait1718.3.16.19withoutliers1$OWA

trait1718.3.16.19withoutliers1 <- trait1718.3.16.19withoutliers1[1:30]
colnames(trait1718.3.16.19withoutliers1)[colnames(trait1718.3.16.19withoutliers1)=="Surv"] <- "OWA"
str(trait1718.3.16.19withoutliers1)
write.csv(trait1718.3.16.19withoutliers1,file="data/trait1718.3.16.19withoutliers1.csv")


