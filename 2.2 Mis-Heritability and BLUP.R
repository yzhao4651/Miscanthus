###import the data
###import the data
#normadata <- read.csv("~/Documents/whole traits/traits1718normalited1.csv",na.strings = c("",".","NA"))
normadata <- read.csv("data/traits1718normalited1.csv",na.strings = c("",".","NA"),row.names=1)
###check the data format
str(normadata)

###change the format of the several variables 
normadata$GS <- as.numeric(as.character(normadata$GS))
normadata$SRD <- as.numeric(as.character(normadata$SRD))
normadata$Entry=as.factor(normadata$Entry)
normadata$Rep=as.factor(normadata$Rep)
normadata$Year=as.factor(normadata$Year)
###check the data format
str(normadata)
###get the heritablity of all the traits using the function in the Function
source("Function/Heritability.R")
Herit <- Heritability(4,25,normadata)
write.csv(Herit, file = "data/heritabilityall.csv", row.names = T, na = ".")

###get the heritablity of all the traits using the function in the Function
source("Function/BLUP.R")
ranefvalueall <- ranefvalue(4,25,normadata)
write.csv(ranefvalueall, file = "data/ranefvalueall.csv", row.names = T, na = ".")

###only for OWA
###only for OWA
#normadata <- read.csv("~/Documents/whole traits/traits1718normalited1.csv",na.strings = c("",".","NA"))
normadata <- read.csv("data/traits1718normalited3.csv",na.strings = c("",".","NA"),row.names=1)
###check the data format
normadata <- normadata[2:5]
str(normadata)
###change the format of the several variables 
normadata$Entry=as.factor(normadata$Entry)
normadata$Rep=as.factor(normadata$Rep)
normadata$Year=as.factor(normadata$Year)
###check the data format
str(normadata)
###get the heritablity of all the traits using the function in the Function
source("Function/Heritability.R")
Herit <- Heritability(4, 4, normadata)
write.csv(Herit, file = "data/heritabilitOWA.csv", row.names = T, na = ".")
###get the heritablity of all the traits using the function in the Function
source("Function/BLUP.R")
ranefvalueall <- ranefvalue(4, 4, normadata)
write.csv(ranefvalueall, file = "data/ranefvalueOWA.csv", row.names = T, na = ".")

###only for flowering traits
###only for flowering traits
normadata <- read.csv("data/flowpc.csv",na.strings = c("",".","NA"))
###check the data format
str(normadata)
###change the format of the several variables 
normadata$Entry=as.factor(normadata$Entry)
normadata$Rep=as.factor(normadata$Rep)
normadata$Year=as.factor(normadata$Year)
###check the data format
str(normadata)
###get the heritablity of all the traits using the function in the Function
source("Function/Heritability.R")
Herit <- Heritability(4, 7, normadata)
write.csv(Herit, file = "data/heritabilityflopc.csv", row.names = T, na = ".")