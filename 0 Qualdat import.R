# Script to generate the `qualdat` object for analysis

qualdat  <- read.csv("data/trait1718.3.16.19.csv" , na.strings = c("",".","NA"))

###change character of Growth Stage to number N=1,B=2,F=3,P=4
levels(qualdat$GS)[levels(qualdat$GS)=="N"] <- "1"
levels(qualdat$GS)[levels(qualdat$GS)=="n"] <- "1"
levels(qualdat$GS)[levels(qualdat$GS)=="B"] <- "2"
levels(qualdat$GS)[levels(qualdat$GS)=="b"] <- "2"
levels(qualdat$GS)[levels(qualdat$GS)=="F"] <- "3"
levels(qualdat$GS)[levels(qualdat$GS)=="p"] <- "4"
levels(qualdat$GS)[levels(qualdat$GS)=="P"] <- "4"
####calculate days 
qualdat$SRD <- as.numeric(as.Date(qualdat$SRD,format = "%m/%d/%Y")-as.Date(qualdat$datest2,format = "%m/%d/%Y"))
qualdat$ADD <-as.numeric(as.Date(qualdat$ADD,format = "%m/%d/%Y")-as.Date(qualdat$datest2,format = "%m/%d/%Y"))
qualdat$HD_1 <- as.numeric(as.Date(qualdat$HD_1,format = "%m/%d/%Y")-as.Date(qualdat$datest1,format = "%m/%d/%Y"))
qualdat$FD_1 <-as.numeric(as.Date(qualdat$FD_1,format = "%m/%d/%Y")-as.Date(qualdat$datest1,format = "%m/%d/%Y"))
qualdat$HD_50. <- as.numeric(as.Date(qualdat$HD_50.,format = "%m/%d/%Y")-as.Date(qualdat$datest1,format = "%m/%d/%Y"))
qualdat$FD_50. <- as.numeric(as.Date(qualdat$FD_50.,format = "%m/%d/%Y")-as.Date(qualdat$datest1,format = "%m/%d/%Y"))

###select the data need for analysis 
qualdat <- qualdat[,c(3,5:6,4,7:27)]

#### change several variables numeric format
# Note from Lindsay -- converting from integer to float is unnecessary, but should not cause problems
indx <- sapply(qualdat[,c(16:17,19:22)], is.integer)
qualdat[,c(16:17,19:22)][indx] <- lapply(qualdat[,c(16:17,19:22)][indx], function(x) as.numeric(as.character(x)))
###change several variables format
qualdat$GS <- as.numeric(as.character(qualdat$GS))
qualdat$Entry=as.factor(qualdat$Entry)
qualdat$Rep=as.factor(qualdat$Rep)
qualdat$Year=as.factor(qualdat$Year)

rm(indx)