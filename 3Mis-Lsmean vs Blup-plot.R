####This codes for plotting the Lsmean anD BLUP of all traits
####import the data
ranefvalueall<- read.csv("data/ranefvalueall.csv",na.strings = c("",".","NA"))
###check data format
str(ranefvalueall)
###rename of the column name 
colnames(ranefvalueall)[colnames(ranefvalueall)=="X"] <- "Entry"
###check data format
str(ranefvalueall)
####firstly using Q-Q plot and histogram to check if the BLUP is normally distribution
####qq plot and histogram ####approximately normally distribution 
####qq plot and histogram 
source("Function/histqq.R")
pdf(paste("Q-Q normal plot and Histogram", 2 ,".pdf",sep=""))
histqq(ranefvalueall[-1])

####plot BLUP VS LineBlup
###firstly load the data of qualdat
load("qualdat.RData")
####check the format of dataset
str(qualdat)
####trying to write function to plot 
####Compare BLUP to line averages on a scatterplot
source("Function/bluplinevslsmean.R")
pdf(paste("Linemean VS LineBlup", 1 ,".pdf",sep="")) 
bluplinevslsmean(ranefvalueall[-1],qualdat)






