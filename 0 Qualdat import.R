# Script to generate the `qualdat` object for analysis

read_qualdat <- function(file){
  qualdat  <- read.csv(file , na.strings = c("",".","NA"))
  ###change character of Growth Stage to number N=1,B=2,F=3,P=4
  levels(qualdat$GS)[levels(qualdat$GS)=="N"] <- "1"
  levels(qualdat$GS)[levels(qualdat$GS)=="n"] <- "1"
  levels(qualdat$GS)[levels(qualdat$GS)=="B"] <- "2"
  levels(qualdat$GS)[levels(qualdat$GS)=="b"] <- "2"
  levels(qualdat$GS)[levels(qualdat$GS)=="F"] <- "3"
  levels(qualdat$GS)[levels(qualdat$GS)=="p"] <- "4"
  levels(qualdat$GS)[levels(qualdat$GS)=="P"] <- "4"
  ###also trying to get the week number for flowering traits
  require(lubridate)
  ###get the week number 
  qualdat$HW_1 <- lubridate::week(mdy(qualdat$HD_1))
  qualdat$FW_1 = lubridate::week(mdy(qualdat$FD_1))
  qualdat$HW_50 = lubridate::week(mdy(qualdat$HD_50.))
  qualdat$FW_50 = lubridate::week(mdy(qualdat$FD_50.))
  ###using this function to get the group of week (group methods: 23 and 24 change into 12, 
  ###and the 25, 26 to 13,27 and 28 to 14 ......)
  group_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("G",colnames(data[i]))] <- (data[[i]]+1) %/% 2
    }
    return(data)
  }
  qualdat <- group_function(qualdat,31,34)
  ##change the data in to month 
  qualdat$HM_1 <- lubridate::month(mdy(qualdat$HD_1))
  qualdat$FM_1 = lubridate::month(mdy(qualdat$FD_1))
  qualdat$HM_50 = lubridate::month(mdy(qualdat$HD_50.))
  qualdat$FM_50 = lubridate::month(mdy(qualdat$FD_50.))
  ####calculate days 
  qualdat$SRD <- as.numeric(as.Date(qualdat$SRD,format = "%m/%d/%Y")-as.Date(qualdat$datest2,format = "%m/%d/%Y"))
  qualdat$ADD <-as.numeric(as.Date(qualdat$ADD,format = "%m/%d/%Y")-as.Date(qualdat$datest2,format = "%m/%d/%Y"))
  qualdat$HD_1 <- as.numeric(as.Date(qualdat$HD_1,format = "%m/%d/%Y")-as.Date(qualdat$datest1,format = "%m/%d/%Y"))
  qualdat$FD_1 <-as.numeric(as.Date(qualdat$FD_1,format = "%m/%d/%Y")-as.Date(qualdat$datest1,format = "%m/%d/%Y"))
  qualdat$HD_50. <- as.numeric(as.Date(qualdat$HD_50.,format = "%m/%d/%Y")-as.Date(qualdat$datest1,format = "%m/%d/%Y"))
  qualdat$FD_50. <- as.numeric(as.Date(qualdat$FD_50.,format = "%m/%d/%Y")-as.Date(qualdat$datest1,format = "%m/%d/%Y"))
  ###checking the format of the data
  #str(qualdat)
  ###select the data need for analysis 
  qualdat <- qualdat[,c(2:27,31:42)]
  ###checking the format of the data
 #str(qualdat)
  ###change several variables format
  qualdat$GS <- as.numeric(as.character(qualdat$GS))
  qualdat$Entry=as.factor(qualdat$Entry)
  qualdat$Rep=as.factor(qualdat$Rep)
  qualdat$Year=as.factor(qualdat$Year)
  #str(qualdat)
  #### change several variables numeric format
  # Note from Lindsay -- converting from integer to float is unnecessary, but should not cause problems
  indx <- sapply(qualdat[,c(5,14,19,21:24)], is.integer)
  qualdat[,c(5,14,19,21:24)][indx] <- lapply(qualdat[,c(5,14,19,21:24)][indx], function(x) as.numeric(as.character(x)))
  ###check the data format 
  #str(qualdat)
  return(qualdat)
}
