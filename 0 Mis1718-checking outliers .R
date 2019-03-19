####import the data####
library(readr)
##qualdat  <- read.csv("~/Documents/whole traits/trait1718.3.16.19.csv" , na.strings = c("",".","NA"))
qualdat  <- read.csv("data/trait1718.3.16.19.csv" , na.strings = c("",".","NA")) # on Yongli's computer
# mywd <- "." # on Lindsay's computer
#qualdat <- read.csv(file.path(mywd, "Copy of trait1718-4.csv"), na.strings = c("",".","NA"))
#check the data formate
str(qualdat)
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
####check the data formate
str(qualdat)
###select the data need for analysis 
qualdat <- qualdat[,c(3,5:6,4,7:27)]
####check the data formate
str(qualdat)
#### change several variables numeric format
# Note from Lindsay -- converting from integer to float is unnecessary, but should not cause problems
indx <- sapply(qualdat[,c(16:17,19:22)], is.integer)
qualdat[,c(16:17,19:22)][indx] <- lapply(qualdat[,c(16:17,19:22)][indx], function(x) as.numeric(as.character(x)))
###change several variables format
qualdat$GS <- as.numeric(as.character(qualdat$GS))
qualdat$Entry=as.factor(qualdat$Entry)
qualdat$Rep=as.factor(qualdat$Rep)
qualdat$Year=as.factor(qualdat$Year)

####elow will for checking outlies and get the new data without outlies
####elow will for checking outlies and get the new data without outlies
###below will for checking outlies and get the new data without outlies 
###elow will for checking outlies and get the new data without outlies
####check the data format
str(qualdat)
###seperate one data set to two data set accodring to the year 
qualdat.17 <- subset(qualdat,qualdat$Year==2017)
str(qualdat.17)
qualdat.18 <- subset(qualdat,qualdat$Year==2018)
str(qualdat.18)
###download the packages need for checking and removing the outliers 
install.packages("dplyr")
install.packages("tidyr")
install.packages("ruler")
install.packages("ggplot2")
library(dplyr)
library(tidyr)
library(ggplot2)
library(ruler)
####function with Z score
isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
}
####function with z-score with Median Absolute Deviation 
isnt_out_mad <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - median(x, na.rm = na.rm)) <= thres * mad(x, na.rm = na.rm)
}
### funtion with tukey method
isnt_out_tukey <- function(x, k = 1.5, na.rm = TRUE) {
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  (quar[1] - k * iqr <= x) & (x <= quar[2] + k * iqr)
}
####combining  three function to one 
isnt_out_funs <- funs(
  z = isnt_out_z,
  mad = isnt_out_mad,
  tukey = isnt_out_tukey
)
###check the outliers for traits in 2017 and also output in one file without outliers
###check the outliers for traits in 2017 and also output in one file without outliers
compute_group_non_outliers <- . %>%
  # Compute per group mean values of columns
  group_by(Entry) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup() %>%
  # Detect outliers among groups
  mutate_if(is.numeric, isnt_out_tukey) %>%
  # Remove unnecessary columns
  select_if(Negate(is.numeric))
outliers.17 <- qualdat.17 %>% compute_group_non_outliers()

###check the outliers for traits in 2018 and also output in one file without outliers
###check the outliers for traits in 2018 and also output in one file without outliers
outliers.18 <- qualdat.18 %>% compute_group_non_outliers()


### Checking which one is the outliner number 
### Checking which one is the outliner number 
group_packs_isnt_out <- group_packs(
  # Non-outliers based on grouping
  group = compute_group_non_outliers,
  .group_vars = "Entry"
)
# Don't remove obeyers to compute total number of applied rules
full_report <- qualdat.18 %>%
  expose(group_packs_isnt_out,
         .remove_obeyers = FALSE) %>%
  get_report()

used_rules <- full_report %>%
  distinct(pack, rule)

breaker_report <- full_report %>%
  filter(!(value %in% TRUE))

group_breakers <- breaker_report %>%
  # Filter group packs
  filter(pack == "Entry") %>%
  # Expand rows by matching group with its rows
  select(-id) %>%
  left_join(
    y = qualdat.18 %>% transmute(var = Entry, id = 1:n()),
    by = "var"
  ) %>%
  select(pack, rule, var, id, value)

outliers <- bind_rows(
  breaker_report %>% filter(pack != "Entry"),
  group_breakers
) %>%
  select(pack, rule, id)

# Not all group based definitions resulted with outliers
outliers %>%
  count(pack, rule) %>%
  filter(pack == "Entry") %>%
  print(n = Inf)

## Note from Lindsay -- the last command seems to indicate that there were no
## outliers, so the code below is unnecessary.
## Additionally, I am not sufficiently familiar with Tidyverse to understand
## all of the code in this section.

###check the dataset is the same to the origianl one or not 
all.equal(qualdat.18.no.outlies,qualdat.18)

### combination both data set together from vertial direction 
# vertical merge
install.packages("lessR")
library(lessR)
qualdat <- Merge(qualdat.17, qualdat.18)
####check the data formate
str(qualdat)
###save this data set 
save(qualdat,file="qualdat.RData")
####check the data formate
str(qualdat)



