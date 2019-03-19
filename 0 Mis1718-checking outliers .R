####import the data####
source("0 Qualdat import.R")

####elow will for checking outlies and get the new data without outlies
####elow will for checking outlies and get the new data without outlies
###below will for checking outlies and get the new data without outlies 
###elow will for checking outlies and get the new data without outlies
####check the data format
str(qualdat)

###download the packages need for checking and removing the outliers 
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("ruler")
#install.packages("ggplot2")
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

### make a grouping based on Entry and Year, since we want to treat years separately
qualdatEY <- unite(qualdat, col = "Entry.Year", Entry, Year)
str(qualdatEY)
qualdatEY$Entry.Year

###check the outliers for traits both years
compute_group_non_outliers <- . %>%
  # Compute per group mean values of columns
  group_by(Entry.Year) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup() %>%
  # Detect outliers among groups
  mutate_if(is.numeric, isnt_out_tukey) %>%
  # Remove unnecessary columns
  select_if(Negate(is.numeric))

my_outliers <- qualdatEY %>% compute_group_non_outliers()
my_outliers
my_outliers[11:20,]

### Checking which one is the outliner number 
### Checking which one is the outliner number 
group_packs_isnt_out <- group_packs(
  # Non-outliers based on grouping
  Entry.Year = compute_group_non_outliers,
  .group_vars = "Entry.Year"
)
# Don't remove obeyers to compute total number of applied rules
full_report <- qualdatEY %>%
  expose(group_packs_isnt_out,
         .remove_obeyers = FALSE) %>%
  get_report()

used_rules <- full_report %>%
  distinct(pack, rule)

breaker_report <- full_report %>%
  filter(!(value %in% TRUE))

group_breakers <- breaker_report %>%
  # Filter group packs
  filter(pack == "Entry.Year") %>%
  # Expand rows by matching group with its rows
  select(-id) %>%
  left_join(
    y = qualdatEY %>% transmute(var = Entry.Year, id = 1:n()),
    by = "var"
  ) %>%
  select(pack, rule, var, id, value)

outliers <- bind_rows(
  breaker_report %>% filter(pack != "Entry.Year"),
  group_breakers
) %>%
  select(pack, rule, id)

# Not all group based definitions resulted with outliers
my_outliers <- outliers %>%
  count(pack, rule) %>%
  filter(pack == "Entry.Year") %>%
  print(n = Inf)

## Note from Lindsay --  I am not sufficiently familiar with Tidyverse to 
## understand all of the code in this section.

### combination both data set together from vertial direction 
# vertical merge
#install.packages("lessR")
library(lessR)
qualdat <- Merge(qualdat.17, qualdat.18)
####check the data formate
str(qualdat)
###save this data set 
save(qualdat,file="qualdat.RData")
####check the data formate
str(qualdat)



