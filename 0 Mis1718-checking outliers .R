####import the data####
source("0 Qualdat import.R")
qualdat <- read_qualdat("data/trait1718.3.16.19.csv")

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
outliers_count <- outliers %>%
  count(pack, rule) %>%
  filter(pack == "Entry.Year") %>%
  print(n = Inf)

## Note from Lindsay --  I am not sufficiently familiar with Tidyverse to 
## understand all of the code in this section.

### Visualize outliers using grouping method
plot_outliers <- function(trait, .qualdat = qualdat, .outliers = outliers){
  # get rows that are outliers for this trait
  outrows <- .outliers$id[.outliers$rule == trait]
  mycol <- rep("black", nrow(.qualdat))
  mycol[outrows] <- "red"
  mypch <- ifelse(.qualdat$Year == 2017, 1, 0)
  plot(.qualdat[[trait]], col = mycol, pch = mypch, ylab = trait,
       xlab = "Row of qualdat")
}

plot_outliers("HD_1")
plot_outliers("FD_1")
plot_outliers("HD_50.")
plot_outliers("FD_50.")
plot_outliers("FNMain")
plot_outliers("FNsmall")
plot_outliers("CmDW_g")
plot_outliers("TFN.")
plot_outliers("Cml_cm")
plot_outliers("CmD_BI_mm")
plot_outliers("CmD_LI_mm")
plot_outliers("CmN.")
plot_outliers("Bcirc_cm")
plot_outliers("Yld_kg")
plot_outliers("SDW_kg")
plot_outliers("Lg")
plot_outliers("GS")
plot_outliers("FD")
plot_outliers("CCirc_cm")
plot_outliers("SRD")
plot_outliers("ADD")

## Note from Lindsay: There are clearly some outliers for some traits, but I
## am not sure that this grouping method is the best for every trait, probably
## due to the small replication within Entry*Year.

### Visualize outliers without grouping

plot_outliers2 <- function(trait, .qualdat = qualdat){
  outliers <- !isnt_out_tukey(.qualdat[[trait]])
  outliers[is.na(outliers)] <- FALSE
  mycol <- rep("black", nrow(.qualdat))
  mycol[outliers] <- "red"
  mypch <- mypch <- ifelse(.qualdat$Year == 2017, 1, 0)
  plot(.qualdat[[trait]], col = mycol, pch = mypch, ylab = trait,
       xlab = "Row of qualdat")
}

plot_outliers2("HD_1")     # Some clear outliers identified
plot_outliers2("FD_1")     # No outliers
plot_outliers2("HD_50.")   # Some clear outliers
plot_outliers2("FD_50.")   # No outliers
plot_outliers2("FNMain")   # Should probably be transformed before looking for outliers
plot_outliers2("FNsmall")  # Should probably be transformed before looking for outliers
plot_outliers2("CmDW_g")   # Should probably be transformed before looking for outliers
plot_outliers2("TFN.")     # Should probably be transformed before looking for outliers
plot_outliers2("Cml_cm")   # Some outliers identified, but don't seem far outside of distribution
plot_outliers2("CmD_BI_mm")# Some outliers identified, but don't seem far outside of distribution
plot_outliers2("CmD_LI_mm")# Some outliers identified, but don't seem far outside of distribution
plot_outliers2("CmN.")     # Some outliers identified, but don't seem far outside of distribution
plot_outliers2("Bcirc_cm") # Should probably be transformed before looking for outliers
plot_outliers2("Yld_kg")   # Should probably be transformed before looking for outliers
plot_outliers2("SDW_kg")   # No outliers
plot_outliers2("Lg")       # Too qualitative to detect outliers
plot_outliers2("GS")       # Too qualitative to detect outliers
plot_outliers2("FD")       # No outliers
plot_outliers2("CCirc_cm") # Should probably be transformed before looking for outliers
plot_outliers2("SRD")      # One clear outlier
plot_outliers2("ADD")      # No outliers

### Note: once you have determined a way to detect outliers for each trait, you will need
### to modify qualdat before saving it.  Your previous code did not modify anything.

# dummy example code
qualdat$FD_50.[is_an_outlier(qualdat$FD_50.)] <- NA

###save this data set 
save(qualdat,file="qualdat.RData")

