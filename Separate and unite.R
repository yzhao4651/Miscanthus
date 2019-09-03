###question:
##I did the gene search and then anonation
###qestion1: many SNPs can not find their genes from Sorghum bicolor, Is it right? total 895 SNPs, I found only 606 has hit. 
###question2: there are some hit with GO terms, and also Best.hit.arabi.name like AT1G01540.2, 
###I used GO number to get the GOterm, but found many function are related to animals. 
###and then I doneload one txt file from Tair web (I named it ATHgenes and also attach it from email, can not push to git, becasue of large file)
###I am trying to get the more information about their function using Best.hit.arabi.name like AT1G01540.2,
###and then i found that file with locus.name like AT1G01540 without ".2", and then the same AT1G01540 matched so many GO number, it is hard to choose which one?
###but one column called "object.name" has names like AT1G01540.2, I tried to used that one to search, but the GO term from that file is not matched from original anotation? 
####i mean from the "Sbicolor_454_v3_1_1_annotation_info.txt", would you please give me some advice? 
###by the way, i did get several genes match with flowering time, i used one webs called http://www.flor-id.org, but not used best. ara.hit.name, but the symbol.
###below is codes question
###below: some codes i am trying to separate columns and them unite across different columns together, 
library(tidyr)
TO <- read.csv("data/Total.SNP.csv")
TO1 <- TO[2:16]
str(TO1)
###seperate mutiple column,
for(name in names(TO1)) {
  TO1 <- separate(TO1, name, into = paste0(name, '_col', 1:17),";")
}
TO <- cbind(TO[1],TO1)
str(TO)
###this one does not work: 
### I am trying to unite X1_col1, X2_col1, X3_col1, X4_col1, X5_col1, .....X15_col1 to one column with ";" separated. 
###and then X1_col2, X2_col2, X3_col2, X4_col2, X5_col2, .....X15_col2. and until X1_col17, X2_col17, X3_col17, X4_col17, X5_col17, .....X15_col17.
###but does not work using below the function, would you please give some advice? thanks 
cols <- unique(substr(names(TO)[-1], 1, 256)) ## this gives a result identical to names(TO)[-1], so I am not sure what you are trying to accomplish.
str(cols)
library(dplyr)
TO_New <- do.call(cbind,lapply(cols, function(x){unite_(TO[2:256], x, grep(x, names(TO[2:256]), value = TRUE), 
                                                        sep = ';', remove = TRUE) %>% select_(x)})
)

str(TO_New)

## I am not very experienced with dplyr, but here is how I would do it in the R base.
placeholder <- rep(NA_character_, nrow(TO))
TO_New <- data.frame(SNP = as.character(TO[[1]]),
                     col1 = placeholder, col2 = placeholder, col3 = placeholder, col4 = placeholder,
                     col5 = placeholder, col6 = placeholder, col7 = placeholder, col8 = placeholder,
                     col9 = placeholder, col10 = placeholder, col11 = placeholder, col12 = placeholder,
                     col13 = placeholder, col14 = placeholder, col15 = placeholder, col16 = placeholder, col17 = placeholder)
for(col in 1:17){
  thesecols <- paste("X", 1:15, "_col", col, sep = "")
  TO_New[[paste("col", col, sep = "")]] <- do.call(paste, c(as.list(TO[,thesecols]), list(sep = ";")))
}

