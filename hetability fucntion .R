
### write funcion 
####heritability
heritability <- function(x, y){
  require(lme4)
  require(Matrix)
  
  Entry = as.factor(y$Entry)
  Rep = as.factor(y$Rep)
  Year = as.factor(y$Year)
  
  Misvarcomp <- lmer(x ~ (1|Entry) + (1|Rep) + (1|Year) + (1|Entry:Year))
  A = unlist(VarCorr(Misvarcomp))[["Entry"]]         # variance attributable to entry
  B = attr(VarCorr(Misvarcomp), "sc") ^ 2            # residual variance
  C = unlist(VarCorr(Misvarcomp))[["Entry:Year"]]    # variance attributable to entry by year
  lvlsB <- length(levels(Rep:Year))       # get the divisor for residual (8)
  lvlsC <- length(levels(Year))           # get the divisor for entry by year (2)
  
  h2 = A/(A + B/lvlsB + C/lvlsC)
  ####question: I used 8 as the number of Var_residual devided, but i am not sure, please let me know if it is not right. Thanks.
  #### yes, 8 is correct. -LVC
  return(h2)
}
# ###import the data 
# qualdat <- read.csv("Copy of alltraitsflowerday.csv", na.strings = ".")
# heritability(qualdat$fday,qualdat)
# 
# # for debugging
# x <- qualdat$fday
# y <- qualdat
