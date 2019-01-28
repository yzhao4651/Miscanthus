###function
###import the dataset
traits <- read.csv("Copy of alltraitsflowerday.csv")
str(traits)
###change factor into numeric 
indx <- sapply(traits[,c(15:25)], is.factor)
traits[,c(15:25)][indx] <- lapply(traits[,c(15:25)][indx], function(x) as.numeric(as.character(x)))
traits <- traits[ ,c(2, 6:7, 15:25)]
traits$Entry= as.factor(traits$Entry)
traits$Rep = as.factor(traits$Rep)
traits$Year = as.factor(traits$Year)
str(traits)
library(lme4)
library(Matrix)
out_start=15
out_end=25
out_nvar=out_end-out_start+1
out_variable=rep(NA, out_nvar)
number=1
for (i in 15:25){
  outcome = colnames(traits)[i]
    model <- lmer(get(outcome)~ (1|Entry)+ (1|Rep) + (1|Year),
                  na.action = na.exclude,
                  data=traits)
    beta <- ranef(model)
    out_beta[number] = as.numeric(beta[2])
    out_variable[number] = outcome
    number = number + 1
}
outcome = data.frame(out_variable, out_beta)



