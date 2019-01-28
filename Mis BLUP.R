###function
###import the dataset
traits <- read.csv("Copy of alltraitsflowerday.csv", na.strings = c("", ".", "NA"))
str(traits)

traits <- traits[ ,c(2, 6:7, 15:25)]
traits$Entry= as.factor(traits$Entry)
traits$Rep = as.factor(traits$Rep)
traits$Year = as.factor(traits$Year)
str(traits)
library(lme4)
library(Matrix)
out_start=4
out_end=14
out_nvar=out_end-out_start+1
out_variable = colnames(traits)[out_start:out_end]
# set up a matrix to contain all BLUP values
out_beta <- matrix(0, nrow = length(levels(traits$Entry)),
                   ncol = out_nvar,
                   dimnames = list(levels(traits$Entry), 
                                   out_variable))
number=1
for (i in out_start:out_end){
  outcome = colnames(traits)[i]
  model <- lmer(get(outcome)~ (1|Entry)+ (1|Rep) + (1|Year) + (1|Entry:Year),
                na.action = na.exclude,
                data=traits)
  beta <- ranef(model)
  out_beta[rownames(beta$Entry), number] = beta$Entry[[1]]
  number = number + 1
}

# out_beta is the new output
