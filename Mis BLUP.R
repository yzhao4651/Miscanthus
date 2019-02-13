###import the data
###import the data
#normadata <- read.csv("~/Documents/whole traits/traits1718normalited1.csv",na.strings = c("",".","NA"))
normadata <- read.csv("traits1718normalited1.csv",na.strings = c("",".","NA"))
###check the data format
str(normadata)
### this one help me to get which variables has more missing values and help to set up the condition 
out_start <- 8
out_end <- 28
for (i in out_start:out_end){
  t <-cbind(colnames(normadata[i]),sum(is.na(normadata[normadata$Year==2017 | normadata$Year==2018,][i])))
  print(t)
}
    
###change the format of the several variables 
normadata$GS <- as.numeric(as.character(normadata$GS))
normadata$SRD <- as.numeric(as.character(normadata$SRD))
normadata$Entry=as.factor(normadata$Entry)
normadata$Rep=as.factor(normadata$Rep)
normadata$Year=as.factor(normadata$Year)
###check the data format
str(normadata)
###function for the ranef of all of traits 
ranefvalue <- function(out_start,out_end,y){
  require(lme4)
  require(Matrix)
  Entry = as.factor(y$Entry)
  Rep = as.factor(y$Rep)
  Year = as.factor(y$Year)
  out_nvar=out_end-out_start+1
  out_variable = colnames(y[out_start:out_end])
  out_beta <- matrix(0, nrow = length(levels(y$Entry)),
                     ncol = out_nvar,
                     dimnames = list(levels(y$Entry), 
                                     out_variable))
  number=1
  for (i in out_start:out_end){
    if(sum(is.na(y[y$Year==2017 | y$Year==2018,][i])) >=900 ){
      outcome = colnames(y)[i]
      model <- lmer(get(outcome)~ (1|Entry)+ (1|Rep),
                    na.action = na.exclude,
                    data=y)
      beta <- ranef(model)
    } else {
  outcome = colnames(y)[i]
  model <- lmer(get(outcome)~ (1|Entry)+ (1|Rep) + (1|Year) + (1|Entry:Year),
                na.action = na.exclude,
                data=y)
  beta <- ranef(model)
    }
    out_beta[rownames(beta$Entry), number] = beta$Entry[[1]]
    number = number + 1
  }
 return(out_beta)
}

ranefvalueall <- ranefvalue(8,28,normadata)
