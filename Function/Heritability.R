####function heritability to get heritability for each variable 
####out_start: the first variable, 
####out_end: the last variable
####heritablity for all of the varialbes 

### this one help me to get which variables has more missing values and help to set up 
###the condition should larger then 900, or other number 
### in the function Heritability 
#out_start <- 4
#out_end <- 25
#for (i in out_start:out_end){
  #t <-cbind(colnames(normadata[i]),sum(is.na(normadata[normadata$Year==2017 | normadata$Year==2018,][i])))
  #print(t)
#}

Heritability <- function(out_start,out_end,y){
  require(lme4)
  require(Matrix)
  Entry = as.factor(y$Entry)
  Rep = as.factor(y$Rep)
  Year = as.factor(y$Year)
  out_nvar=out_end-out_start+1
  out_variable = colnames(y[out_start:out_end])
  outcome <- matrix(NA_real_, nrow=ncol(y[out_start:out_end]),
                    ncol =1, dimnames = list(out_variable, "Heritability"))
  number=1
  for (i in out_start:out_end){
    if(sum(is.na(y[y$Year==2017 | y$Year==2018,][i])) >= 900){
      Misvarcomp<- lmer(y[,i] ~ (1|Entry) + (1|Rep))
      A= unlist(VarCorr(Misvarcomp))[["Entry"]]         # variance attributable to entry
      B= attr(VarCorr(Misvarcomp), "sc") ^ 2            # residual variance
      lvlsB <- length(levels(Rep))       # get the divisor for residual (8)
      h2 = A/(A+ B/lvlsB)
    } else {
    Misvarcomp<- lmer(y[,i] ~ (1|Entry) + (1|Rep) + (1|Year) + (1|Entry:Year))
    A= unlist(VarCorr(Misvarcomp))[["Entry"]]         # variance attributable to entry
    B= attr(VarCorr(Misvarcomp), "sc") ^ 2            # residual variance
    C = unlist(VarCorr(Misvarcomp))[["Entry:Year"]]    # variance attributable to entry by year
    lvlsB <- length(levels(Rep:Year))       # get the divisor for residual (8)
    lvlsC <- length(levels(Year))           # get the divisor for entry by year (2)
    h2 = A/(A+ B/lvlsB + C/lvlsC)
    }
    outcome[number] <-h2 
    number=number+1
 }
  return(outcome)
}

## Note from Lindsay: Right now you have two identically named functions, with
## very similar code.  If you have to fix something in one function, you have 
## to fix it in both, which leaves a lot of room for error in the future. It
## would be better to consolidate this into one function, perhaps with a 
## logical argument to control behavior, like

# heritability <- function(out_start,out_end,y, include_year = TRUE){
# ...
# if(!include_year || sum(is.na(y[y$Year==2017 | y$Year==2018,][i])) >= 900){ # see addtl edit below
# ...

## I think previously we loaded this function into another script using `source`.
## Since you have added actual data analysis to this script, in addition to the
## function definition(s), you probably do not want to load it with `source` 
## any more.

## The code `y$Year==2017 | y$Year==2018` seems unnecessary, since it evaluates
## to TRUE for every line of the data frame.  Are you instead trying to find
## traits that were only measured in one year?  If so, maybe do:

# all(is.na(y[y$Year == "2017", i])) || all(is.na(y[y$Year == "2018", i])) 
