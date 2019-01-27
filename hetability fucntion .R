
### write funcion 
####heritability
heritability <- function(x, y){
  require(lme4)
  require(Matrix)
  Entry = as.factor(y$Entry)
  Rep = as.factor(y$Rep)
  Year = as.factor(y$Year)
  x = as.numeric(as.character(x))
  Misvarcomp <- lmer(x ~ (1|Entry) + (1|Rep) + (1|Year))
  A = unlist(VarCorr(Misvarcomp))[[1]]
  B = attr(VarCorr(Misvarcomp), "sc") ^ 2
  h2 = A/(A+B/8)
  ####question: I used 8 as the number of Var_residual devided, but i am not sure, please let me know if it is not right. Thanks.
  return(h2)
}
###import the data 
qualdat <- read.csv("Copy of alltraitsflowerday.csv")
heritability(qualdat$fday,qualdat)



