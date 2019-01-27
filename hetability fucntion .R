
### write funcion 
####heritability
heritability <- function(x,y){
  library(lme4)
  library(Matrix)
  Entry= as.factor(y$Entry)
  REP = as.factor(y$Rep)
  x=as.numeric(x)
  Misvarcomp <- lmer(x~ (1|Entry)+ (1|Rep) + (1|Year), data=y)
  A=unlist(VarCorr(Misvarcomp))[[1]]
  B=attr(VarCorr(Misvarcomp), "sc")^2
  h2=A/(A+B/8)
  ####question: I used 8 as the number of Var_residual devided, but i am not sure, please let me know if it is not right. Thanks.
  return(h2)
}
###import the data 
qualdat <- read.csv("~/Documents/whole traits/Copy of alltraitsflowerday.csv")
heritability(qualdat$fday,qualdat)



