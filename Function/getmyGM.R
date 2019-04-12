getmyGM <- function(myGM,myGD){
  n <- data.frame(t(data.frame(myGD)))
  myGM2 <- myGM[match(dimnames(n)[[1]],myGM$Name, nomatch=0),]
  snporder <- order(myGM2$Chromosome, myGM2$Position)
  myGM2 <- myGM2[snporder,]
  myGD <- myGD[,snporder]
  return(myGM2)
}