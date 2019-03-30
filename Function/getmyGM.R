getmyGM <- function(myGM,myGD){
  n <- data.frame(t(data.frame(myGD,row.names = 1)))
  myGM2 <- myGM[match(dimnames(n)[[1]],myGM$Name, nomatch=0),]
  snporder <- order(myGM2$Chromosome, myGM2$Position)
  myGM2 <- myGM2[snporder,]
  myGD <- myGD[,c(1, snporder + 1)]
  return(myGM2)
}