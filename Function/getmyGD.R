getmyGD <- function(myGM,myGD){
  n <- data.frame(t(data.frame(myGD.clum.107.3707,row.names = 1)))
  myGM2 <- myGM[match(dimnames(n)[[1]],myGM$Name, nomatch=0),]
  snporder <- order(myGM2$Chromosome, myGM2$Position)
  myGM2 <- myGM2[snporder,]
  myGD <- myGD[,c(1, snporder + 1)]
  return(myGD)
}