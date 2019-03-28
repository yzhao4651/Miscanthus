colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("trait",i-1,sep="")
  }  
  return(x)
}