subsetNONmisingSNP.GAPIT <- function(datacomb,allblup){
  datacomb <- data.frame(datacomb)
  Taxa <- rownames(datacomb)
  rownames(datacomb) <- NULL
  datacomb <- cbind(Taxa,datacomb)
  ##############step3 select matched SNP 
  ##############step3 select matched SNP 
  ### checking Taxa if match in SNP and in phenotype data 
  allblup$Taxa %in% datacomb$Taxa
  ### submet the SNP data set match with the phenotype with total missing value
  SNP <- datacomb[match(allblup$Taxa, datacomb$Taxa,nomatch=0),]
  ###change colomn to row for datacomb2
  SNP <- data.frame(SNP,row.names=1)
  return(SNP)
}