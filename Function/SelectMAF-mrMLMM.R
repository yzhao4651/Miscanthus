Select.MAF <- function(SNP){
  SNPmeans <- colMeans(SNP, na.rm = TRUE)
  #range(SNPmeans)
  #levels(as.factor(SNPmeans))
  # determine for which SNPs 0 is the major allele
  ZeroIsMajor <- SNPmeans < 1
  # set up an empty vector to indicate which SNPs to keep
  SNPsToKeep <- logical(length(ZeroIsMajor))
  # define a minimum number of accessions with minor allele
  minWithMinor <- 3
  # figure out which to keep based on how many individuals have minor allele
  SNPsToKeep[ZeroIsMajor] <- colSums(SNP[, ZeroIsMajor] == 1 | SNP[, ZeroIsMajor] == 2, na.rm = TRUE) >= minWithMinor
  SNPsToKeep[!ZeroIsMajor] <- colSums(SNP[, !ZeroIsMajor] == 1 | SNP[, !ZeroIsMajor] == 0, na.rm = TRUE) >= minWithMinor
  # this is how I would filter based on MAF
  MAF <- SNPmeans/2
  MAF[!ZeroIsMajor] <- 1 - MAF[!ZeroIsMajor] # flip the frequencies if allele 2 is the common one
  
  hist(MAF) # most near zero, max 0.5
  SNPsToKeep <- SNPsToKeep & MAF > 0.01 # update SNPsToKeep; 3 ind. with minor allele AND MAF > 0.01. 35,279 SNPs.
  # subset genotype matrix 
  ###using this mymat as new genotype for the next step. 
  subgenotranonly <- SNP[, SNPsToKeep]-1
  ###transform the data back again with SNP in the row and individual in the column
  subgenotranonly <- data.frame(t(subgenotranonly))
  #############changing the rowname to the first columne####################
  rn <- rownames(subgenotranonly)
  rownames(subgenotranonly) <- NULL
  subgenotranonly <- cbind(rn,subgenotranonly)
  return(subgenotranonly)
}