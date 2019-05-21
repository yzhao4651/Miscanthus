# skeleton of a function to do what '10 multipl imageflow.R' does

source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/R/CMplot.r")

# copy of Yongli's function for correcting p-values
adj_P_function <- function(data, start_var, end_var){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}

## Some considerations as you edit the main function:
## You might standardize the naming and locations of your results files (e.g. for
## programs other than GAPIT) to make it easier to access them programmatically.
## If the function needs to behave differently for the single locus vs. multi
## locus methods, use the `multi` argument to control that behavior.

# Example use
# for(tr in c("HD_1", "FD_1", "HD_50", "FD_50", "CM_N", "CmDW_g")){
#   multiple_plots(trait = tr)
# }

# Main function
multiple_plots <- function(trait = "HD_1", methods = c("GLM", "MLM", "CMLM", "SUPER"),
                           files = c("Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.GLM.TRAIT.GWAS.Results.csv",
                                     "Result GAPIT1/GLMMLMMLMMFarmCPUCV/GAPIT.MLM.TRAIT.GWAS.Results.csv",
                                     "Result GAPIT1/ECMMLCV/GAPIT.CMLM.TRAIT.GWAS.Results.csv",
                                     "Result GAPIT1/MLMSUPPER/GAPIT.SUPER.TRAIT.GWAS.Results.csv"),
                           colors = c("blue", "orange", "cyan","magenta"),
                           imagedir = file.path("Allimages", trait),
                           multi = FALSE){
  nmeth <- length(methods)
  if(length(files) != nmeth){
    stop("Need one file for each method.")
  }
  if(length(colors) != nmeth){
    stop("Need one color for each method.")
  }
  
  origwd <- getwd() # current working directory
  
  # Wherever there is "TRAIT" in the file name, substitute in the name of the trait we want.
  files <- gsub("TRAIT", trait, files)
  
  # read in files and join into one table
  if(multi){
    ### Here, insert any code needed for the multi-locus methods
  } else {
    for(i in 1:nmeth){
      newcolname <- paste(methods[i], trait, sep = ".") # e.g. "GLM.HD_1"
      if(i == 1){
        allp <- read.csv(files[i], stringsAsFactors = FALSE)[,1:4]
        colnames(allp)[colnames(allp)=="P.value"] <- newcolname
      } else {
        nextp <- read.csv(files[i], stringsAsFactors = FALSE)[,c(1,4)]
        allp[[newcolname]] <- nextp[[2]][match(allp$SNP, nextp$SNP)]
      }
    }
  }
  
  # adjust p-values (`allp` was called `plot` in the script but that is a poor name)
  allp <- adj_P_function(allp, 4, ncol(allp))
  
  ## Make plots
  setwd(imagedir)
  # QQ plot
  CMplot(allp, plot.type = "q", col = colors, threshold=1.13e-5,
         signal.pch=19,cex=0.5,signal.cex=1.0, signal.col="blue",conf.int.col="grey",box=FALSE,multracks=
           TRUE,file="jpg",memo="",dpi=300,
         file.output = TRUE, verbose=TRUE)
  
  ### Etc.... add in other calls to CMplot or anything else you want to happen
  
  # return to the original working directory
  setwd(origwd)
}
