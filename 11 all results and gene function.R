##Merge all Results 

multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
}
full_data = multmerge("allresults")
str(full_data)
write.csv(full_data,file="allresults/full_data.csv")

###using blast to get the gene functions
##install packages
#install.packages("BiocManager")
#BiocManager::install("rtracklayer")
#BiocManager::install("GenomicFeatures")
###import GFF3 with  makeTxDbFromGFF function in the package ("GenomicFeatures")
library(rtracklayer)
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("Sbicolor_454_v3.1.1.gene.gff3",
                        format = "gff3",
                        dataSource = "Phytosome 12",
                        organism = "Sorghum bicolor")
seqinfo(txdb)
#To save the txdb database for later and avoid having to recreate it every time we use it, we can use saveDb() and, later, loadDb()
BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
saveDb(txdb, 'txdb.Sbicolor_454_v3.1.1')
#Once the database is stored to disk, it can be reloaded later, almost instantly.
library(AnnotationDbi)
txdb = loadDb(file = 'txdb.Sbicolor_454_v3.1.1')
genes(txdb)
seqinfo(txdb)
###install.packags library(IRanges)
#BiocManager::install("IRanges")
### fucntion to get the neargenes 

Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)
    saved_genes[[i]] <- mygenes$tx_name
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        Gene = unlist(saved_genes))
  return(gene_df)
}

gene_df <- Findneargenes(full_data, 1e4)
str(gene_df)

###but i did not get any information from gene_df?


### I did not doing the rest of them 
write.csv(gene_df, file = "test_genes.csv", row.names = FALSE)

annotations <- read.delim("Osativa_323_v7.0.annotation_info.txt",
                          stringsAsFactors = FALSE)
head(annotations)

myrows <- match(gene_df$Gene, annotations$transcriptName)
write.csv(annotations[myrows,], file = "test_annotations.csv")
