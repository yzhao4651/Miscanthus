##Merge all Results 

multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
  Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
}
full_data = multmerge("allresults")
full_data <- full_data[full_data$P.value != 0,]
str(full_data)
write.csv(full_data,file="data/full_data.csv")
str(full_data)
full_data.FarCPU <- subset(full_data,full_data$Method =="FarmCPU" & full_data$P.value >= 0.05)
tail(full_data.FarCPU$Trait.name)

max(full_data$P.value)

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

full_data <- read.csv("data/full_data.csv",row.names=1)
str(full_data)
gene_df <- Findneargenes(full_data, 1e4)
write.csv(gene_df, file = "data/test_genes.csv",row.names = FALSE)
str(gene_df)
head(gene_df)

###
### get the annotations for each genes 
annotations <- read.delim("Sbicolor_454_v3_1_1_annotation_info.txt",
                          stringsAsFactors = FALSE)
head(annotations)
str(annotations)
myrows <- match(gene_df$Gene, annotations$transcriptName)
allhitsgenes <- annotations[myrows,]
colnames(allhitsgenes)[colnames(allhitsgenes)=="transcriptName"] <- "Gene"
str(allhitsgenes)
head(allhitsgenes)

allSNPSGenes.2 <- plyr::join_all(list(allhitsgenes,gene_df), by="Gene")
allSNPSGenes.trt.mtd <- plyr::join_all(list(allSNPSGenes.2,full_data[-7]), by="SNP")
str(allSNPSGenes.trt.mtd)

allSNPSGenes.trt.mtd <- allSNPSGenes.trt.mtd[!duplicated(allSNPSGenes.trt.mtd[,c("Trait.name","Method","Chromosome","Position","P.value",
                                                                                 "SNP","locusName","arabi.symbol","arabi.defline","rice.defline","Pfam")]),]
str(allSNPSGenes.trt.mtd)

head(allSNPSGenes.trt.mtd)

#allSNPSGenes.trt.mtd$GO[allSNPSGenes.trt.mtd$GO==""]<- "AT"

library(dplyr)
library(tidyr)
allSNPSGenes.trt.mtd_new <- allSNPSGenes.trt.mtd %>% tidyr::separate(GO, c("GO1","GO2","GO3","GO4","GO5","GO6","GO7","GO8","GO9"), ",")
head(allSNPSGenes.trt.mtd_new)
str(allSNPSGenes.trt.mtd_new)
levels(as.factor(allSNPSGenes.trt.mtd_new$GO10))
allSNPSGenes.trt.mtd_new.2 <- allSNPSGenes.trt.mtd_new[,-c(11:18)]
colnames(allSNPSGenes.trt.mtd_new.2)[colnames(allSNPSGenes.trt.mtd_new.2)=="GO1"] <- "GO"
str(allSNPSGenes.trt.mtd_new.2)
allSNPSGenes.trt.mtd_new.3 <- allSNPSGenes.trt.mtd_new[,-c(10,12:18)][which(allSNPSGenes.trt.mtd_new$GO2 != "NA"),]
colnames(allSNPSGenes.trt.mtd_new.3)[colnames(allSNPSGenes.trt.mtd_new.3)=="GO2"] <- "GO"
str(allSNPSGenes.trt.mtd_new.3)
allSNPSGenes.trt.mtd_new.4 <- allSNPSGenes.trt.mtd_new[,-c(10:11,13:18)][which(allSNPSGenes.trt.mtd_new$GO3 != "NA"),]
colnames(allSNPSGenes.trt.mtd_new.4)[colnames(allSNPSGenes.trt.mtd_new.4)=="GO3"] <- "GO"
str(allSNPSGenes.trt.mtd_new.4)
allSNPSGenes.trt.mtd_new.5 <- allSNPSGenes.trt.mtd_new[,-c(10:12,14:18)][which(allSNPSGenes.trt.mtd_new$GO4 != "NA"),]
colnames(allSNPSGenes.trt.mtd_new.5)[colnames(allSNPSGenes.trt.mtd_new.5)=="GO4"] <- "GO"
str(allSNPSGenes.trt.mtd_new.5)
allSNPSGenes.trt.mtd_new.6 <- allSNPSGenes.trt.mtd_new[,-c(10:13,15:18)][which(allSNPSGenes.trt.mtd_new$GO5 != "NA"),]
colnames(allSNPSGenes.trt.mtd_new.6)[colnames(allSNPSGenes.trt.mtd_new.6)=="GO5"] <- "GO"
str(allSNPSGenes.trt.mtd_new.6)
allSNPSGenes.trt.mtd_new.7 <- allSNPSGenes.trt.mtd_new[,-c(10:14,16:18)][which(allSNPSGenes.trt.mtd_new$GO6 != "NA"),]
colnames(allSNPSGenes.trt.mtd_new.7)[colnames(allSNPSGenes.trt.mtd_new.7)=="GO6"] <- "GO"
str(allSNPSGenes.trt.mtd_new.7)
allSNPSGenes.trt.mtd_new.8 <- allSNPSGenes.trt.mtd_new[,-c(10:15,17:18)][which(allSNPSGenes.trt.mtd_new$GO7 != "NA"),]
colnames(allSNPSGenes.trt.mtd_new.8)[colnames(allSNPSGenes.trt.mtd_new.8)=="GO7"] <- "GO"
str(allSNPSGenes.trt.mtd_new.8)
allSNPSGenes.trt.mtd_new.9 <- allSNPSGenes.trt.mtd_new[,-c(10:16,18)][which(allSNPSGenes.trt.mtd_new$GO8 != "NA"),]
colnames(allSNPSGenes.trt.mtd_new.9)[colnames(allSNPSGenes.trt.mtd_new.9)=="GO8"] <- "GO"
str(allSNPSGenes.trt.mtd_new.9)
allSNPSGenes.trt.mtd_new.10 <- allSNPSGenes.trt.mtd_new[,-c(10:17)][which(allSNPSGenes.trt.mtd_new$GO9 != "NA"),]
colnames(allSNPSGenes.trt.mtd_new.10)[colnames(allSNPSGenes.trt.mtd_new.10)=="GO9"] <- "GO"
str(allSNPSGenes.trt.mtd_new.10)

library(lessR)
total <- Merge(allSNPSGenes.trt.mtd_new.2, allSNPSGenes.trt.mtd_new.3)
total <- Merge(total, allSNPSGenes.trt.mtd_new.4)
total <- Merge(total, allSNPSGenes.trt.mtd_new.5)
total <- Merge(total, allSNPSGenes.trt.mtd_new.6)
total <- Merge(total, allSNPSGenes.trt.mtd_new.7)
total <- Merge(total, allSNPSGenes.trt.mtd_new.8)
total <- Merge(total, allSNPSGenes.trt.mtd_new.9)
total <- Merge(total, allSNPSGenes.trt.mtd_new.10)
str(total)
total$GO
total$GO[total$GO ==""]<- NA
total <- total[,c(18:22,17,2:3,10:14,16)]
str(total)
total.1 <- total[!duplicated(total[,c("Trait.name","Method","Chromosome","Position","P.value",
                                                                                 "SNP","locusName","arabi.symbol","arabi.defline","rice.defline","GO")]),]
###
head(total.1)



###import the txt file from tair webs 
ATHgenes <- read.delim("/Users/yonglizhao/Documents/Genes function/ATH_GO_GOSLIM.txt",header=F,
                       stringsAsFactors = FALSE)
names(ATHgenes) <- c("locus.name","TAIR.accession",
                     "object.name","relationship.type","GO.term","GO","TAIR.Keyword.ID","Aspect","GOslim.term","Evidence.code",
                     "Evidence.description","Evidence","Reference","Annotator","Date.annotated")

head(ATHgenes)
ATHgenes <- unique(ATHgenes)

index <- match(total.1$GO, ATHgenes$GO)
ATHgenes_go <- ATHgenes[index,]

total$GO %in% ATHgenes_go$GO
#total.2 <- plyr::join(total, ATHgenes_go, by="GO", type = "left", match = "all")
rm(t_1) rm(total.1) rm(total.2) rm(t)
t_1 <- merge(total.1,ATHgenes_go,by.x = "GO", by.y="GO")
str(t_1)
t_2 <- t_1[,-c(15,17,20)]
str(t_2)
t_3 <- t_2[!duplicated(t_2[,c("Trait.name","Method","Chromosome","Position","P.value",
                                      "SNP","locusName","GO")]),]
str(t_3)
head(t_3)

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
write.csv(genefucntion,file="data/genefunctiontair.csv",row.names = F)
head(allSNPSGenes.trt.mtd)
str(allSNPSGenes.trt.mtd)



#library(stringi)
#ATHgenes$Function = stri_join(ATHgenes$V4,ATHgenes$V5,sep=" ")


###
head(allSNPSGenes.trt.mtd_new)

index <- match(as.factor(allSNPSGenes.trt.mtd_new[i]), ATHgenes$GO.ID)



###merge all the gene function





allSNPSGenes.trt.mtd.fun <- allSNPSGenes.trt.mtd.fun[!duplicated(allSNPSGenes.trt.mtd.fun[,c("Trait.name","Method","Chromosome","Position",
                                                                                             "P.value","SNP","locusName","Pfam",
                                                                                             "Best.hit.arabi.name","Function")]),]
str(allSNPSGenes.trt.mtd.fun)
allSNPSGenes.trt.mtd.fun <- allSNPSGenes.trt.mtd.fun[,c(11:12,10,13:15,1:7,18,27,19:22,24:26,8:9)]
head(allSNPSGenes.trt.mtd.fun)
str(allSNPSGenes.trt.mtd.fun)
###separate the file into small file depending the Traits. name
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/Alltraitsgenes")
for (i in levels(allSNPSGenes.trt.mtd.fun$Trait.name)){
  Trait <- allSNPSGenes.trt.mtd.fun[which(allSNPSGenes.trt.mtd.fun$Trait.name== i),]
  Trait <- Trait[order(Trait$Method,Trait$Chromosome,Trait$Position),]
  write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/Alltraitsgenes")
for (i in levels(allSNPSGenes.trt.mtd.fun$Trait.name)){
  Trait <- allSNPSGenes.trt.mtd.fun[which(allSNPSGenes.trt.mtd.fun$Trait.name== i),]
  Trait <- Trait[order(Trait$Method,Trait$Chromosome,Trait$Position),]
  write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
}
