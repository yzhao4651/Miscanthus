###this one used for search the genes around SNPs

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
flo.maf.t.7.14 <- read.csv("data/flo.T.U.7.14.csv")
flo.maf.t.7.14 <- read.csv("data/flo.T.U.7.14.csv")
flo.maf.t.7.27 <- read.csv("flo.7.27/flo7.27.unique.csv")
flo.maf.t.7.31 <- read.csv("flo.7.31/flo.7.31.unique.csv")
flo.maf.t.8.4 <- read.csv("flo.8.4/flo.8.4.unique.csv")
str(flo.maf.t.8.4)
##select the MAF < 0.01

Chr11.gff3.txt

flo.maf.t.8.4.MAF <- flo.maf.t.8.4[flo.maf.t.8.4$MAF >= 0.05,]

flo.maf.t.8.4.MAF.01 <- flo.maf.t.8.4[flo.maf.t.8.4$MAF < 0.05,]
str(flo.maf.t.7.14)
str(flo.maf.t.8.4)
library(rtracklayer)
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("Sbicolor_454_v3.1.1.gene.gff3",
                        format = "gff3",
                        dataSource = "Phytosome 12",
                        organism = "Sorghum bicolor")




mySNPs <- flo.maf.t.8.4
source("Function/Findneargenes.R")
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

gene_df.genes <- Findneargenes(flo.maf.t.8.4, 1e4)
str(gene_df.genes)
nlevels(gene_df.genes$Gene)
#seqinfo(txdb)
#To save the txdb database for later and avoid having to recreate it every time we use it, we can use saveDb() and, later, loadDb()
#BiocManager::install("AnnotationDbi")
#library(AnnotationDbi)
#saveDb(txdb, 'txdb.Sbicolor_454_v3.1.1')
#Once the database is stored to disk, it can be reloaded later, almost instantly.
#library(AnnotationDbi)
#txdb = loadDb(file = 'txdb.Sbicolor_454_v3.1.1')
#genes(txdb)
#seqinfo(txdb)
###install.packags library(IRanges)
#BiocManager::install("IRanges")
### fucntion to get the neargenes 
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
### this one for flowering traits
#full_data.flo <- read.csv("data/flo-106.116.122.106im.csv")
#str(full_data.flo)
#levels(full_data.flo$SNP)
#nlevels(full_data.flo$SNP)
##this function need input dataset and search_radius 
source("https://bioconductor.org/biocLite.R")
biocLite("Repitools")
library(Repitools)
columns(txdb)
mySNPs <- flo.maf.t.8.4
source("Function/Findneargenes.R")
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
      t <- annoGR2DF(mygenes) 
    saved_genes[[i]] <- t$chr
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        chr = unlist(saved_genes))
  return(gene_df)
}


gene_df.chr <- Findneargenes(flo.maf.t.8.4, 1e4)
str(gene_df.chr)

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

    t <- annoGR2DF(mygenes)
    saved_genes[[i]] <- t$start
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        start = unlist(saved_genes))
  return(gene_df)
}

gene_df.start <- Findneargenes(flo.maf.t.8.4, 1e4)
str(gene_df.start)


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
    
    
    t <- annoGR2DF(mygenes)
    saved_genes[[i]] <- t$end
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        end = unlist(saved_genes))
  return(gene_df)
}

gene_df.end <- Findneargenes(flo.maf.t.8.4, 1e4)
str(gene_df.end)
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
    
    
    t <- annoGR2DF(mygenes)
    saved_genes[[i]] <- t$width
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        width = unlist(saved_genes))
  return(gene_df)
}

gene_df.width <- Findneargenes(flo.maf.t.8.4, 1e4)
str(gene_df.width)

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
    
    
    t <- annoGR2DF(mygenes)
    saved_genes[[i]] <- t$strand
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        strand = unlist(saved_genes))
  return(gene_df)
}

gene_df.strand <- Findneargenes(flo.maf.t.8.4, 1e4)
str(gene_df.strand)

###combine all of them together 
gene_df8.4 <- do.call("cbind",list(gene_df.genes,gene_df.chr, gene_df.start, gene_df.end, gene_df.width, gene_df.strand))
str(gene_df8.4)

gene_df8.4.a <- gene_df8.4[,c(1,2,4,6,8,10,12)]
gene_df8.4.a <- droplevels(gene_df8.4.a)

nlevels(gene_df8.4.a$SNP) ##1175
nlevels(gene_df8.4.a$Gene)###4603

write.csv(gene_df8.4.a,file="flo.8.4/AHG.flo.8.4.T.hitandnohit.csv")

##merge all of data withunique

flo.8.4.T <- merge(flo.maf.t.8.4, gene_df8.4.a,by="SNP",all.x=TRUE)
flo.8.4.T <- droplevels(flo.8.4.T)
nlevels(flo.maf.t.8.4$SNP)##1857
nlevels(flo.8.4.T$Gene)##4603
nlevels(flo.8.4.T$SNP)###1857
###get all of the no any hit SNP
flo.8.4.T.nohhit <- flo.8.4.T[is.na(flo.8.4.T$Gene),]
flo.8.4.T.nohhit <- droplevels(flo.8.4.T.nohhit)
nlevels(flo.8.4.T.nohhit$SNP)###682
write.csv(flo.8.4.T.nohhit,file="flo.8.4/AHG.flo.8.4.T.nohit.csv")

###all of them with hit gene
flo.8.4.T.hit <- flo.8.4.T[!(is.na(flo.8.4.T$Gene)),]
str(flo.8.4.T.hit)
flo.8.4.T.hit$dis_start.positon <- flo.8.4.T.hit$start - flo.8.4.T.hit$Position
flo.8.4.T.hit$dis_positon.end <- flo.8.4.T.hit$Position- flo.8.4.T.hit$end
flo.8.4.T.hit <- droplevels(flo.8.4.T.hit)
nlevels(flo.8.4.T.hit$SNP) ###1175
nlevels(flo.8.4.T.hit$Gene)###4603

write.csv(flo.8.4.T.hit,file="flo.8.4/AHG.flo.8.4.T.hit.csv")
###and then get gene anotation 
annotations <- read.delim("Sbicolor_454_v3_1_1_annotation_info.txt",
                          stringsAsFactors = FALSE)
head(annotations)
str(annotations)
myrows <- match(flo.8.4.T.hit$Gene, annotations$transcriptName)
allhitsgenes.flo <- annotations[myrows,]
head(allhitsgenes.flo)
###combine two dataset in order to get the SNP with gene anotation in one file. 
AHG.flo.1 <- cbind(flo.8.4.T.hit,allhitsgenes.flo)
colnames(AHG.flo.1)[colnames(AHG.flo.1)=="Gene"] <- "Gene"

##seperate the Gene with number and then delete the duplicated genes
AHG.flo.2 <- separate(AHG.flo.1, Gene, c("Gene1","Gene2","Gene.num"), "\\.") 
AHG.flo.3 <- unite(AHG.flo.2,"SNP.Gene",SNP,Gene1,Gene2,sep="/")
AHG.flo.4<- AHG.flo.3[!(duplicated(AHG.flo.3$SNP.Gene)),]
AHG.flo.5 <- separate(AHG.flo.4, SNP.Gene, c("SNP","Gene1","Gene2"), "/")
AHG.flo.5 <- unite(AHG.flo.5, Gene, c("Gene1","Gene2"), sep=".")
AHG.flo.5 <- droplevels(AHG.flo.5)
nlevels(as.factor(AHG.flo.5$SNP)) ###1175
nlevels(as.factor(AHG.flo.5$Gene))###3315

write.csv(AHG.flo.5,file="flo.8.4/AHG.flo.8.4.anota.csv")
###select the Genes without studying: which means the X.pac is NA 
flo.hit.gene <- read.csv("flo.8.4/AHG.flo.8.4.anota.csv")
str(flo.hit.gene)
nlevels(as.factor(flo.hit.gene$Gene))

###this one works even this is hit with sobic genes, but many of them has no any study 80 has no any function found 
flo.NOX.pacId <- flo.hit.gene[is.na(flo.hit.gene$X.pacId),]
flo.NOX.pacId <- droplevels(flo.NOX.pacId)
nlevels(flo.NOX.pacId$SNP)###148
nlevels(flo.NOX.pacId$Gene)###415
write.csv(flo.NOX.pacId,file="flo.8.4/AHG.flo.NOX.pacId.csv")

## with some hit 
flo.X.pacId <- flo.hit.gene[!is.na(flo.hit.gene$X.pacId),]
flo.X.pacId <- droplevels(flo.X.pacId)
nlevels(flo.X.pacId$SNP)###1028
nlevels(flo.X.pacId$Gene)###2900
write.csv(flo.X.pacId,file="flo.8.4/AHG.flo.X.pacId.csv")

###check how many hit with ara 
flo.NO.ara <- flo.X.pacId[flo.X.pacId$Best.hit.arabi.name=="",]
flo.NO.ara <- droplevels(flo.NO.ara)
names(flo.NO.ara)
nlevels(flo.NO.ara$SNP)###425
nlevels(flo.NO.ara$Gene)###527
write.csv(flo.NO.ara,file="flo.8.4/AHG.flo.X.pacId.NO.ara.csv")


###check how many hit with ara 
flo.X.pacId.ara <- flo.X.pacId[!(flo.X.pacId$Best.hit.arabi.name==""),]
flo.X.pacId.ara <- droplevels(flo.X.pacId.ara)
names(flo.X.pacId.ara)
nlevels(flo.X.pacId.ara$SNP)###992
nlevels(flo.X.pacId.ara$Gene)###2373
write.csv(flo.X.pacId.ara,file="flo.8.4/AHG.flo.X.pacId.ara.csv")
###check the GO information
### and then using the data file to find the SNPs hit GO term 
flo.X.pacId.GO <- flo.X.pacId[!(flo.X.pacId$GO==""),]
flo.X.pacId.GO <- droplevels(flo.X.pacId.GO)
names(flo.X.pacId.GO)
nlevels(flo.X.pacId.GO$SNP)###837
nlevels(flo.X.pacId.GO$Gene)###1429
write.csv(flo.X.pacId.GO,file="flo.8.4/AHG.flo.X.pacI.GO.csv")
#### rice 
flo.X.pacId.rice <- flo.X.pacId[!(flo.X.pacId$Best.hit.rice.name==""),]
flo.X.pacId.rice <- droplevels(flo.X.pacId.rice)
names(flo.X.pacId.rice)
nlevels(flo.X.pacId.rice$SNP)###1015
nlevels(flo.X.pacId.rice$Gene)###2549
write.csv(flo.X.pacId.rice,file="flo.8.4/AHG.flo.X.pacI.rice.csv")


###separate GO 
allSNPSGenes.trt.mtd <- flo.X.pacId.GO
source("Function/get GO.R")
flo.X.pacId.GO.1 <- GetGO(allSNPSGenes.trt.mtd)
write.csv(flo.X.pacId.GO.1,file="flo.8.4/AHG.flo.X.pacId.GO.1.csv")

###trying to check the GO term with ATH 
ATHgenes <- read.delim("/Users/yonglizhao/Documents/Genes function/ATH_GO_GOSLIM.txt",header=F,
                       stringsAsFactors = FALSE)
names(ATHgenes) <- c("locus.name","TAIR.accession",
                     "object.name","relationship.type","GO.term","GO","TAIR.Keyword.ID","Aspect","GOslim.term","Evidence.code",
                     "Evidence.description","Evidence","Reference","Annotator","Date.annotated")

head(ATHgenes)
head(ATHgenes)
ATHgenes <- unique(ATHgenes)

##using the hit ara. name to get the GOterm
index <- match(flo.X.pacId.GO.1$GO, ATHgenes$GO)
ATHgenes_go <- ATHgenes[index,]
###write out the data set with hit to GO
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
write.csv(ATHgenes_go,file="flo.8.4/ATHgenes_go.csv")

str(ATHgenes_go)
levels(as.factor(ATHgenes_go$GOslim.term))


###this one for searching the flowering time genes 
###sepate the gene with .1 using function 
library(tidyr)
AHG.flo.2 <- separate(AHG.flo.1, Gene, c("Gene1","Gene2","Gene.num"), "\\.") 
### and then remove the same genes and same SNP
AHG.flo.3 <- AHG.flo.2[!duplicated(AHG.flo.2[,c("SNP","Gene2")]),]

AHG.flo.3 <- unite(AHG.flo.3,"Gene",Gene1,Gene2,sep=".")
str(AHG.flo.3)
AHG.flo.4 <- separate(AHG.flo.3, Best.hit.arabi.name, c("Best.hit.arabi.name","Best.hit.arabi.name.num"), "\\.") 
AHG.flo.5 <- separate(AHG.flo.4, Best.hit.rice.name, c("Best.hit.rice.name","Best.hit.rice.name.num"), "\\.") 
str(AHG.flo.5)
###import the candidate gene
Ara.gene.306 <- read.csv("data/ARA.candita.genes.csv")
str(Ara.gene.306)

##hit with flowering time genes with ara

AHG.flo.5.ara.hit <- AHG.flo.5[AHG.flo.5$Best.hit.arabi.name %in% Ara.gene.306$ARA.can.gene,]
AHG.flo.5.ara.hit$hit.gene <- paste("arabidopsis.306")


###import the candidategene list 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
candidategene1 <- read.csv("data/flo.time.Sorghum.Rice.ARA.csv")
str(candidategene1)
###hit with Sorgu

AHG.flo.5$Gene %in% candidategene$Homologous.Gene

AHG.flo.5.Sor.hit <- AHG.flo.5[AHG.flo.5$Gene %in% candidategene1$Sorghum,]

AHG.flo.5.Sor.hit$hit.gene <- paste("sorghum")
write.csv(AHG.flo.5.Sor.hit,file="data/AHG.flo.5.Sor.hit.csv")
###with rice 

candidategene2 <- read.csv("data/rice.flo.t.csv")
str(candidategene2)
AHG.flo.5.Rice.hit.2 <- AHG.flo.5[AHG.flo.5$Best.hit.rice.name %in% candidategene2$Rice,]
head(AHG.flo.5.Rice.hit.2)

AHG.flo.5.Rice.hit.2$hit.gene <- paste("oryza(rice)")
write.csv(AHG.flo.5.Rice.hit.2,file="data/AHG.flo.5.Rice.hit.csv")
###with ara from candidata list 

AHG.flo.5.ara.hit.2 <- AHG.flo.5[AHG.flo.5$Best.hit.arabi.name %in% candidategene1$ARA,]
AHG.flo.5.ara.hit.2$hit.gene <- paste("arabidopsis.alternative")
write.csv(AHG.flo.5.ara.hit.2, file="flo.8.4/AHG.flo.5.ara.hit.2.csv")
##combine all of the gens together 

AHG.flo.hit.T <- do.call("rbind",list(AHG.flo.5.ara.hit,AHG.flo.5.Sor.hit,AHG.flo.5.Rice.hit.2,AHG.flo.5.ara.hit.2))

write.csv(AHG.flo.hit.T, file="flo.8.4/AHG.flo.hit.T.8.8.2.csv")



